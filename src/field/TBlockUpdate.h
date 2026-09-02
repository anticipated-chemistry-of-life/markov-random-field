//
// The traversal that draws the field and both tree fields, one leaf pair at a time.
//
// A leaf pair carries three variables. They are drawn together, from all eight combinations at
// once, because single-variable draws leave the triple metastable. ADR-0005 carries the argument,
// and field/TFieldMath.h carries the arithmetic.
//
// This file is the loop and nothing else. The probabilities and the counter move come from the
// kernel. The cells come and go through windows. Everything else comes from a model the loop asks
// one row at a time (field/TBlockModel.h). What is left is a loop a test can state properties
// against: every leaf pair visited once, the right neighbours read, the writes landing, and one
// chain at any thread count.
//
// A thread takes a species leaf. A leaf pair's Markov blanket holds no cell of another leaf pair.
// It holds the two tree parents and the two data terms. So the rows are conditionally independent.
// A cell's uniform is hashed from its position (ADR-0007), so the thread that reaches a cell does
// not decide what it gets.
//

#pragma once

#include "cli.h"
#include "constants.h"
#include "coretools/Main/TError.h"
#include "coretools/Types/probability.h"
#include "field/TFieldMath.h"
#include "omp.h"
#include "random/TCellUniforms.h"
#include "tree/TPhylogeny.h"
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <vector>

namespace block_update {

// The update names the two trees. A thread takes a species leaf and walks the molecule leaves of
// its row. A third tree would be a third window and a third factor, so this says so here rather
// than looping over a dimension count that cannot change.
static_assert(NUMBER_OF_TREES == 2,
              "The block update is written for one species tree and one molecule tree.");

/// The factors at one leaf pair that the update's own windows do not hold: what each tree's
/// process says about its own tree field cell, and what each data source says about the field cell.
struct TLeafPairFactors {
	/// P(Z_s = 1 | the species parent's state), from the species tree's transition grid.
	coretools::Probability prob_z_s_is_one;
	/// P(Z_m = 1 | the molecule parent's state), from the molecule tree's transition grid.
	coretools::Probability prob_z_m_is_one;
	/// {P(L | Y = 0), P(L | Y = 1)} for this cell.
	std::array<coretools::Probability, 2> lotus;
	/// {P(D | Y = 0), P(D | Y = 1)} for this cell.
	std::array<coretools::Probability, 2> simple_error;
};

/// The model at one row of the field: one species leaf against every molecule leaf.
///
/// `factors` is asked once per leaf pair, and is handed the two parent states the traversal read
/// through its own windows. `record` is told what the leaf pair was given, so each data source can
/// carry its own bookkeeping forward. The factors come back with it, because a source that scored
/// the cell scored both field states and only now knows which one to keep.
template<typename T>
concept BlockRow = requires(T &row, size_t molecule_leaf, bool species_parent, bool molecule_parent,
                            const TLeafPairFactors &factors,
                            const field_math::TBlockStates &drawn) {
	{ row.factors(molecule_leaf, species_parent, molecule_parent) } -> std::same_as<TLeafPairFactors>;
	{ row.record(molecule_leaf, factors, drawn) } -> std::same_as<void>;
};

/// Everything the traversal reads that is not a state in one of its own windows. It hands out one
/// row per species leaf, which is a thread's whole share of the update.
template<typename T>
concept BlockModel = requires(T &model, size_t species_leaf) {
	{ model.open_row(species_leaf) } -> BlockRow;
};

/// What one thread's share of a block update added up to.
///
/// The counters are the link's sufficient statistic, `n(bucket, field state)`. A block update
/// visits every leaf pair exactly once, so what one update tallies is the whole configuration
/// rather than a delta. A thread can therefore hold its own tally without a count going negative.
/// The caller merges them after the parallel region.
struct TThreadTally {
	field_math::TLinkCounters counters;
	/// The two tree factors at the drawn states, summed over this thread's cells. The link's own
	/// term is not here: it comes from the merged counters, in one call, for the whole field.
	double log_density = 0.0;
};

/// The two tree factors at the states one draw assigned. Log of `P(Z_s = z_s | parent)` plus log of
/// `P(Z_m = z_m | parent)`.
[[nodiscard]] inline double log_tree_factors(const field_math::TBlockStates &drawn,
                                             const TLeafPairFactors &factors) {
	const double p_s = factors.prob_z_s_is_one.get();
	const double p_m = factors.prob_z_m_is_one.get();
	return std::log(drawn.z_s ? p_s : 1.0 - p_s) + std::log(drawn.z_m ? p_m : 1.0 - p_m);
}

/// Draws every leaf pair of one species leaf's row.
///
/// The row is a thread's whole share of the update. It opens four windows on state: the field's
/// row, the species tree field's row, the species parent's row, and the molecule node state's row
/// over every molecule node. The last of those covers both that tree field and every molecule
/// parent the row reads, because a molecule parent is a column of that one row.
///
/// The three written windows hand out the cells they could not write in place. That is the only
/// exit a window inside a parallel region may take (ADR-0006).
template<field_math::LinkPolicy Policy, typename Field, typename NodeState, BlockModel Model>
void update_row(size_t species_leaf, Field &Y, NodeState &Z_species, NodeState &Z_molecule,
                const TPhylogeny &species, const TPhylogeny &molecule,
                const field_math::TErrorProbability &omega, Model &model,
                const TCellUniforms &uniforms, TThreadTally &tally,
                std::vector<size_t> &field_inserts, std::vector<size_t> &species_inserts,
                std::vector<size_t> &molecule_inserts) {
	const size_t n_molecule_leaves = molecule.n_leaves();
	const IndexArray row_start{species_leaf, 0};

	// The field's row, and the species tree field's row beside it. A leaf pair sits at the same
	// (row, column) in both (ADR-0005), so one subscript addresses both.
	auto field_row   = Y.open_window(row_start, n_molecule_leaves, /*stride=*/1);
	auto species_row = Z_species.open_window(row_start, n_molecule_leaves, /*stride=*/1);

	// The species parent's row. A leaf is never a root, so its parent is an internal node, and an
	// internal node is never a leaf -- so no thread writes the row another thread reads here.
	auto species_parent_row = Z_species.open_window(IndexArray{species.parent_of(species_leaf), 0},
	                                                n_molecule_leaves, /*stride=*/1);

	// The molecule node state's row, over every molecule node. A cell's molecule parent is a
	// column of this one row, so one window covers the tree field and every parent read.
	auto molecule_row = Z_molecule.open_window(row_start, molecule.n_nodes(), /*stride=*/1);

	auto row = model.open_row(species_leaf);

	for (size_t molecule_leaf = 0; molecule_leaf < n_molecule_leaves; ++molecule_leaf) {
		const bool species_parent  = species_parent_row.is_one(molecule_leaf);
		const bool molecule_parent = molecule_row.is_one(molecule.parent_of(molecule_leaf));
		const TLeafPairFactors factors = row.factors(molecule_leaf, species_parent, molecule_parent);

		const field_math::TBlockStates current{.y   = field_row.is_one(molecule_leaf),
		                                       .z_s = species_row.is_one(molecule_leaf),
		                                       .z_m = molecule_row.is_one(molecule_leaf)};

		// The cell names the uniform it draws, so the thread that happens to reach this cell does
		// not decide what it gets.
		const auto draw = field_math::draw_block<Policy>(
		    factors.prob_z_s_is_one, factors.prob_z_m_is_one, omega, factors.lotus,
		    factors.simple_error, current, uniforms.at(field_row.linear_index(molecule_leaf)));

		// Write each cell back through the window it was read from. A window that holds the cell
		// writes it in place, and one that does not buffers the insert for the caller.
		if (current.y != draw.drawn.y) { field_row.set_state(molecule_leaf, draw.drawn.y); }
		if (current.z_s != draw.drawn.z_s) { species_row.set_state(molecule_leaf, draw.drawn.z_s); }
		if (current.z_m != draw.drawn.z_m) {
			molecule_row.set_state(molecule_leaf, draw.drawn.z_m);
		}

		tally.counters.add(draw.to.bucket, draw.to.y);
		tally.log_density += log_tree_factors(draw.drawn, factors);
		row.record(molecule_leaf, factors, draw.drawn);
	}

	field_inserts    = field_row.take_buffered_inserts();
	species_inserts  = species_row.take_buffered_inserts();
	molecule_inserts = molecule_row.take_buffered_inserts();

	// The parent's row is read and never written, so it buffers nothing. It is ended here all the
	// same, so that no window of the four reaches its storage on the way out.
	const std::vector<size_t> parent_inserts = species_parent_row.take_buffered_inserts();
	DEBUG_ASSERT(parent_inserts.empty());
}

/// One block update: every leaf pair, one species leaf per thread.
///
/// `tallies` is one entry per thread, and the caller merges them. It is not cleared here: a caller
/// that hands over used tallies gets their counts added to, which is what makes the merge its
/// business rather than this one's.
template<field_math::LinkPolicy Policy, typename Field, typename NodeState, BlockModel Model>
void run(Field &Y, NodeState &Z_species, NodeState &Z_molecule, const TPhylogeny &species,
         const TPhylogeny &molecule, const field_math::TErrorProbability &omega, Model &model,
         const TCellUniforms &uniforms, std::vector<TThreadTally> &tallies) {
	// One list per species leaf, filled by the window that leaf's row is written through. No two
	// leaves share a list, so the lists need neither a thread index nor a lock, and nothing has to
	// be drained inside the region.
	// Not const: OpenMP made a const variable predetermined shared before version 4.0 and does not
	// now, so a `default(none)` clause that names one is right under some compilers and wrong under
	// others.
	size_t n_species_leaves = species.n_leaves();
	std::vector<std::vector<size_t>> field_inserts(n_species_leaves);
	std::vector<std::vector<size_t>> species_inserts(n_species_leaves);
	std::vector<std::vector<size_t>> molecule_inserts(n_species_leaves);

	// The species tree's leaf count is the widest this team can run. Splitting the other way would
	// give a thread a column, and a column of a sparse container is what no window may insert into
	// (ADR-0006).
#pragma omp parallel for num_threads(ProgramOptions::NUMBER_OF_THREADS)                            \
    schedule(static) default(none)                                                                 \
    shared(Y, Z_species, Z_molecule, species, molecule, omega, model, uniforms, tallies,           \
               field_inserts, species_inserts, molecule_inserts, n_species_leaves)
	for (size_t species_leaf = 0; species_leaf < n_species_leaves; ++species_leaf) {
		update_row<Policy>(species_leaf, Y, Z_species, Z_molecule, species, molecule, omega, model,
		                   uniforms, tallies[static_cast<size_t>(omp_get_thread_num())],
		                   field_inserts[species_leaf], species_inserts[species_leaf],
		                   molecule_inserts[species_leaf]);
	}

	// The inserts every row deferred, in one batch per container. A dense container is handed
	// empty lists and does nothing with them.
	Y.insert_in_Y(field_inserts);
	Z_species.insert_in_Z(species_inserts);
	Z_molecule.insert_in_Z(molecule_inserts);
}

} // namespace block_update
