//
// The traversal that draws the field and both tree fields, one leaf pair at a time.
//
// A leaf pair carries three variables. They are drawn together, from all eight combinations at
// once, because single-variable draws leave the triple metastable. ADR-0005 carries the argument,
// and field/TFieldMath.h carries the arithmetic.
//
// This file is the loop and nothing else. The probabilities and the counter move come from the
// kernel. The cells are read and written one at a time, straight from the storages
// (storages/cell_write.h). Everything else comes from a model the loop asks one leaf pair at a
// time (field/TBlockModel.h). What is left is a loop a test can state properties against: every
// leaf pair visited once, the right neighbours read, the writes landing, and one chain at any
// thread count.
//
// A thread takes a species leaf. A leaf pair's Markov blanket holds no cell of another leaf pair.
// It holds the two tree parents and the two data terms. So the rows are conditionally independent.
// A cell's uniform is hashed from its position (ADR-0007), so the thread that reaches a cell does
// not decide what it gets.
//

#pragma once

#include "cli.h"
#include "constants.h"
#include "coretools/Types/probability.h"
#include "field/TFieldMath.h"
#include "omp.h"
#include "random/TCellUniforms.h"
#include "storages/cell_write.h"
#include "tree/TPhylogeny.h"
#include <array>
#include <concepts>
#include <cstddef>
#include <vector>

namespace block_update {

// The update names the two trees. A thread takes a species leaf and walks the molecule leaves of
// its row. A third tree would be a third storage and a third factor, so this says so here rather
// than looping over a dimension count that cannot change.
static_assert(NUMBER_OF_TREES == 2,
              "The block update is written for one species tree and one molecule tree.");

/// The factors at one leaf pair that the update does not read from its storages: what each tree's
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

/// Everything the traversal reads that is not a state in one of its three storages.
///
/// `factors` is asked once per leaf pair, and is handed the two parent states the traversal read
/// from the node states. `record` is told what the leaf pair was given, so each data source can
/// carry its own bookkeeping forward. The factors come back with it, because a source that scored
/// the cell scored both field states and only now knows which one to keep.
///
/// The model answers one leaf pair at a time. Threads ask it at once, so `factors` reads and
/// writes nothing a second thread also touches.
template<typename T>
concept BlockModel = requires(T &model, size_t species_leaf, size_t molecule_leaf,
                              bool species_parent, bool molecule_parent,
                              const TLeafPairFactors &factors,
                              const field_math::TBlockStates &drawn) {
	{
		model.factors(species_leaf, molecule_leaf, species_parent, molecule_parent)
	} -> std::same_as<TLeafPairFactors>;
	{ model.record(species_leaf, molecule_leaf, factors, drawn) } -> std::same_as<void>;
};

/// What one thread's share of a block update added up to.
///
/// The counters are the link's sufficient statistic, `n(bucket, field state)`. A block update
/// visits every leaf pair exactly once, so what one update tallies is the whole configuration
/// rather than a delta. A thread can therefore hold its own tally without a count going negative.
/// The caller merges them after the parallel region.
/// The traversal keeps no running density. What one draw was worth is not what the configuration
/// is worth: the leaf states move again in the same iteration, and so do the parameters. The joint
/// density is asked of the configuration the iteration leaves behind instead
/// (tree/node_state_density.h).
struct TThreadTally {
	field_math::TLinkCounters counters;
};

/// Draws every leaf pair of one species leaf's row.
///
/// The row is a thread's whole share of the update. It reads five cells per leaf pair -- the three
/// the draw moves, and the two tree parents -- and writes back the ones the draw changed. A leaf
/// pair sits at the same multidimensional index in all three containers (ADR-0005), so one
/// subscript names its cell in each; the linear indices differ, because the container shapes do.
///
/// A write reaches the storage in place, or comes back as a deferred insert. That is the only exit
/// a write inside a parallel region may take (ADR-0006).
template<field_math::LinkPolicy Policy, typename Field, typename NodeState, BlockModel Model>
void update_row(size_t species_leaf, Field &Y, NodeState &Z_species, NodeState &Z_molecule,
                const TPhylogeny &species, const TPhylogeny &molecule,
                const field_math::TErrorProbability &omega, Model &model,
                const TCellUniforms &uniforms, TThreadTally &tally,
                std::vector<size_t> &field_inserts, std::vector<size_t> &species_inserts,
                std::vector<size_t> &molecule_inserts) {
	const size_t n_molecule_leaves = molecule.n_leaves();

	// A leaf is never a root, so its parent is an internal node, and an internal node is never a
	// leaf -- so no thread writes the cells this row reads there.
	const size_t species_parent_node = species.parent_of(species_leaf);

	for (size_t molecule_leaf = 0; molecule_leaf < n_molecule_leaves; ++molecule_leaf) {
		const IndexArray leaf_pair{species_leaf, molecule_leaf};
		const size_t field_cell    = Y.get_linear_index_in_container_space(leaf_pair);
		const size_t species_cell  = Z_species.get_linear_index_in_container_space(leaf_pair);
		const size_t molecule_cell = Z_molecule.get_linear_index_in_container_space(leaf_pair);

		// The species parent is the cell above this one in the species node state. The molecule
		// parent is the cell beside it in the molecule node state.
		const bool species_parent = Z_species.is_one(Z_species.get_linear_index_in_container_space(
		    IndexArray{species_parent_node, molecule_leaf}));
		const bool molecule_parent =
		    Z_molecule.is_one(Z_molecule.get_linear_index_in_container_space(
		        IndexArray{species_leaf, molecule.parent_of(molecule_leaf)}));

		const TLeafPairFactors factors =
		    model.factors(species_leaf, molecule_leaf, species_parent, molecule_parent);

		const field_math::TBlockStates current{.y   = Y.is_one(field_cell),
		                                       .z_s = Z_species.is_one(species_cell),
		                                       .z_m = Z_molecule.is_one(molecule_cell)};

		// The cell names the uniform it draws, so the thread that happens to reach this cell does
		// not decide what it gets.
		const auto draw = field_math::draw_block<Policy>(
		    factors.prob_z_s_is_one, factors.prob_z_m_is_one, omega, factors.lotus,
		    factors.simple_error, current, uniforms.at(field_cell));

		// Write each cell back where it was read from. A storage that holds the cell takes the
		// write in place, and one that does not hands the cell back for a later insert.
		if (current.y != draw.drawn.y) {
			write_or_defer(Y, field_cell, draw.drawn.y, field_inserts);
		}
		if (current.z_s != draw.drawn.z_s) {
			write_or_defer(Z_species, species_cell, draw.drawn.z_s, species_inserts);
		}
		if (current.z_m != draw.drawn.z_m) {
			write_or_defer(Z_molecule, molecule_cell, draw.drawn.z_m, molecule_inserts);
		}

		tally.counters.add(draw.to.bucket, draw.to.y);
		model.record(species_leaf, molecule_leaf, factors, draw.drawn);
	}
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
	// One list per species leaf, filled by the writes of that leaf's row. No two leaves share a
	// list, so the lists need neither a thread index nor a lock, and nothing has to be drained
	// inside the region.
	// Not const: OpenMP made a const variable predetermined shared before version 4.0 and does not
	// now, so a `default(none)` clause that names one is right under some compilers and wrong under
	// others.
	size_t n_species_leaves = species.n_leaves();
	std::vector<std::vector<size_t>> field_inserts(n_species_leaves);
	std::vector<std::vector<size_t>> species_inserts(n_species_leaves);
	std::vector<std::vector<size_t>> molecule_inserts(n_species_leaves);

	// The species tree's leaf count is the widest this team can run.
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
