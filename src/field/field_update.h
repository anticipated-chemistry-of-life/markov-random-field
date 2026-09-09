//
// The field's own pass over its cells.
//
// The pass visits every cell of the field's container space. It draws each cell from the two tree
// field cells at that leaf pair, and from every data source that observes the field. It retallies
// the six link counters as it goes, so what one pass adds up to is a whole configuration and not
// a delta. Clearing the tallies between passes is the caller's, as the merge is.
//
// The tree fields are not the field's to draw. Each is drawn by its own tree, as the leaf block of
// that tree's node state, so this pass reads them and writes neither.
//
// This file is the loop and nothing else. The probability comes from the link
// (field/TFieldMath.h), the draw from random/two_state_draw.h, and the write from the storage
// (storages/cell_handle.h). Everything else comes from a model the loop asks one cell at a time.
// What is left is a loop a test can state properties against: every cell visited once, the right
// tree field cells read, the writes landing, one chain at any thread count, and a tally that
// matches a naive recount.
//
// A cell's Markov blanket holds no other cell of the field. It holds the two tree field cells and
// the data terms. So the cells are conditionally independent and the pass is one flat parallel
// loop. A cell's uniform is hashed from its position (ADR-0007), so the thread that reaches a cell
// does not decide what it gets. The stream is the caller's to name, and the field draws from the
// field stream at the cell's own linear index, one uniform per cell.
//

#pragma once

#include "cli.h"
#include "constants.h"
#include "coretools/Main/TError.h"
#include "coretools/Types/probability.h"
#include "field/TFieldMath.h"
#include "omp.h"
#include "random/TCellUniforms.h"
#include "random/two_state_draw.h"
#include "storages/cell_handle.h"
#include "storages/storage_concepts.h"
#include <array>
#include <concepts>
#include <cstddef>
#include <vector>

namespace field_update {

// The field's container space has one dimension per tree, and the link reads one tree field cell
// per tree. A third tree would be a third storage and a third state in the bucket, so this says so
// here rather than looping over a dimension count that cannot change.
static_assert(NUMBER_OF_TREES == 2,
              "The field update is written for one species tree and one molecule tree.");

/// What the data says about one field cell, which the pass does not read from its storages.
///
/// The neutral value is 1.0, which adds log(1) = 0. A build that compiled a source out leaves that
/// source's pair at it, so the loop draws from the link alone.
struct TCellFactors {
	/// {P(L | Y = 0), P(L | Y = 1)} for this cell.
	std::array<coretools::Probability, 2> lotus{coretools::P(1.0), coretools::P(1.0)};
	/// {P(D | Y = 0), P(D | Y = 1)} for this cell.
	std::array<coretools::Probability, 2> simple_error{coretools::P(1.0), coretools::P(1.0)};
};

/// Everything the pass reads that is not a state in one of its storages.
///
/// `factors` is asked once per cell. `record` is told what that cell was given, so each data source
/// can carry its own likelihood bookkeeping forward. The factors come back with it, because a
/// source that scored the cell scored both field states and only now knows which one to keep.
///
/// A cell is named by its multidimensional index, which is the leaf pair the data is indexed by.
///
/// The model answers one cell at a time. Threads ask it at once, so `factors` reads and writes
/// nothing a second thread also touches.
template<typename T>
concept FieldModel =
    requires(T &model, const IndexArray &cell, const TCellFactors &factors, bool drawn) {
	    { model.factors(cell) } -> std::same_as<TCellFactors>;
	    { model.record(cell, factors, drawn) } -> std::same_as<void>;
    };

/// Draws one field cell, and adds it to the tally.
///
/// It reads three cells -- the field cell, and the two tree field cells at that leaf pair -- and
/// writes back the one the draw moved. A leaf pair sits at the same multidimensional index in all
/// three containers (ADR-0004), so one subscript names its cell in each; the linear indices differ,
/// because the container shapes do.
///
/// The write reaches the storage in place, or comes back as a deferred insert. That is the only
/// exit a write inside a parallel region may take (ADR-0006).
template<field_math::LinkPolicy Policy, BinaryStorage Field, BinaryStorage NodeState,
         FieldModel Model, CellUniforms Uniforms>
void update_cell(size_t field_cell, Field &Y, const NodeState &Z_species,
                 const NodeState &Z_molecule, const field_math::TErrorProbability &omega,
                 Model &model, const Uniforms &uniforms, field_math::TLinkCounters &tally,
                 std::vector<size_t> &field_inserts) {
	const IndexArray leaf_pair = Y.get_multi_dimensional_index(field_cell);

	const bool z_s = Z_species.is_one(Z_species.get_linear_index_in_container_space(leaf_pair));
	const bool z_m = Z_molecule.is_one(Z_molecule.get_linear_index_in_container_space(leaf_pair));

	const TCellFactors factors                      = model.factors(leaf_pair);
	const coretools::Probability probability_of_one = field_math::prob_field_cell_is_one<Policy>(
	    z_s, z_m, omega, factors.lotus, factors.simple_error);

	// The cell names the uniform it draws, so the thread that happens to reach this cell does not
	// decide what it gets.
	const bool drawn = two_state_draw::sample(probability_of_one, uniforms.at(field_cell));

	// A storage that holds the cell takes the write in place, and one that does not hands the cell
	// back for a later insert.
	const auto handle = Y.locate(field_cell);
	if (handle.is_one != drawn) { write_or_defer(handle, drawn, field_inserts); }

	// The counters describe the configuration the pass leaves, so the cell counts in the bucket it
	// is in now.
	tally.add(Policy::bucket(z_s, z_m), drawn);
	model.record(leaf_pair, factors, drawn);
}

/// One field update: every cell of the field's container space, split over the threads.
///
/// `tallies` is one entry per thread, and the caller merges them. It is not cleared here: a caller
/// that hands over used tallies gets their counts added to, which is what makes the merge its
/// business rather than this one's.
///
/// A pass visits every cell exactly once, so what one thread tallies is a share of a whole
/// configuration rather than a delta. A thread can therefore hold its own tally without a count
/// going negative.
///
/// The pass keeps no running density. What one draw was worth is not what the configuration is
/// worth: the parameters move again in the same iteration. The joint density is asked of the
/// configuration the iteration leaves behind instead (field/joint_density.h).
template<field_math::LinkPolicy Policy, BinaryStorage Field, BinaryStorage NodeState,
         FieldModel Model, CellUniforms Uniforms>
void run(Field &Y, const NodeState &Z_species, const NodeState &Z_molecule,
         const field_math::TErrorProbability &omega, Model &model, const Uniforms &uniforms,
         std::vector<field_math::TLinkCounters> &tallies) {
	// One list per thread, filled by the writes that thread deferred. No two threads share a list,
	// so the lists need no lock, and nothing has to be drained inside the region.
	// Not const: OpenMP made a const variable predetermined shared before version 4.0 and does not
	// now, so a `default(none)` clause that names one is right under some compilers and wrong under
	// others.
	if (tallies.size() < ProgramOptions::NUMBER_OF_THREADS) {
		// A thread writes the tally at its own index, from inside the region, where a short vector
		// would be a silent write past the end. The caller sizes this one because the caller
		// merges it.
		throw coretools::TDevError("The field update was handed ", tallies.size(), " tallies for ",
		                           ProgramOptions::NUMBER_OF_THREADS,
		                           " threads. It needs one per thread.");
	}

	size_t n_cells = Y.total_size_of_container_space();
	std::vector<std::vector<size_t>> field_inserts(ProgramOptions::NUMBER_OF_THREADS);

#pragma omp parallel for num_threads(ProgramOptions::NUMBER_OF_THREADS)                            \
    schedule(static) default(none)                                                                 \
    shared(Y, Z_species, Z_molecule, omega, model, uniforms, tallies, field_inserts, n_cells)
	for (size_t field_cell = 0; field_cell < n_cells; ++field_cell) {
		const auto thread = static_cast<size_t>(omp_get_thread_num());
		update_cell<Policy>(field_cell, Y, Z_species, Z_molecule, omega, model, uniforms,
		                    tallies[thread], field_inserts[thread]);
	}

	// The inserts every thread deferred, in one batch. A dense field is handed empty lists and
	// does nothing with them.
	Y.insert_in_Y(field_inserts);
}

} // namespace field_update
