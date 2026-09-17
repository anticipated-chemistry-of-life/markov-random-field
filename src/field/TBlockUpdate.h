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
// A thread takes a row: every leaf pair of one species leaf. A leaf pair's Markov blanket holds no
// cell of another leaf pair -- it holds the two tree parents and the two data terms -- so the leaf
// pairs are conditionally independent and the split between threads is free. The pass is a
// parallel loop over the field's rows, and a row walks its own columns in ascending order. A
// cell's uniform is hashed from its position (ADR-0007), so neither the thread that reaches a cell
// nor the order it is reached in decides what it gets: the chain is the same at any thread count,
// and the same as the flat per-cell loop this replaced.
//
// The row is the unit because that is the granularity at which the model can answer once instead
// of per cell. A clique of one tree is named by the leaves of every other tree (ADR-0011), so the
// molecule tree's clique -- and with it its whole transition grid -- is the species leaf, constant
// down a row, while the species tree's clique is the molecule leaf and varies along it. The row
// asks the model for what it can be told once (`begin_row`), and the cell asks for the rest.
//

#pragma once

#include "cache_line_slot.h"
#include "cli.h"
#include "constants.h"
#include "coretools/Main/TError.h"
#include "coretools/Types/probability.h"
#include "field/TFieldMath.h"
#include "omp.h"
#include "random/TCellUniforms.h"
#include "storages/cell_write.h"
#include "tree/TPhylogeny.h"
#include <array>
#include <concepts>
#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

namespace block_update {

// The update names the two trees. A leaf pair is one coordinate in each. A third tree would be a
// third storage and a third factor, so this says so here rather than looping over a dimension
// count that cannot change.
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

/// The type a model's `begin_row` hands back. Named here so that the concept below can ask for an
/// lvalue of it: a row is written as it is walked, because what the row accumulates lives in it.
template<typename T>
using BlockModelRow = std::remove_cvref_t<decltype(std::declval<T &>().begin_row(size_t{},
                                                                                 size_t{}))>;

/// Everything the traversal reads that is not a state in one of its three storages.
///
/// The model is asked at two granularities. `begin_row` is asked once per species leaf, and
/// answers with whatever that row lets it settle once -- the molecule tree's clique is the species
/// leaf, so its whole transition grid is one of those things -- and with whatever that row will
/// accumulate. It is told the row's first linear index as well as the species leaf, because a
/// source whose point query walks a cursor forward seeds that cursor there. `factors` is asked per
/// leaf pair, is handed the row it belongs to and the two parent states the traversal read from the
/// node states, and answers what is left. The row's type is the model's own: the traversal only
/// carries it from one call to the next.
///
/// `record` is told what the leaf pair was given, so each data source can carry its own
/// bookkeeping forward, and it writes that bookkeeping into the row. The factors come back with
/// it, because a source that scored the cell scored both field states and only now knows which one
/// to keep. `end_row` is asked once the row's last leaf pair is drawn, and is where a row hands
/// over what it added up.
///
/// A thread takes whole rows, and everything a row accumulates is in the row. So no two threads
/// touch one object, whatever the schedule, and nothing the model sums depends on how the rows
/// were shared out.
///
/// `prepare_for_traversal` is asked once, single-threaded, before the parallel region starts
/// (`run`, below): a source readies there whatever cannot be readied concurrently -- LOTUS sorts
/// its record cache, so that each row may then seed a cursor of its own. It is told nothing about
/// the schedule. A model with nothing to ready answers with an empty function; `TBlockModel`'s
/// does, in a build without LOTUS.
template<typename T>
concept BlockModel =
    requires(T &model, BlockModelRow<T> &row, size_t species_leaf, size_t first_field_cell,
	         size_t molecule_leaf, bool species_parent, bool molecule_parent,
	         const TLeafPairFactors &factors, const field_math::TBlockStates &drawn) {
	    { model.begin_row(species_leaf, first_field_cell) } -> std::same_as<BlockModelRow<T>>;
	    {
		    model.factors(row, molecule_leaf, species_parent, molecule_parent)
	    } -> std::same_as<TLeafPairFactors>;
	    { model.record(row, molecule_leaf, factors, drawn) } -> std::same_as<void>;
	    { model.end_row(row) } -> std::same_as<void>;
	    { model.prepare_for_traversal() } -> std::same_as<void>;
    };

/// What every leaf pair of one species leaf shares, settled once before the row is walked.
///
/// The traversal's own half is here; the model's half is `model_row`, whose type is the model's
/// business and not this file's. Between them they are what a cell no longer has to work out for
/// itself.
template<typename ModelRow> struct TRow {
	/// The species leaf this row is, which is the first coordinate of every leaf pair in it.
	size_t species_leaf = 0;
	/// The species node whose state is the species parent of every leaf pair in this row. A leaf
	/// is never a root, so it has one.
	size_t species_parent_node = 0;
	/// The field's linear index of this row's first cell. The field is row-major over
	/// (species leaf, molecule leaf), so a cell's index is this plus its molecule leaf.
	size_t first_field_cell = 0;
	/// What the model settled for this row.
	ModelRow model_row;
};

/// Draws one leaf pair: the field cell and both tree field cells at it.
///
/// It reads five cells -- the three the draw moves, and the two tree parents -- and writes back the
/// ones the draw changed. A leaf pair sits at the same multidimensional index in all three
/// containers (ADR-0005), so one subscript names its cell in each; the linear indices differ,
/// because the container shapes do.
///
/// The leaf pair is named by its row and its molecule leaf, and the field's linear index follows
/// from the row rather than from a division: a linear index is `row * n_columns + column`. It
/// names no cell of either node state -- the molecule node state has one column per molecule
/// *node* where the field has one per leaf -- so the two node states are addressed by the leaf
/// pair's subscript instead, and each linearises it against its own dimensions.
///
/// A write reaches the storage in place, or comes back as a deferred insert. That is the only exit
/// a write inside a parallel region may take (ADR-0006).
template<field_math::LinkPolicy Policy, typename Field, typename NodeState, BlockModel Model,
         typename ModelRow>
void update_cell(TRow<ModelRow> &row, size_t molecule_leaf, Field &Y, NodeState &Z_species,
                 NodeState &Z_molecule, const TPhylogeny &molecule,
                 const field_math::TErrorProbability &omega, Model &model,
                 const TCellUniforms &uniforms, field_math::TLinkCounters &counters,
                 std::vector<size_t> &field_inserts, std::vector<size_t> &species_inserts,
                 std::vector<size_t> &molecule_inserts) {
	const IndexArray leaf_pair{row.species_leaf, molecule_leaf};
	const size_t field_cell = row.first_field_cell + molecule_leaf;
	// The row's arithmetic is the field's own linearisation, or the cell drawn here is not the
	// cell the uniform was hashed for.
	DEBUG_ASSERT(field_cell == Y.get_linear_index_in_container_space(leaf_pair));

	// The species parent is the cell above this one in the species node state -- its node is the
	// row's, and only the column moves. The molecule parent is the cell beside it in the molecule
	// node state. A leaf is never a root, so a parent is an internal node, and an internal node is
	// never a leaf -- so no leaf pair of this update writes what this one reads there.
	const bool species_parent =
	    Z_species.is_one(IndexArray{row.species_parent_node, molecule_leaf});
	const bool molecule_parent =
	    Z_molecule.is_one(IndexArray{row.species_leaf, molecule.parent_of(molecule_leaf)});

	const TLeafPairFactors factors =
	    model.factors(row.model_row, molecule_leaf, species_parent, molecule_parent);

	const field_math::TBlockStates current{.y   = Y.is_one(field_cell),
	                                       .z_s = Z_species.is_one(leaf_pair),
	                                       .z_m = Z_molecule.is_one(leaf_pair)};

	// The cell names the uniform it draws, so the thread that happens to reach this cell does not
	// decide what it gets.
	const auto draw = field_math::draw_block<Policy>(
	    factors.prob_z_s_is_one, factors.prob_z_m_is_one, omega, factors.lotus,
	    factors.simple_error, current, uniforms.at(field_cell));

	// Write each cell back where it was read from. A storage that holds the cell takes the write in
	// place, and one that does not hands the cell back for a later insert.
	if (current.y != draw.drawn.y) { write_or_defer(Y, field_cell, draw.drawn.y, field_inserts); }
	if (current.z_s != draw.drawn.z_s) {
		write_or_defer(Z_species, leaf_pair, draw.drawn.z_s, species_inserts);
	}
	if (current.z_m != draw.drawn.z_m) {
		write_or_defer(Z_molecule, leaf_pair, draw.drawn.z_m, molecule_inserts);
	}

	// The counters and the model's bookkeeping are the row's own, on this thread's stack. They
	// used to be a slot in a vector indexed by `omp_get_thread_num()`, written here once per leaf
	// pair: two threads' slots then sat in one cache line, and a write by one took the line from
	// the other. That carried 45% of this function's profile samples.
	counters.add(draw.to.bucket, draw.to.y);
	model.record(row.model_row, molecule_leaf, factors, draw.drawn);
}

/// Draws one row: every leaf pair of one species leaf, in ascending molecule-leaf order.
///
/// What the row settles once is what a per-cell loop was asking for over and over: the species
/// leaf's parent node and its branch, the molecule tree's clique -- which *is* the species leaf,
/// so its transition grid is constant down the row -- and the row's first linear index. The
/// ascending order is not cosmetic: it is what lets a data source answer a point query from a
/// forward cursor (`TLotus::holds_a_record`) rather than by hashing every cell.
///
/// The row is also the unit the tally and the model's bookkeeping are accumulated in. Both live
/// here, on the calling thread's stack, for as long as the row lasts, and are handed over once
/// when it ends -- `row_counters` into this row's slot, the model's into the model's own.
template<field_math::LinkPolicy Policy, typename Field, typename NodeState, BlockModel Model>
void update_row(size_t species_leaf, size_t n_molecule_leaves, Field &Y, NodeState &Z_species,
                NodeState &Z_molecule, const TPhylogeny &species, const TPhylogeny &molecule,
                const field_math::TErrorProbability &omega, Model &model,
                const TCellUniforms &uniforms, field_math::TLinkCounters &row_counters,
                std::vector<size_t> &field_inserts, std::vector<size_t> &species_inserts,
                std::vector<size_t> &molecule_inserts) {
	const size_t first_field_cell = species_leaf * n_molecule_leaves;
	TRow<BlockModelRow<Model>> row{.species_leaf        = species_leaf,
	                               .species_parent_node = species.parent_of(species_leaf),
	                               .first_field_cell    = first_field_cell,
	                               .model_row = model.begin_row(species_leaf, first_field_cell)};

	field_math::TLinkCounters counters;
	for (size_t molecule_leaf = 0; molecule_leaf < n_molecule_leaves; ++molecule_leaf) {
		update_cell<Policy>(row, molecule_leaf, Y, Z_species, Z_molecule, molecule, omega, model,
		                    uniforms, counters, field_inserts, species_inserts, molecule_inserts);
	}
	row_counters.merge(counters);
	model.end_row(row.model_row);
}

/// One block update: every leaf pair of the field's container space, one row to a thread.
///
/// `per_row_counters` is one entry per field row, and the caller merges them in row order. It is
/// not cleared here: a caller that hands over used counters gets their counts added to, which is
/// what makes the merge its business rather than this one's.
///
/// One slot per row, and not one per thread. A row writes its slot once, when it ends, so two
/// rows sharing a cache line costs one write per row rather than one per leaf pair -- and a merge
/// in row order does not depend on how the rows were shared out. The counters are integers and so
/// would sum the same either way; the model's likelihood is not, and it takes the same treatment
/// for the reason `TTree::log_node_state_density` gives.
///
/// Because nothing here is indexed by thread and nothing the model sums depends on the schedule,
/// the schedule is this file's own business again. It is a bare `schedule(static)`: every row
/// costs the same, so an even split is a balanced one, whatever size team the runtime gives us.
/// The chunked schedule this replaces had to be computed from `ProgramOptions::NUMBER_OF_THREADS`
/// because LOTUS seeded one cursor per thread from it -- and a runtime that handed back fewer
/// threads than were asked for then ran the loop in two rounds, leaving the second half idle.
template<field_math::LinkPolicy Policy, typename Field, typename NodeState, BlockModel Model>
void run(Field &Y, NodeState &Z_species, NodeState &Z_molecule, const TPhylogeny &species,
         const TPhylogeny &molecule, const field_math::TErrorProbability &omega, Model &model,
         const TCellUniforms &uniforms,
         std::vector<field_math::TLinkCounters> &per_row_counters) {
	// The field is row-major over (species leaf, molecule leaf): one row per species leaf, one
	// column per molecule leaf.
	// Not const: OpenMP made a const variable predetermined shared before version 4.0 and does not
	// now, so a `default(none)` clause that names one is right under some compilers and wrong under
	// others.
	size_t n_rows    = Y.dimensions()[0];
	size_t n_columns = Y.dimensions()[1];
	if (per_row_counters.size() != n_rows) {
		throw coretools::TDevError("The block update was given ", per_row_counters.size(),
		                           " counter slots for a field of ", n_rows, " rows.");
	}

	// Whatever a source cannot ready from inside the parallel region. It is told nothing about the
	// schedule: a row seeds its own cursor.
	model.prepare_for_traversal();

	// One list per thread and per container, filled by the writes that thread deferred. No two
	// threads share a list, so the lists need no lock, and nothing has to be drained inside the
	// region. These stay per thread: they are an output buffer and not an accumulator, so nothing
	// about them is order-sensitive, and only a storage that refuses a write appends to one -- a
	// dense one never does. Padded all the same, because the sparse backend appends often and the
	// three `std::vector` headers of two threads would otherwise share a line.
	std::vector<TCacheLineSlot<std::vector<size_t>>> field_inserts(
	    ProgramOptions::NUMBER_OF_THREADS);
	std::vector<TCacheLineSlot<std::vector<size_t>>> species_inserts(
	    ProgramOptions::NUMBER_OF_THREADS);
	std::vector<TCacheLineSlot<std::vector<size_t>>> molecule_inserts(
	    ProgramOptions::NUMBER_OF_THREADS);

#pragma omp parallel for num_threads(ProgramOptions::NUMBER_OF_THREADS) schedule(static)           \
    default(none) shared(Y, Z_species, Z_molecule, species, molecule, omega, model, uniforms,      \
	                         per_row_counters, field_inserts, species_inserts, molecule_inserts,   \
	                         n_rows, n_columns)
	for (size_t species_leaf = 0; species_leaf < n_rows; ++species_leaf) {
		const auto thread = static_cast<size_t>(omp_get_thread_num());
		update_row<Policy>(species_leaf, n_columns, Y, Z_species, Z_molecule, species, molecule,
		                   omega, model, uniforms, per_row_counters[species_leaf],
		                   field_inserts[thread].value, species_inserts[thread].value,
		                   molecule_inserts[thread].value);
	}

	// The inserts every thread deferred, in one batch per container. A dense container is handed
	// empty lists and does nothing with them.
	std::vector<std::vector<size_t>> field_batch(field_inserts.size());
	std::vector<std::vector<size_t>> species_batch(species_inserts.size());
	std::vector<std::vector<size_t>> molecule_batch(molecule_inserts.size());
	for (size_t thread = 0; thread < field_inserts.size(); ++thread) {
		field_batch[thread]    = std::move(field_inserts[thread].value);
		species_batch[thread]  = std::move(species_inserts[thread].value);
		molecule_batch[thread] = std::move(molecule_inserts[thread].value);
	}
	Y.insert_in_Y(field_batch);
	Z_species.insert_in_Z(species_batch);
	Z_molecule.insert_in_Z(molecule_batch);
}

} // namespace block_update
