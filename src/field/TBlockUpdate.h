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
// A thread takes a band of rows: every leaf pair of a few dozen species leaves. A leaf pair's
// Markov blanket holds no cell of another leaf pair -- it holds the two tree parents and the two
// data terms -- so the leaf pairs are conditionally independent and the split between threads is
// free. The pass is a parallel loop over bands; a band is drawn a tile of columns at a time, row
// after row, and each row walks its own columns in ascending order. A cell's uniform is hashed from
// its position (ADR-0007), so neither the thread that reaches a cell nor the order it is reached in
// decides what it gets: the chain is the same at any thread count, and the same as the flat
// per-cell loop this replaced.
//
// The row is the unit of what is settled and added up because that is the granularity at which the
// model can answer once instead of per cell. A clique of one tree is named by the leaves of every
// other tree (ADR-0011), so the molecule tree's clique -- and with it its whole transition grid --
// is the species leaf, constant down a row, while the species tree's clique is the molecule leaf
// and varies along it. The row asks the model for what it can be told once (`begin_row`), and the
// cell asks for the rest. The band and the tile are only an order of visits: the species node
// state keeps its cliques contiguous (ADR-0013), so a field row reads it across a stride, and the
// tile is what lets the rows of a band share the cache lines that stride lands on.
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
#include <algorithm>
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
/// process says about its own tree field cell, and what the data sources say about the field cell.
struct TLeafPairFactors {
	/// P(Z_s = 1 | the species parent's state), from the species tree's transition grid.
	coretools::Probability prob_z_s_is_one;
	/// P(Z_m = 1 | the molecule parent's state), from the molecule tree's transition grid.
	coretools::Probability prob_z_m_is_one;
	/// The cell data factor: {P(obs | Y = 0), P(obs | Y = 1)}, every compiled-in data source's
	/// answer for this cell multiplied together.
	field_math::TCellDataFactor cell_data;
};

/// The type a model's `begin_row` hands back. Named here so that the concept below can ask for an
/// lvalue of it: a row is written as it is walked, because what the row accumulates lives in it.
template<typename T>
using BlockModelRow =
    std::remove_cvref_t<decltype(std::declval<T &>().begin_row(size_t{}, size_t{}))>;

/// The type a model's `begin_column` hands back. A column is only read once it is made.
template<typename T>
using BlockModelColumn = std::remove_cvref_t<decltype(std::declval<T &>().begin_column(size_t{}))>;

/// Everything the traversal reads that is not a state in one of its three storages.
///
/// The model is asked at three granularities. `begin_row` is asked once per species leaf, and
/// answers with whatever that row lets it settle once -- the molecule tree's clique is the species
/// leaf, so its whole transition grid is one of those things -- and with whatever that row will
/// accumulate. It is told the row's first linear index as well as the species leaf, because a
/// source whose point query walks a cursor forward seeds that cursor there. `begin_column` is the
/// same for a molecule leaf: the species tree's clique is the molecule leaf, so its transition grid
/// is the column's, and so is the molecule leaf's branch. It is asked once per column of a tile,
/// and the column it answers is read by every row of the band. `factors` is asked per leaf pair,
/// is handed the row and the column it belongs to and the two parent states the traversal read
/// from the node states, and answers what is left. The row's and the column's types are the
/// model's own: the traversal only carries them from one call to the next.
///
/// `record` is told what the leaf pair was given, so each data source can carry its own
/// bookkeeping forward, and it writes that bookkeeping into the row. The factors come back with
/// it, because a source that scored the cell scored both field states and only now knows which one
/// to keep. `end_row` is asked once the row's last leaf pair is drawn, and is where a row hands
/// over what it added up.
///
/// A thread takes whole rows, and everything a row accumulates is in the row. So no two threads
/// touch one object, whatever the schedule, and nothing the model sums depends on how the rows
/// were shared out. A column accumulates nothing, which is what lets the rows of a band share it.
///
/// `prepare_for_traversal` is asked once, single-threaded, before the parallel region starts
/// (`run`, below): a source readies there whatever cannot be readied concurrently -- LOTUS sorts
/// its record cache, so that each row may then seed a cursor of its own. It is told nothing about
/// the schedule. A model with nothing to ready answers with an empty function; `TBlockModel`'s
/// does, in a build without LOTUS.
template<typename T>
concept BlockModel =
    requires(T &model, BlockModelRow<T> &row, const BlockModelColumn<T> &column,
	         size_t species_leaf, size_t first_field_cell, size_t molecule_leaf,
	         bool species_parent, bool molecule_parent, const TLeafPairFactors &factors,
	         const field_math::TBlockStates &drawn) {
	    { model.begin_row(species_leaf, first_field_cell) } -> std::same_as<BlockModelRow<T>>;
	    { model.begin_column(molecule_leaf) } -> std::same_as<BlockModelColumn<T>>;
	    {
		    model.factors(row, column, molecule_leaf, species_parent, molecule_parent)
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
///
/// Where a cell of a node state sits is a sum of one term per subscript -- the subscript times
/// that dimension's stride, under any layout (ADR-0013) -- so the row keeps its own terms and a
/// cell adds the column's to them. The row's term of a node state is where its cell at molecule
/// leaf (or node) 0 sits.
template<typename ModelRow> struct TRow {
	/// The species leaf this row is, which is the first coordinate of every leaf pair in it.
	size_t species_leaf            = 0;
	/// The field's linear index of this row's first cell. The field is row-major over
	/// (species leaf, molecule leaf), so a cell's index is this plus its molecule leaf.
	size_t first_field_cell        = 0;
	/// The row's term of a position in the species node state: the species leaf's.
	size_t species_position        = 0;
	/// The same for the species leaf's parent node, whose state is the species parent of every
	/// leaf pair in this row. A leaf is never a root, so it has one.
	size_t species_parent_position = 0;
	/// The row's term of a position in the molecule node state.
	size_t molecule_position       = 0;
	/// What the model settled for this row.
	ModelRow model_row;
};

/// What every leaf pair of one molecule leaf shares, settled once per tile and read by every row of
/// the band. The traversal's own half is here, and is `TRow`'s turned around; the model's half is
/// `model_column`.
template<typename ModelColumn> struct TColumn {
	/// The molecule leaf this column is, which is the second coordinate of every leaf pair in it.
	size_t molecule_leaf            = 0;
	/// The column's term of a position in the species node state.
	size_t species_position         = 0;
	/// The column's term of a position in the molecule node state: the molecule leaf's.
	size_t molecule_position        = 0;
	/// The same for the molecule leaf's parent node, whose state is the molecule parent of every
	/// leaf pair in this column.
	size_t molecule_parent_position = 0;
	/// What the model settled for this column.
	ModelColumn model_column;
};

/// Draws one leaf pair: the field cell and both tree field cells at it.
///
/// It reads five cells -- the three the draw moves, and the two tree parents -- and writes back the
/// ones the draw changed. A leaf pair sits at the same multidimensional index in all three
/// containers (ADR-0005), so one subscript names its cell in each; the linear indices differ,
/// because the container shapes do.
///
/// The leaf pair is named by its row and its column. The field's linear index follows from the row
/// rather than from a division: a linear index is `row * n_columns + column`. It names no cell of
/// either node state -- the molecule node state has one column per molecule *node* where the field
/// has one per leaf -- so the four node-state cells read here are found by adding the row's term of
/// their position to the column's.
///
/// A write reaches the storage in place, or comes back as a deferred insert. That is the only exit
/// a write inside a parallel region may take (ADR-0006).
template<field_math::LinkPolicy Policy, typename Field, typename NodeState, BlockModel Model,
         typename ModelRow, typename ModelColumn>
void update_cell(TRow<ModelRow> &row, const TColumn<ModelColumn> &column, Field &Y,
                 NodeState &Z_species, NodeState &Z_molecule,
                 const field_math::TErrorProbability &omega, Model &model,
                 const TCellUniforms &uniforms, field_math::TLinkCounters &counters,
                 std::vector<size_t> &field_inserts, std::vector<size_t> &species_inserts,
                 std::vector<size_t> &molecule_inserts) {
	const size_t molecule_leaf = column.molecule_leaf;
	const IndexArray leaf_pair{row.species_leaf, molecule_leaf};
	const size_t field_cell = row.first_field_cell + molecule_leaf;
	// The row's arithmetic is the field's own linearisation, or the cell drawn here is not the
	// cell the uniform was hashed for. The same goes for the node states' positions.
	DEBUG_ASSERT(field_cell == Y.get_linear_index_in_container_space(leaf_pair));
	DEBUG_ASSERT(row.species_position + column.species_position ==
	             Z_species.position_of(leaf_pair));
	DEBUG_ASSERT(row.molecule_position + column.molecule_position ==
	             Z_molecule.position_of(leaf_pair));

	// The species parent is the cell above this one in the species node state -- its node is the
	// row's, and only the column moves. The molecule parent is the cell beside it in the molecule
	// node state -- its node is the column's, and only the row moves. A leaf is never a root, so a
	// parent is an internal node, and an internal node is never a leaf -- so no leaf pair of this
	// update writes what this one reads there.
	const bool species_parent =
	    Z_species.is_one_at(row.species_parent_position + column.species_position);
	const bool molecule_parent =
	    Z_molecule.is_one_at(row.molecule_position + column.molecule_parent_position);

	const TLeafPairFactors factors = model.factors(row.model_row, column.model_column,
	                                               molecule_leaf, species_parent, molecule_parent);

	const field_math::TBlockStates current{
	    .y   = Y.is_one(field_cell),
	    .z_s = Z_species.is_one_at(row.species_position + column.species_position),
	    .z_m = Z_molecule.is_one_at(row.molecule_position + column.molecule_position)};

	// The cell names the uniform it draws, so the thread that happens to reach this cell does not
	// decide what it gets.
	const auto draw =
	    field_math::draw_block<Policy>(factors.prob_z_s_is_one, factors.prob_z_m_is_one, omega,
		                               factors.cell_data, current, uniforms.at(field_cell));

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

/// How many species leaves a band of the update holds, and how many molecule leaves one tile of a
/// band spans.
///
/// A tile is drawn row after row, and every row of it reads the same columns of the species node
/// state. That node state keeps its own cliques contiguous (ADR-0013), so a column is a molecule
/// leaf and the band's leaves at one column sit side by side, a cache line or two for the whole
/// band. Drawn one whole row at a time, the next row would come back to a line long gone; drawn a
/// tile at a time, it finds the line the row above it just used. The tile is sized so that those
/// lines, and the lines holding the band's parents, stay in a core's L2 between rows.
inline constexpr size_t BAND_ROWS    = 64;
inline constexpr size_t TILE_COLUMNS = 128;

/// Draws one band: every leaf pair of the species leaves `[first_row, end_row)`, a tile at a time.
///
/// Every row still walks its own columns in ascending order, and still settles what it can once
/// and accumulates what it adds up in a row of its own -- the band only interleaves the rows.
/// What the row settles once is what a per-cell loop was asking for over and over: the species
/// leaf's parent node and its branch, the molecule tree's clique -- which *is* the species leaf,
/// so its transition grid is constant down the row -- and the row's first linear index. A tile's
/// columns are settled the same way, once per tile rather than once per row: the molecule leaf's
/// parent node and its branch, and the species tree's clique -- which *is* the molecule leaf. The
/// ascending order is not cosmetic: it is what lets a data source answer a point query from a
/// forward cursor (LOTUS's pass, lotus/block_pass.h) rather than by hashing every cell.
///
/// No leaf pair reads a cell another one writes, so the interleaving draws what the rows drawn one
/// after another would: each row's tally and bookkeeping see the same cells in the same order, and
/// each is handed over once when the band ends -- the tally into the row's slot, the model's into
/// the model's own.
template<field_math::LinkPolicy Policy, typename Field, typename NodeState, BlockModel Model>
void update_band(size_t first_row, size_t end_row, size_t n_molecule_leaves, Field &Y,
                 NodeState &Z_species, NodeState &Z_molecule, const TPhylogeny &species,
                 const TPhylogeny &molecule, const field_math::TErrorProbability &omega,
                 Model &model, const TCellUniforms &uniforms,
                 std::vector<field_math::TLinkCounters> &per_row_counters,
                 std::vector<size_t> &field_inserts, std::vector<size_t> &species_inserts,
                 std::vector<size_t> &molecule_inserts) {
	// A row's and a column's term of a position: where the cell at subscript 0 in the other
	// dimension sits.
	const auto row_term = [](const NodeState &Z, size_t row) {
		return Z.position_of(IndexArray{row, 0});
	};
	const auto column_term = [](const NodeState &Z, size_t column) {
		return Z.position_of(IndexArray{0, column});
	};

	std::vector<TRow<BlockModelRow<Model>>> rows;
	rows.reserve(end_row - first_row);
	for (size_t species_leaf = first_row; species_leaf < end_row; ++species_leaf) {
		const size_t first_field_cell = species_leaf * n_molecule_leaves;
		rows.push_back(
		    {.species_leaf            = species_leaf,
			 .first_field_cell        = first_field_cell,
			 .species_position        = row_term(Z_species, species_leaf),
			 .species_parent_position = row_term(Z_species, species.parent_of(species_leaf)),
			 .molecule_position       = row_term(Z_molecule, species_leaf),
			 .model_row               = model.begin_row(species_leaf, first_field_cell)});
	}
	std::vector<field_math::TLinkCounters> counters(rows.size());
	std::vector<TColumn<BlockModelColumn<Model>>> columns;
	columns.reserve(TILE_COLUMNS);

	for (size_t first_column = 0; first_column < n_molecule_leaves; first_column += TILE_COLUMNS) {
		const size_t end_column = std::min(first_column + TILE_COLUMNS, n_molecule_leaves);
		columns.clear();
		for (size_t molecule_leaf = first_column; molecule_leaf < end_column; ++molecule_leaf) {
			columns.push_back({.molecule_leaf     = molecule_leaf,
			                   .species_position  = column_term(Z_species, molecule_leaf),
			                   .molecule_position = column_term(Z_molecule, molecule_leaf),
			                   .molecule_parent_position =
			                       column_term(Z_molecule, molecule.parent_of(molecule_leaf)),
			                   .model_column = model.begin_column(molecule_leaf)});
		}
		for (size_t r = 0; r < rows.size(); ++r) {
			for (const auto &column : columns) {
				update_cell<Policy>(rows[r], column, Y, Z_species, Z_molecule, omega, model,
				                    uniforms, counters[r], field_inserts, species_inserts,
				                    molecule_inserts);
			}
		}
	}

	for (size_t r = 0; r < rows.size(); ++r) {
		per_row_counters[rows[r].species_leaf] = counters[r];
		model.end_row(rows[r].model_row);
	}
}

/// One block update: every leaf pair of the field's container space, one band of rows to a thread
/// at a time. It returns the link counters it tallied over every leaf pair. The update covers
/// every leaf pair, so the tally is that of the configuration it leaves.
///
/// Because nothing here is indexed by thread and nothing the model sums depends on the schedule,
/// the schedule is this file's own business. It is `schedule(dynamic)` over bands: there are a
/// few hundred of them, too few for a static split to come out even, and a band holding a species
/// with mass-spec runs costs more than one without.
template<field_math::LinkPolicy Policy, typename Field, typename NodeState, BlockModel Model>
[[nodiscard]] field_math::TLinkCounters run(Field &Y, NodeState &Z_species, NodeState &Z_molecule,
                                            const TPhylogeny &species, const TPhylogeny &molecule,
                                            const field_math::TErrorProbability &omega,
                                            Model &model, const TCellUniforms &uniforms) {
	// The field is row-major over (species leaf, molecule leaf): one row per species leaf, one
	// column per molecule leaf.
	// Not const: OpenMP made a const variable predetermined shared before version 4.0 and does not
	// now, so a `default(none)` clause that names one is right under some compilers and wrong under
	// others.
	size_t n_rows    = Y.dimensions()[0];
	size_t n_columns = Y.dimensions()[1];

	// The link counters keep one slot per field row, not one per thread. A row writes its slot
	// once, when it ends. So two rows that share a cache line cost one write per row, not one per
	// leaf pair. This function merges the slots in row order after the region. Nothing in the
	// merge depends on the schedule, so the counters are the same at any thread count.
	std::vector<field_math::TLinkCounters> per_row_counters(n_rows);

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

	size_t n_bands = (n_rows + BAND_ROWS - 1) / BAND_ROWS;
#pragma omp parallel for num_threads(ProgramOptions::NUMBER_OF_THREADS)                            \
    schedule(dynamic) default(none)                                                                \
    shared(Y, Z_species, Z_molecule, species, molecule, omega, model, uniforms, per_row_counters,  \
	           field_inserts, species_inserts, molecule_inserts, n_rows, n_columns, n_bands)
	for (size_t band = 0; band < n_bands; ++band) {
		const auto thread      = static_cast<size_t>(omp_get_thread_num());
		const size_t first_row = band * BAND_ROWS;
		update_band<Policy>(first_row, std::min(first_row + BAND_ROWS, n_rows), n_columns, Y,
		                    Z_species, Z_molecule, species, molecule, omega, model, uniforms,
		                    per_row_counters, field_inserts[thread].value,
		                    species_inserts[thread].value, molecule_inserts[thread].value);
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

	field_math::TLinkCounters counters;
	for (const auto &row : per_row_counters) { counters.merge(row); }
	return counters;
}

} // namespace block_update
