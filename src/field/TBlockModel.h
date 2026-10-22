//
// What the block update asks about a leaf pair: the two trees' processes, and the data sources.
//
// The traversal in field/TBlockUpdate.h reads three states and two tree parents from its storages.
// It asks this for the rest. Each tree says what its own tree field cell should be, given that
// tree's parent. The data says what it makes of the field cell, as one cell data factor
// (GLOSSARY.md). Every answer is about one leaf pair, because that is the unit the traversal works
// in.
//
// This is the whole of what the traversal must not know about. So the traversal keeps no tree, no
// clique and no data source, and a test can run its loop against a model of its own.
//
// The data sources that bring their own pass through the block update are folded by
// field/TBlockData.h, which this holds (ADR-0017). LOTUS and mass spectrometry do not bring one
// yet: their factors are multiplied in here, each at its own position in the fixed order (LOTUS,
// the passes, mass spectrometry), and LOTUS's present-cell bookkeeping is kept here too.
//

#pragma once

#include "TDataModel.h"
#include "Types.h"
#include "constants.h"
#include "coretools/Types/probability.h"
#include "field/TBlockUpdate.h"
#include "field/TFieldMath.h"
#include "field/TBlockData.h"
#include "omp.h"
#include "tree/TTree.h"
#include "tree/branch/TTransitionGrid.h"
#include "tree/branch/TTransitionGridTable.h"
#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

//-----------------------------------
// LOTUS's bookkeeping
//-----------------------------------

/// What the rows of a block update hand LOTUS, committed to it when the update ends.
///
/// What a row owes LOTUS is not a number but a row of the field's present-cell counts: a slot per
/// row would be a copy of the counts per species leaf. Rows add into the one set of counts instead,
/// entry by entry and atomically. Its entries are integers, so the counts the update ends with are
/// the same whichever rows got there first (ADR-0015).
///
/// A build without LOTUS holds nothing here, so that the caller is the same in every build.
class TDataUpdateAccumulator {
private:
#ifdef USE_LOTUS
	lotus_math::TPresentCells _lotus_present_cells;
#endif

public:
	explicit TDataUpdateAccumulator([[maybe_unused]] const TDataModel &data_model)
#ifdef USE_LOTUS
	    : _lotus_present_cells(data_model.get_lotus().no_present_cells())
#endif
	{
	}

#ifdef USE_LOTUS
	/// Adds one row's present cells to the field's: `present[reported * n_molecule_classes +
	/// molecule_class]`, for a row whose species leaf is in `species_class`. Called once per row,
	/// from inside the parallel region, by rows that may share a species class.
	void add_lotus_row(uint32_t species_class, const std::vector<uint32_t> &present) {
		auto &cells                     = _lotus_present_cells;
		const size_t n_molecule_classes = cells.n_molecule_classes;
		for (size_t reported = 0; reported < 2; ++reported) {
			for (uint32_t m = 0; m < n_molecule_classes; ++m) {
				const uint32_t n = present[reported * n_molecule_classes + m];
				if (n == 0) { continue; }
				std::atomic_ref<uint64_t>(cells.ones[reported][cells.slot(species_class, m)])
				    .fetch_add(n, std::memory_order_relaxed);
			}
		}
	}
#endif

	/// Installs what the rows added up in LOTUS. Called once, after the parallel region.
	void commit([[maybe_unused]] TDataModel &data_model) {
#ifdef USE_LOTUS
		data_model.get_lotus().install_present_cells(std::move(_lotus_present_cells));
#endif
	}
};

//-----------------------------------
// The model
//-----------------------------------

/// The two trees and the data sources, as the block update wants them.
///
/// Both trees have something to say from the first update on, because the chain is given a
/// starting node state before it (TMarkovField::_start_the_chain).
///
/// An inferred chain is the only one that has a model to ask. A simulated one draws its whole
/// configuration forward and runs no update (TMarkovField::simulate).
///
/// Every thread of the update asks this one object. `begin_row` hands back a row that carries
/// everything the walk down that row accumulates, so `factors` and `record` touch the caller's own
/// row and nothing another thread also holds; `end_row` hands the row's tallies over, once.
class TBlockModel {
private:
	const TTree &_species_tree;
	const TTree &_molecule_tree;
	TDataModel &_data_model;
	TCompiledBlockData &_data;
	TDataUpdateAccumulator &_accumulator;

public:
	TBlockModel(const std::vector<std::unique_ptr<TTree>> &trees, TDataModel &data_model,
	            TCompiledBlockData &data, TDataUpdateAccumulator &accumulator)
	    : _species_tree(*trees.front()), _molecule_tree(*trees.back()), _data_model(data_model),
	      _data(data), _accumulator(accumulator) {}

	/// Told once, before the parallel region starts, that a traversal is about to begin. Every
	/// data source readies whatever cannot be readied concurrently -- LOTUS sorts its record cache
	/// here, so that every row may then seed a cursor of its own. It is told nothing about how the
	/// rows are shared out: that is what used to tie the traversal's schedule to
	/// `ProgramOptions::NUMBER_OF_THREADS`.
	void prepare_for_traversal() {
#ifdef USE_LOTUS
		_data_model.get_lotus().prepare_for_block_update();
#endif
		_data.begin_pass();
	}

	/// P(a tree field cell = 1 | its parent's state), under one clique's process on one branch.
	template<TransitionGridLike Process>
	[[nodiscard]] static coretools::Probability
	prob_of_one(const Process &process, TypeBinnedBranchLengths branch, bool parent_state) {
		return coretools::P(process.probability(branch, parent_state, /*to=*/true));
	}

	/// What one field row -- one species leaf, every molecule leaf -- settles before it is walked.
	///
	/// A clique of one tree is named by a leaf of every other tree (ADR-0011), so the molecule
	/// tree's clique *is* the species leaf: its process is the same for every cell of the row, and
	/// asking the table for it once per row rather than once per cell is the whole point of this
	/// type. The species tree's branch is the row's too -- a species leaf sits on one branch,
	/// whichever molecule leaf the cell pairs it with. What is left per cell is the species tree's
	/// clique (the molecule leaf) and the molecule tree's branch (likewise).
	/// It also carries what the row accumulates: each data source's row state, and LOTUS's forward
	/// cursor into its record cache. Both used to be one slot per thread, indexed by
	/// `omp_get_thread_num()` and written once per leaf pair, which put two threads' slots in one
	/// cache line and made the likelihood's rounding a function of the schedule.
	struct TRow {
		/// The molecule tree's process for this row.
		TTransitionGridView molecule_process;
		/// The bin the species leaf's branch sits in.
		TypeBinnedBranchLengths species_branch;
		/// The species leaf the row is, which `factors` needs to name a cell.
		size_t species_leaf;
		/// The molecule tree's clique this row is, kept so a debug build can check that it really
		/// does not move down the row.
		size_t molecule_clique;
		/// Every pass's state for this row.
		TCompiledBlockData::TRow data;
#ifdef USE_LOTUS
		/// This row's forward cursor into LOTUS's sorted record cache, seeded at the row's first
		/// cell. The row walks its columns in ascending order, which is all the cursor asks.
		TLotus::RecordCursor lotus_cursor;
		/// The paper-count class of the row's species leaf: the row of the counts it adds to.
		uint32_t lotus_species_class;
		/// The row's present cells, `[reported * n_molecule_classes + molecule_class]`. A row
		/// has at most a molecule leaf's worth, so a 32-bit count holds it.
		std::vector<uint32_t> lotus_present;
#endif
	};

	[[nodiscard]] TRow begin_row(size_t species_leaf, size_t first_field_cell) {
		// The molecule tree drops the molecule coordinate, so any molecule leaf names this row's
		// clique; 0 is the one every field row has.
		const IndexArray first_cell_of_row{species_leaf, 0};
		return TRow{
		    .molecule_process = _molecule_tree.transition_grid_of_cell(first_cell_of_row),
		    .species_branch   = _species_tree.get_binned_branch_length(species_leaf),
		    .species_leaf     = species_leaf,
		    .molecule_clique  = _molecule_tree.clique_of_cell(first_cell_of_row),
		    .data             = _data.begin_row(species_leaf, first_field_cell),
#ifdef USE_LOTUS
		    .lotus_cursor        = _data_model.get_lotus().cursor_for_row(first_field_cell),
		    .lotus_species_class = _data_model.get_lotus().species_class(species_leaf),
		    .lotus_present = std::vector<uint32_t>(2 * _data_model.get_lotus().n_molecule_classes()),
#endif
		};
	}

	/// What one tile column -- one molecule leaf, every species leaf of a band -- settles before
	/// the band's rows are walked across it. The species tree's clique is the molecule leaf, so its
	/// process is the column's, and so is the molecule leaf's branch. Nothing here accumulates, so
	/// every row of the band reads the one column.
	struct TColumn {
		/// The species tree's process for this column.
		TTransitionGridView species_process;
		/// The bin the molecule leaf's branch sits in.
		TypeBinnedBranchLengths molecule_branch;
		/// The species tree's clique this column is, kept so a debug build can check that it really
		/// does not move down the column.
		size_t species_clique;
	};

	[[nodiscard]] TColumn begin_column(size_t molecule_leaf) const {
		// The species tree drops the species coordinate, so any species leaf names this column's
		// clique; 0 is the one every column has.
		const IndexArray first_cell_of_column{0, molecule_leaf};
		return TColumn{
		    .species_process = _species_tree.transition_grid_of_cell(first_cell_of_column),
		    .molecule_branch = _molecule_tree.get_binned_branch_length(molecule_leaf),
		    .species_clique  = _species_tree.clique_of_cell(first_cell_of_column),
		};
	}

	[[nodiscard]] block_update::TLeafPairFactors factors(TRow &row, const TColumn &column,
	                                                     size_t molecule_leaf, bool species_parent,
	                                                     bool molecule_parent) const {
		[[maybe_unused]] const IndexArray cell{row.species_leaf, molecule_leaf};
		// What the row and the column settled once has to be what this cell would have asked for.
		DEBUG_ASSERT(_molecule_tree.clique_of_cell(cell) == row.molecule_clique);
		DEBUG_ASSERT(_species_tree.clique_of_cell(cell) == column.species_clique);

		block_update::TLeafPairFactors leaf_pair;
		// Each tree's process comes from where its clique is settled: the species tree's from the
		// column, the molecule tree's from the row. Each branch comes from the other one -- a leaf
		// sits on one branch, whichever leaf of the other tree it is paired with (ADR-0011).
		leaf_pair.prob_z_s_is_one =
		    prob_of_one(column.species_process, row.species_branch, species_parent);
		leaf_pair.prob_z_m_is_one =
		    prob_of_one(row.molecule_process, column.molecule_branch, molecule_parent);

		// The cell data factor, multiplied left to right in the fixed order from 1.0: LOTUS, the
		// passes, mass spectrometry (ADR-0017).
		std::array<double, 2> cell_data{1.0, 1.0};
#ifdef USE_LOTUS
		const TLotus &lotus = _data_model.get_lotus();
		lotus.calculate_LL_update_Y(cell, lotus.holds_a_record(row.lotus_cursor, cell), cell_data);
#endif
		cell_data = _data.cell_factor(row.data, molecule_leaf, cell_data);
#ifdef USE_MS_DATA
		// Exactly {1, 1} for a species without runs, so its draws match a build without mass spec.
		const auto prob_mass_spec =
		    _data_model.get_msms_data().per_cell_factor(row.species_leaf, molecule_leaf);
		cell_data[0] *= prob_mass_spec[0];
		cell_data[1] *= prob_mass_spec[1];
#endif
		leaf_pair.cell_data = {coretools::P(cell_data[0]), coretools::P(cell_data[1])};
		return leaf_pair;
	}

	/// Each source scored both field states before the draw, and only now knows which one to keep.
	/// This writes the row, so the row is not const.
	void record(TRow &row, size_t molecule_leaf,
	            const block_update::TLeafPairFactors & /*factors*/,
	            const field_math::TBlockStates &drawn) const {
#ifdef USE_LOTUS
		// Only a present cell is tallied; an absent one's term does not depend on where it is.
		// Asking the cursor about the cell `factors` asked about moves it nowhere.
		if (drawn.y) {
			const IndexArray cell{row.species_leaf, molecule_leaf};
			const TLotus &lotus = _data_model.get_lotus();
			const auto reported = static_cast<size_t>(lotus.holds_a_record(row.lotus_cursor, cell));
			++row.lotus_present[reported * lotus.n_molecule_classes() +
			                    lotus.molecule_class(molecule_leaf)];
		}
#endif
		_data.record(row.data, molecule_leaf, drawn.y);
	}

	/// The row is done: every pass folds the row's tally into its own, and the row adds its LOTUS
	/// present cells to the counts, one atomic add per kind of cell it holds.
	void end_row(TRow &row) {
		_data.end_row(row.data);
#ifdef USE_LOTUS
		_accumulator.add_lotus_row(row.lotus_species_class, row.lotus_present);
#endif
	}
};
