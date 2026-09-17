//
// What the block update asks about a leaf pair: the two trees' processes, and the data.
//
// The traversal in field/TBlockUpdate.h reads three states and two tree parents from its storages.
// It asks this for the rest. Each tree says what its own tree field cell should be, given that
// tree's parent. Each compiled-in data source says what it makes of the field cell. Every answer
// is about one leaf pair, because that is the unit the traversal works in.
//
// This is the whole of what the traversal must not know about. So the traversal keeps no tree, no
// clique and no data source, and a test can run its loop against a model of its own.
//
// The per-source likelihood bookkeeping lives here too. It is the data's business, and the
// traversal only hands each leaf pair back once.
//

#pragma once

#include "TDataModel.h"
#include "Types.h"
#include "constants.h"
#include "coretools/Math/TSumLog.h"
#include "coretools/Types/probability.h"
#include "field/TBlockUpdate.h"
#include "field/TFieldMath.h"
#include "omp.h"
#include "tree/TTree.h"
#include "tree/branch/TTransitionGrid.h"
#include "tree/branch/TTransitionGridTable.h"
#include <array>
#include <cstddef>
#include <memory>
#include <vector>

//-----------------------------------
// Data-source bookkeeping
//-----------------------------------

/// Per-cell outcome of one block update. Each data source keeps its own likelihood bookkeeping
/// (they are independent terms), so the results are handed back separately instead of merged.
struct TCellOutcome {
#ifdef USE_LOTUS
	/// P(L_cell | Y = the drawn state). Neutral value 1.0 (log 0).
	double prob_lotus_new_state = 1.0;
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
	/// Whether the observed D cell contradicts the state the field was just given.
	bool simple_model_disagrees = false;
#endif
};

/// What one row of a block update owes each data source: the running bookkeeping of one species
/// leaf, kept on the row's own stack and handed over once when the row ends.
///
/// The accumulators are bundled into a single object on purpose: `#ifdef` cannot appear inside a
/// `#pragma omp` line, and a `default(none) shared(...)` clause has to name every variable it
/// touches. One object keeps that clause identical in every build configuration.
struct TRowOutcome {
#ifdef USE_LOTUS
	/// sum over the row's leaf pairs of log P(L_cell | Y = the drawn state).
	coretools::TSumLogProbability lotus_LL;
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
	/// How many of the row's leaf pairs the observed D cell contradicts.
	size_t n_disagree = 0;
#endif

	/// Hot path: called once per updated leaf pair, on the calling row's own object.
	void add([[maybe_unused]] const TCellOutcome &outcome) {
#ifdef USE_LOTUS
		lotus_LL.add(outcome.prob_lotus_new_state);
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
		n_disagree += static_cast<size_t>(outcome.simple_model_disagrees);
#endif
	}
};

/// One slot per field row, committed to the data sources when the update ends.
///
/// One slot per row, and not one per thread. A thread-indexed accumulator adds the rows the
/// schedule handed it, in that order, so the last bits of the likelihood it sums move with the
/// thread count -- and that likelihood is what the gamma and error-rate Metropolis ratios read
/// (`TLotus::ll_ratio_after_parameter_move`). Row slots summed in row order make it a function of
/// the seed alone. `TTree::log_node_state_density` keeps one slot per clique for the same reason.
///
/// A row writes its slot once, when it ends, so two rows sharing a cache line costs one write per
/// 141,619 leaf pairs rather than one per leaf pair.
class TDataUpdateAccumulator {
private:
	std::vector<TRowOutcome> _per_row;

public:
	explicit TDataUpdateAccumulator(size_t n_rows) : _per_row(n_rows) {}

	/// Takes what one row added up to. Called once per row, at the end of the row, from inside
	/// the parallel region. Only ever touches that row's slot.
	void flush_row(size_t species_leaf, const TRowOutcome &row) { _per_row[species_leaf] = row; }

	/// Sums the row slots in row order and installs the results in the data sources. Called once,
	/// after the parallel region.
	void commit([[maybe_unused]] TDataModel &data_model) {
#ifdef USE_LOTUS
		double sum_new_LL = 0.0;
		for (auto &row : _per_row) { sum_new_LL += row.lotus_LL.getSum(); }
		data_model.get_lotus().update_cur_LL(sum_new_LL);
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
		size_t total_disagree = 0;
		for (const auto &row : _per_row) { total_disagree += row.n_disagree; }
		// The update visits every leaf pair exactly once, so this is the complete disagreement
		// count.
		data_model.get_simple_error_model().set_n_disagree(total_disagree);
#endif
	}
};

//-----------------------------------
// The model
//-----------------------------------

/// The two trees and every compiled-in data source, as the block update wants them.
///
/// Both trees have something to say from the first update on, because the chain is given a
/// starting node state before it (TMarkovField::_start_the_chain).
///
/// An inferred chain is the only one that has a model to ask. A simulated one draws its whole
/// configuration forward and runs no update (TMarkovField::simulate).
///
/// Every thread of the update asks this one object. `begin_row` hands back a row that carries
/// everything the walk down that row accumulates, so `factors` and `record` touch the caller's own
/// row and nothing another thread also holds; `end_row` writes that row's slot, once.
class TBlockModel {
private:
	const TTree &_species_tree;
	const TTree &_molecule_tree;
	TDataModel &_data_model;
	TDataUpdateAccumulator &_accumulator;

public:
	TBlockModel(const std::vector<std::unique_ptr<TTree>> &trees, TDataModel &data_model,
	            TDataUpdateAccumulator &accumulator)
	    : _species_tree(*trees.front()), _molecule_tree(*trees.back()), _data_model(data_model),
	      _accumulator(accumulator) {}

	/// Told once, before the parallel region starts, that a traversal is about to begin. A source
	/// whose point query wants the leaf pairs in ascending order readies whatever cannot be
	/// readied concurrently -- LOTUS sorts its record cache here, so that every row may then seed
	/// a cursor of its own. It is told nothing about how the rows are shared out: that is what
	/// used to tie the traversal's schedule to `ProgramOptions::NUMBER_OF_THREADS`. Only LOTUS
	/// uses this today; a build without it is an empty function.
	void prepare_for_traversal() {
#ifdef USE_LOTUS
		_data_model.get_lotus().prepare_for_block_update();
#endif
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
	/// It also carries what the row accumulates: the data sources' running bookkeeping, and
	/// LOTUS's forward cursor into its record cache. Both used to be one slot per thread, indexed
	/// by `omp_get_thread_num()` and written once per leaf pair, which put two threads' slots in
	/// one cache line and made the likelihood's rounding a function of the schedule.
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
		/// What this row has added up for each data source so far.
		TRowOutcome outcome;
#ifdef USE_LOTUS
		/// This row's forward cursor into LOTUS's sorted record cache, seeded at the row's first
		/// cell. The row walks its columns in ascending order, which is all the cursor asks.
		TLotus::RecordCursor lotus_cursor;
#endif
	};

	[[nodiscard]] TRow begin_row(size_t species_leaf, [[maybe_unused]] size_t first_field_cell) {
		// The molecule tree drops the molecule coordinate, so any molecule leaf names this row's
		// clique; 0 is the one every field row has.
		const IndexArray first_cell_of_row{species_leaf, 0};
		return TRow{
		    .molecule_process = _molecule_tree.transition_grid_of_cell(first_cell_of_row),
		    .species_branch   = _species_tree.get_binned_branch_length(species_leaf),
		    .species_leaf     = species_leaf,
		    .molecule_clique  = _molecule_tree.clique_of_cell(first_cell_of_row),
		    .outcome          = {},
#ifdef USE_LOTUS
		    .lotus_cursor = _data_model.get_lotus().cursor_for_row(first_field_cell),
#endif
		};
	}

	[[nodiscard]] block_update::TLeafPairFactors factors(TRow &row, size_t molecule_leaf,
	                                                     bool species_parent,
	                                                     bool molecule_parent) const {
		const IndexArray cell{row.species_leaf, molecule_leaf};
		// What the row settled once has to be what this cell would have asked for.
		DEBUG_ASSERT(_molecule_tree.clique_of_cell(cell) == row.molecule_clique);

		block_update::TLeafPairFactors leaf_pair;
		// The species tree's clique is named by the molecule leaf, so it moves along the row and is
		// read here. The molecule tree's is the species leaf, so it came from the row. Each tree
		// drops the dimension it owns from the cell it is handed, so both take the whole leaf pair
		// and neither caller names a slot to blank (ADR-0011).
		leaf_pair.prob_z_s_is_one = prob_of_one(_species_tree.transition_grid_of_cell(cell),
		                                        row.species_branch, species_parent);
		leaf_pair.prob_z_m_is_one =
		    prob_of_one(row.molecule_process,
			            _molecule_tree.get_binned_branch_length(molecule_leaf), molecule_parent);

		// 1.0 is the neutral value, adding log(1) = 0. A build that left a source out keeps it.
		leaf_pair.lotus        = {coretools::P(1.0), coretools::P(1.0)};
		leaf_pair.simple_error = {coretools::P(1.0), coretools::P(1.0)};
#ifdef USE_LOTUS
		const TLotus &lotus = _data_model.get_lotus();
		std::array<double, 2> prob_lotus{1.0, 1.0};
		lotus.calculate_LL_update_Y(cell, lotus.holds_a_record(row.lotus_cursor, cell),
		                            prob_lotus);
		leaf_pair.lotus = {coretools::P(prob_lotus[0]), coretools::P(prob_lotus[1])};
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
		const TSimpleErrorModel &simple = _data_model.get_simple_error_model();
		std::array<double, 2> prob_simple{};
		simple.probabilities_for_Y_update(simple.observed_state_of(cell), prob_simple);
		leaf_pair.simple_error = {coretools::P(prob_simple[0]), coretools::P(prob_simple[1])};
#endif
		return leaf_pair;
	}

	/// Each source scored both field states before the draw, and only now knows which one to keep.
	/// This writes the row, so the row is not const.
	void record(TRow &row, [[maybe_unused]] size_t molecule_leaf,
	            [[maybe_unused]] const block_update::TLeafPairFactors &factors,
	            [[maybe_unused]] const field_math::TBlockStates &drawn) const {
		TCellOutcome outcome;
#ifdef USE_LOTUS
		outcome.prob_lotus_new_state = factors.lotus[static_cast<size_t>(drawn.y)].get();
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
		const TSimpleErrorModel &simple = _data_model.get_simple_error_model();
		const IndexArray cell{row.species_leaf, molecule_leaf};
		outcome.simple_model_disagrees = simple.observed_state_of(cell) != drawn.y;
#endif
		row.outcome.add(outcome);
	}

	/// The row is done: hand what it added up to the accumulator's slot for that row. One write
	/// per row, so two rows sharing a cache line costs nothing.
	void end_row(TRow &row) { _accumulator.flush_row(row.species_leaf, row.outcome); }
};
