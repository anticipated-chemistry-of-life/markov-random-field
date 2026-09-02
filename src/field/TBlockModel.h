//
// What the block update asks about a leaf pair: the two trees' processes, and the data.
//
// The traversal in field/TBlockUpdate.h reads three states and two tree parents through windows of
// its own. It asks this for the rest. Each tree says what its own tree field cell should be, given
// that tree's parent. Each compiled-in data source says what it makes of the field cell.
//
// This is the whole of what the traversal must not know about. So the traversal keeps no tree, no
// clique and no data source, and a test can run its loop against a model of its own.
//
// The per-source likelihood bookkeeping lives here too. It is the data's business, and the
// traversal only hands each leaf pair back once.
//

#pragma once

#include "TClique.h"
#include "TDataModel.h"
#include "Types.h"
#include "constants.h"
#include "coretools/Math/TSumLog.h"
#include "coretools/Types/probability.h"
#include "field/TBlockUpdate.h"
#include "field/TFieldMath.h"
#include "omp.h"
#include "storages/storage_backend.h"
#include "tree/TTree.h"
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

/// Per-thread accumulators for one full block update, committed to the data sources at the end.
///
/// The accumulators are bundled into a single object on purpose: `#ifdef` cannot appear inside a
/// `#pragma omp` line, and a `default(none) shared(...)` clause has to name every variable it
/// touches. One object keeps that clause identical in every build configuration.
class TDataUpdateAccumulator {
private:
#ifdef USE_LOTUS
	std::vector<coretools::TSumLogProbability> _lotus_LL;
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
	std::vector<size_t> _n_disagree;
#endif

public:
	/// Sizing happens in the body rather than in a member-initializer list, so that adding or
	/// removing a source does not require rebalancing the commas of a #ifdef'd init list.
	explicit TDataUpdateAccumulator([[maybe_unused]] size_t n_threads) {
#ifdef USE_LOTUS
		_lotus_LL.resize(n_threads);
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
		_n_disagree.assign(n_threads, 0);
#endif
	}

	/// Hot path: called once per updated leaf pair, from inside the parallel region. Only ever
	/// touches the slot of the calling thread.
	void add([[maybe_unused]] size_t thread, [[maybe_unused]] const TCellOutcome &outcome) {
#ifdef USE_LOTUS
		_lotus_LL[thread].add(outcome.prob_lotus_new_state);
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
		_n_disagree[thread] += static_cast<size_t>(outcome.simple_model_disagrees);
#endif
	}

	/// Sums the per-thread slots and installs the results in the data sources. Called once, after
	/// the parallel region.
	void commit([[maybe_unused]] TDataModel &data_model) {
#ifdef USE_LOTUS
		double sum_new_LL = 0.0;
		for (auto &i : _lotus_LL) { sum_new_LL += i.getSum(); }
		data_model.get_lotus().update_cur_LL(sum_new_LL);
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
		size_t total_disagree = 0;
		for (const auto &i : _n_disagree) { total_disagree += i; }
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
/// @tparam IsSimulation   a simulated chain draws from the prior, so every data term is neutral.
/// @tparam InitYFromData  the chain's first update runs before any internal node has a state, so
///                        neither tree has anything to say yet. Both tree factors become a coin
///                        flip, and the data and the link decide the leaf pair on their own. The
///                        node states are built from the tree fields this leaves behind
///                        (TTree::initialize_Z_from_children).
template<bool IsSimulation, bool InitYFromData> class TBlockModel {
private:
	const TTree &_species_tree;
	const TTree &_molecule_tree;
	TDataModel &_data_model;
	TDataUpdateAccumulator &_accumulator;

public:
	/// One species leaf's row: the data sources' cells for it, and the two things that are the
	/// same for every leaf pair in it -- the branch down to the species leaf, and the molecule
	/// tree's clique, which a species leaf names.
	class TRow {
	private:
		TBlockModel *_model;
		size_t _species_leaf;
		size_t _thread;
		TypeBinnedBranchLengths _species_branch;
		const TClique *_molecule_clique;
#ifdef USE_LOTUS
		TFieldStorage::TWindow _lotus_row;
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
		TFieldStorage::TWindow _simple_data_row;
#endif

	public:
		/// Both data sources have the field's dimensions, so the field's index is already theirs.
		TRow(TBlockModel &model, size_t species_leaf)
		    : _model(&model), _species_leaf(species_leaf),
		      _thread(static_cast<size_t>(omp_get_thread_num())),
		      _species_branch(model._species_tree.get_binned_branch_length(species_leaf)),
		      _molecule_clique(&model._molecule_tree.get_clique(IndexArray{species_leaf, 0}))
#ifdef USE_LOTUS
		      ,
		      _lotus_row(model._data_model.get_lotus().open_row(
		          IndexArray{species_leaf, 0}, model._molecule_tree.get_number_of_leaves()))
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
		      ,
		      _simple_data_row(model._data_model.get_simple_error_model().open_row(
		          IndexArray{species_leaf, 0}, model._molecule_tree.get_number_of_leaves()))
#endif
		{
		}

		/// P(a tree field cell = 1 | its parent's state), under one clique's process on one branch.
		[[nodiscard]] static coretools::Probability
		prob_of_one(const TClique &clique, TypeBinnedBranchLengths branch, bool parent_state) {
			return coretools::P(clique.transition_grid().probability(branch, parent_state,
			                                                         /*to=*/true));
		}

		[[nodiscard]] block_update::TLeafPairFactors
		factors(size_t molecule_leaf, bool species_parent, bool molecule_parent) const {
			const IndexArray cell{_species_leaf, molecule_leaf};

			block_update::TLeafPairFactors leaf_pair;
			if constexpr (InitYFromData) {
				leaf_pair.prob_z_s_is_one = coretools::P(0.5);
				leaf_pair.prob_z_m_is_one = coretools::P(0.5);
			} else {
				// The species tree's clique is named by the molecule leaf, so it changes along the
				// row. The molecule tree's is named by the species leaf, so it does not.
				leaf_pair.prob_z_s_is_one = prob_of_one(_model->_species_tree.get_clique(cell),
				                                        _species_branch, species_parent);
				leaf_pair.prob_z_m_is_one = prob_of_one(
				    *_molecule_clique, _model->_molecule_tree.get_binned_branch_length(molecule_leaf),
				    molecule_parent);
			}

			// 1.0 is the neutral value, adding log(1) = 0. A simulated chain keeps it, and so does
			// a build that left a source out.
			leaf_pair.lotus        = {coretools::P(1.0), coretools::P(1.0)};
			leaf_pair.simple_error = {coretools::P(1.0), coretools::P(1.0)};
			if constexpr (!IsSimulation) {
#ifdef USE_LOTUS
				std::array<double, 2> prob_lotus{1.0, 1.0};
				_model->_data_model.get_lotus().calculate_LL_update_Y(
				    cell, _lotus_row.is_one(molecule_leaf), prob_lotus);
				leaf_pair.lotus = {coretools::P(prob_lotus[0]), coretools::P(prob_lotus[1])};
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
				std::array<double, 2> prob_simple{};
				_model->_data_model.get_simple_error_model().probabilities_for_Y_update(
				    _simple_data_row.is_one(molecule_leaf), prob_simple);
				leaf_pair.simple_error = {coretools::P(prob_simple[0]),
				                          coretools::P(prob_simple[1])};
#endif
			}
			return leaf_pair;
		}

		/// Each source scored both field states before the draw, and only now knows which one to
		/// keep. This writes the accumulator, so it is not const.
		void record([[maybe_unused]] size_t molecule_leaf,
		            [[maybe_unused]] const block_update::TLeafPairFactors &factors,
		            [[maybe_unused]] const field_math::TBlockStates &drawn) {
			if constexpr (!IsSimulation) {
				TCellOutcome outcome;
#ifdef USE_LOTUS
				outcome.prob_lotus_new_state = factors.lotus[static_cast<size_t>(drawn.y)].get();
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
				outcome.simple_model_disagrees = _simple_data_row.is_one(molecule_leaf) != drawn.y;
#endif
				_model->_accumulator.add(_thread, outcome);
			}
		}
	};

	TBlockModel(const std::vector<std::unique_ptr<TTree>> &trees, TDataModel &data_model,
	            TDataUpdateAccumulator &accumulator)
	    : _species_tree(*trees.front()), _molecule_tree(*trees.back()), _data_model(data_model),
	      _accumulator(accumulator) {}

	/// A row is returned as a prvalue, and never moved: it holds windows, and a window owns writes
	/// that are not in its storage yet.
	[[nodiscard]] TRow open_row(size_t species_leaf) { return TRow(*this, species_leaf); }
};
