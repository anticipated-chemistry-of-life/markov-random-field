//
// Created by madleina on 10.01.25.
//

#ifndef ACOL_TMARKOVFIELD_H
#define ACOL_TMARKOVFIELD_H

#include "TClique.h"
#include "TCurrentState.h"

#include "Types.h"
#include "cli.h"
#include "constants.h"
#include "coretools/Files/TOutputFile.h"
#include "coretools/Main/TError.h"
#include "coretools/algorithms.h"
#include "mass_spec/msms_data.h"
#include "omp.h"
#include "tree/TTree.h"
#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

//-----------------------------------
// Y sweep bookkeeping
//-----------------------------------

class TDataModel; // forward declaration

/// Per-cell outcome of one Y update. Each data source keeps its own likelihood bookkeeping (they
/// are independent terms), so the results are handed back separately instead of merged.
struct TYUpdateResult {
	int diff_counter_1_in_last_dim = 0;
#ifdef USE_LOTUS
	/// P(L_cell | x = new_state). Neutral value 1.0 (log 0), used when collapsing makes the LOTUS
	/// term identical for both Y states and calculate_LL_update_Y leaves it untouched.
	double prob_lotus_new_state = 1.0;
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
	/// Whether the observed D cell contradicts the state Y was just set to.
	bool simple_model_disagrees = false;
#endif
};

/// Per-thread accumulators for one full Y sweep, committed to the data sources at the end.
///
/// The accumulators are bundled into a single object on purpose: `#ifdef` cannot appear inside a
/// `#pragma omp` line, and the sweep's `default(none) shared(...)` clause has to name every
/// variable it touches. One object keeps that clause identical in every build configuration.
class TDataSweepAccumulator {
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
	explicit TDataSweepAccumulator(size_t n_threads) {
#ifdef USE_LOTUS
		_lotus_LL.resize(n_threads);
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
		_n_disagree.assign(n_threads, 0);
#endif
	}

	/// Hot path: called once per updated Y cell, from inside the parallel region. Only ever touches
	/// the slot of the calling thread.
	void add(size_t thread, const TYUpdateResult &result) {
#ifdef USE_LOTUS
		_lotus_LL[thread].add(result.prob_lotus_new_state);
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
		_n_disagree[thread] += static_cast<size_t>(result.simple_model_disagrees);
#endif
	}

	/// Sums the per-thread slots and installs the results in the data sources. Called once, after
	/// the parallel region.
	void commit(TDataModel &data_model);
};

//-----------------------------------
// TMarkovField
//-----------------------------------

class TMarkovField {
private:
	// trees and Y
	std::vector<std::unique_ptr<TTree>> &_trees;
	TStorageYMatrix _Y;
	std::string _prefix;

	// stuff for updating Y
	size_t _K;
	size_t _num_outer_loops;
	IndexArray _num_leaves_per_dim_except_last{};
	std::vector<TSheet> _sheets;
	TCurrentState _clique_last_dim;

	// fix values?
	bool _fix_Y = false;
	bool _fix_Z = false;

	// mass spectrometry data and the dimension indices of molecules/species within _trees
	std::optional<TMSMSData> _ms_data;

	// complete joint density of the markov random field
	std::vector<double> _complete_log_density;

	/// Was Z initialized from children ?
	bool _z_initialized_from_children = false;

	// output files
	coretools::TOutputFile _Y_trace_file;
	std::vector<coretools::TOutputFile> _Z_trace_files;
	coretools::TOutputFile _joint_density_file;

	// functions for updating Y
	void _update_sheets(bool first, IndexArray &start_index_in_leaves_space,
	                    IndexArray &previous_ix, size_t K_cur_sheet);
	void _fill_clique_along_last_dim(IndexArray start_index_in_leaves_space);
	void _calculate_log_prob_field(const IndexArray &index_in_leaves_space,
	                               std::array<coretools::TSumLogProbability, 2> &sum_log) const;
	[[nodiscard]] bool _need_to_update_sheet(size_t sheet_ix,
	                                         const IndexArray &start_index_in_leaves_space,
	                                         const IndexArray &previous_ix) const;
	int _set_new_Y(bool new_state, const IndexArray &index_in_leaves_space,
	               std::vector<size_t> &linear_indices_in_Y_space_to_insert);
	void _update_counter_1_cliques(bool new_state, bool old_state,
	                               const IndexArray &index_in_leaves_space);

	void _simulate_Y();
	// These forward into the data sources. They are declared here and defined in TMarkovField.cpp
	// (which includes TDataModel.h) because TDataModel is only forward-declared in this header:
	// _update_Y is a template, but `data_model` is not a dependent type, so any member access on it
	// would be checked right here against an incomplete type.
#ifdef USE_LOTUS
	void _calc_lotus_LL(const IndexArray &index_in_leaves_space, size_t index_for_tmp_state,
	                    size_t leaf_index_last_dim, std::array<double, 2> &prob,
	                    const TDataModel &data_model);
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
	static void _calc_simple_error_model_LL(size_t index_for_tmp_state, std::array<double, 2> &prob,
	                                        const TDataModel &data_model);
	[[nodiscard]] static bool _simple_error_model_disagrees(size_t index_for_tmp_state,
	                                                        bool new_state,
	                                                        const TDataModel &data_model);
#endif
	/// Per-sheet preparation: every compiled-in source caches the slice of its data that the sweep
	/// is about to walk over.
	static void _prepare_data_LL(const IndexArray &start_index_in_leaves_space, size_t K_cur_sheet,
	                             TDataModel &data_model);
	double _calculate_complete_joint_density();
	void _reset_log_joint_density() {
		_complete_log_density.clear();
		_complete_log_density.resize(ProgramOptions::NUMBER_OF_THREADS);
	}

	template<bool IsSimulation, bool initYFromData>
	TYUpdateResult _update_Y(const IndexArray &index_in_leaves_space, size_t leaf_index_last_dim,
	                         size_t index_for_tmp_state,
	                         std::vector<size_t> &linear_indices_in_Y_space_to_insert,
	                         const TDataModel &data_model) {
		auto index_copy   = index_in_leaves_space;
		index_copy.back() = leaf_index_last_dim;

		// prepare log probabilities for the two possible states
		std::array<coretools::TSumLogProbability, 2> sum_log;

		// calculate probabilities in Markov random field
		if constexpr (!initYFromData) { _calculate_log_prob_field(index_copy, sum_log); }
		std::array<coretools::TSumLogProbability, 2> sum_log_field = sum_log;

		// Declared outside the IsSimulation branch so the simulation instantiation does not warn
		// about an unused variable. 1.0 is the neutral value: calculate_LL_update_Y leaves prob
		// untouched when collapsing makes the LOTUS term identical for both states, and adding
		// log(1) = 0 is exactly the no-op that case needs.
#ifdef USE_LOTUS
		std::array<double, 2> prob_lotus{1.0, 1.0};
#endif
		if constexpr (!IsSimulation) {
			// calculate log likelihood (lotus)
#ifdef USE_LOTUS
			_calc_lotus_LL(index_copy, index_for_tmp_state, leaf_index_last_dim, prob_lotus,
			               data_model);
			for (size_t i = 0; i < 2; ++i) { sum_log[i].add(prob_lotus[i]); }
#endif
			// calculate log likelihood (simple error model)
#ifdef USE_SIMPLE_ERROR_MODEL
			std::array<double, 2> prob_simple{};
			_calc_simple_error_model_LL(index_for_tmp_state, prob_simple, data_model);
			for (size_t i = 0; i < 2; ++i) { sum_log[i].add(prob_simple[i]); }
#endif
			// calculate log likelihood mass spec data
			if (_ms_data.has_value()) { _ms_data->add_log_likelihood(index_copy, sum_log); }
		}

		// sample state
		const bool new_state = sample(sum_log);

		// update Y accordingly
		TYUpdateResult result;
		result.diff_counter_1_in_last_dim =
		    _set_new_Y(new_state, index_copy, linear_indices_in_Y_space_to_insert);
#ifdef USE_LOTUS
		result.prob_lotus_new_state = prob_lotus[new_state];
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
		if constexpr (!IsSimulation) {
			result.simple_model_disagrees =
			    _simple_error_model_disagrees(index_for_tmp_state, new_state, data_model);
		}
#endif

		_complete_log_density[omp_get_thread_num()] +=
		    sum_log_field[static_cast<size_t>(new_state)].getSum();

		return result;
	}

	template<bool IsSimulation, bool initYFromData>
	void _update_all_Y(TDataModel &data_model, size_t iteration) {
		_reset_log_joint_density();

		if (iteration == 0 && ProgramOptions::WRITE_Y_TRACE && !_Y_trace_file.isOpen() && !_fix_Y) {
			std::vector<size_t> Y_trace_header;
			Y_trace_header.reserve(_Y.total_size_of_container_space());
			for (size_t i = 0; i < _Y.total_size_of_container_space(); ++i) {
				Y_trace_header.push_back(i);
			}
			if constexpr (IsSimulation) {
				_Y_trace_file.open(_prefix + "_simulated_Y_trace.txt", Y_trace_header, "\t");
			} else {
				_Y_trace_file.open(_prefix + "_Y_trace.txt", Y_trace_header, "\t");
			}
		}

		if (_fix_Y) {
			// keep the two ifs separate because if Y is not empty, then we just return
			if (_Y.empty()) {
				throw coretools::TUserError(
				    "Y is currently empty and fixed. Was Y read from a file ? "
				    "(--set_Y)");
			}
			return;
		}

		// loop over sheets in last dimension
		TDataSweepAccumulator acc(ProgramOptions::NUMBER_OF_THREADS);
		std::vector<std::vector<size_t>> linear_indices_in_Y_space_to_insert(
		    ProgramOptions::NUMBER_OF_THREADS);

		// Persistent thread team for the whole sweep: the team is created ONCE here instead of once
		// per inner iteration (the old `omp parallel for` sat inside the k x i loop, paying a
		// fork/join every inner iteration). The k/i loops are now executed redundantly by all
		// threads (SPMD) and the work is shared via `omp for`/`omp single`, turning the per-inner
		// fork/joins into cheap barriers on a warm team.
		// previous_ix and diff_counter_1_in_last_dim are shared across the team: previous_ix is
		// read by all threads in _need_to_update_sheet and written in the post `single`; the
		// reduction combines into diff_counter_1_in_last_dim (reset to 0 in the prep `single` each
		// iteration).
		IndexArray previous_ix;
		int diff_counter_1_in_last_dim = 0;
#pragma omp parallel num_threads(ProgramOptions::NUMBER_OF_THREADS) default(none)                  \
    shared(acc, linear_indices_in_Y_space_to_insert, previous_ix, diff_counter_1_in_last_dim,      \
               data_model)
		{
			for (size_t k = 0; k < _num_outer_loops; ++k) {
				const size_t start_ix_in_leaves_last_dim = k * _K; // 0, _K, 2*_K, ...

				// loop over all dimensions except last (linearized)
				const size_t num_inner_loops =
				    coretools::containerProduct(_num_leaves_per_dim_except_last);
				for (size_t i = 0; i < num_inner_loops; ++i) {
					// get multi-dimensional index from linear coordinate and set the start of the
					// last dimension
					auto start_index_in_leaves_space =
					    coretools::getSubscriptsAsArray(i, _num_leaves_per_dim_except_last);
					start_index_in_leaves_space.back() = start_ix_in_leaves_last_dim;
					// calculate size of current sheet (make sure not to overshoot)
					const size_t K_cur_sheet = std::min(_K, _trees.back()->get_number_of_leaves() -
					                                            start_ix_in_leaves_last_dim);
					// update sheet(s), if necessary. Called by ALL threads: TSheet::fill uses
					// worksharing (omp for) and thus distributes over this team.
					_update_sheets(i == 0, start_index_in_leaves_space, previous_ix, K_cur_sheet);

					// serial prep that writes shared state, done by one thread (implicit barrier)
#pragma omp single
					{
						// fill clique along last dimension
						_fill_clique_along_last_dim(start_index_in_leaves_space);
						_prepare_data_LL(start_index_in_leaves_space, K_cur_sheet, data_model);
						diff_counter_1_in_last_dim = 0; // reset before the reduction below
					}

					// now loop along all leaves of the last dimension for updating (only K leaves
					// for which we have everything)
					const size_t end_ix_in_leaves_last_dim =
					    start_ix_in_leaves_last_dim + K_cur_sheet;
#pragma omp for schedule(static) reduction(+ : diff_counter_1_in_last_dim)
					for (size_t j = start_ix_in_leaves_last_dim; j < end_ix_in_leaves_last_dim;
					     ++j) {
						const auto result = _update_Y<IsSimulation, initYFromData>(
						    start_index_in_leaves_space, j, j - start_ix_in_leaves_last_dim,
						    linear_indices_in_Y_space_to_insert[omp_get_thread_num()], data_model);
						diff_counter_1_in_last_dim += result.diff_counter_1_in_last_dim;
						if constexpr (!IsSimulation) {
							acc.add(static_cast<size_t>(omp_get_thread_num()), result);
						}
					}

					// insert new 1-valued indices into Y
					// Note: indices of where Y is one in _sheets is not accurate anymore, but we
					// don't use them, so it's ok
#pragma omp single
					{
						_trees.back()
						    ->get_clique(start_index_in_leaves_space)
						    .update_counter_leaves_state_1(diff_counter_1_in_last_dim);
						previous_ix = start_index_in_leaves_space;
					}
				}
			}
		}

		_Y.insert_in_Y(linear_indices_in_Y_space_to_insert);
		// at the very end: sum the per-thread accumulators and store them in the data sources
		if constexpr (!IsSimulation) { acc.commit(data_model); }
		if (ProgramOptions::WRITE_Y_TRACE && (iteration % _Y.get_thinning_factor() == 0) &&
		    !_fix_Y) {
			_Y_trace_file.writeln(_Y.get_full_Y_binary_vector());
		}
	}

	void _read_Y_from_file(const std::string &filename);

	template<bool IsSimulation, bool FixZ> void _update_all_Z(size_t iteration) {
		if (iteration == 0 && ProgramOptions::WRITE_Z_TRACE && _Z_trace_files.empty() && !_fix_Z) {
			for (const auto &tree : _trees) {
				std::vector<size_t> Z_trace_header;
				Z_trace_header.reserve(tree->get_Z().total_size_of_container_space());
				for (size_t i = 0; i < tree->get_Z().total_size_of_container_space(); ++i) {
					Z_trace_header.push_back(i);
				}
				_Z_trace_files.emplace_back(_prefix + "_" + tree->get_tree_name() + "_Z_trace.txt",
				                            Z_trace_header, "\t");
			}
		}

		// Sequential on purpose: each tree's nu/alpha update reads the *other* trees' Z, alpha, nu
		// and branch lengths (for the leaf terms of its pseudo-likelihood), so they must be frozen
		// while it runs. Parallelism lives inside each call, over that tree's cliques.
		for (auto &_tree : _trees) {
			_tree->update_Z_and_nus_and_alphas_and_branch_lengths<IsSimulation, FixZ>(_Y, _trees);
		}
		if (_fix_Z) { return; }
		if (iteration % _Y.get_thinning_factor() == 0 && ProgramOptions::WRITE_Z_TRACE) {
			for (size_t tree_idx = 0; tree_idx < _trees.size(); ++tree_idx) {
				const auto &tree = _trees[tree_idx];
				_Z_trace_files[tree_idx].writeln(tree->get_Z().get_full_Z_binary_vector());
			}
		}
	}

	template<bool WriteFullY> void _write_Y_to_file(const std::string &filename) const {
		std::vector<std::string> header;
		header.emplace_back("position");
		header.emplace_back("Y_state");
		for (const auto &tree : _trees) { header.push_back(tree->get_tree_name()); }
		header.emplace_back("fraction_of_one");

		coretools::TOutputFile file(filename, header, "\t");
		if constexpr (WriteFullY) {
			this->_write_full_Y(file);
		} else {
			this->_write_only_values_in_Y_vector(file);
		}
	}

	void _write_full_Y(coretools::TOutputFile &file) const {
		std::array<size_t, 2> line{};
		for (size_t i = 0; i < _Y.total_size_of_container_space(); ++i) {
			// a missing cell reads as state 0, so a direct point lookup covers both cases
			line                 = {i, _Y.is_one(i)};
			auto leaf_index_of_Y = _Y.get_multi_dimensional_index(i);
			std::vector<std::string> node_names;
			for (size_t idx = 0; idx < leaf_index_of_Y.size(); ++idx) {
				size_t node_idx = _trees[idx]->get_node_index_from_leaf_index(leaf_index_of_Y[idx]);
				node_names.push_back(_trees[idx]->get_node_id(node_idx));
			};
			const double fraction = _Y.get_fraction_of_ones(i);
			file.writeln(line, node_names, fraction);
		}
	}

	void _write_only_values_in_Y_vector(coretools::TOutputFile &file) const {
		std::array<size_t, 2> line{};
		// iterate only the stored (non-default) cells, in ascending linear-index order
		for (const auto &[linear_index_in_y, storage] : _Y.get_stored_entries()) {
			line                 = {linear_index_in_y, storage.is_one()};
			auto leaf_index_of_Y = _Y.get_multi_dimensional_index(linear_index_in_y);
			std::vector<std::string> node_names;
			for (size_t idx = 0; idx < leaf_index_of_Y.size(); ++idx) {
				size_t node_idx = _trees[idx]->get_node_index_from_leaf_index(leaf_index_of_Y[idx]);
				node_names.push_back(_trees[idx]->get_node_id(node_idx));
			};
			const double fraction = _Y.get_fraction_of_ones(linear_index_in_y);
			file.writeln(line, node_names, fraction);
		}
	}

public:
	TMarkovField(size_t n_iterations, std::vector<std::unique_ptr<TTree>> &Trees,
	             std::string _prefix);
	~TMarkovField() = default;

	// updates
	void update(TDataModel &data_model, size_t iteration);

	// simulation
	void simulate(TDataModel &data_model);

	// get Y
	[[nodiscard]] const TStorageYMatrix &get_Y_matrix() const;

	// functions to perform stuff on Y after burnin / MCMC finished
	void burninHasFinished();
	void MCMCHasFinished();
	void oneBurninHasFinished();

	static size_t get_num_iterations_simulation() {
		return coretools::instances::parameters().get("num_iterations", 5000);
	}
};

#endif // ACOL_TMARKOVFIELD_H
