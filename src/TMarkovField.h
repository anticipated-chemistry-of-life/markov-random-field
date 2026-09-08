#pragma once

#include "constants.h"
#include "storages/storage_concepts.h"
#include "tree/TTree.h"

class TDataModel; // forward declaration

/// Per-cell outcome of one Y update. Each data source keeps its own likelihood bookkeeping (they
/// are independent terms), so the results are handed back separately instead of merged.
struct TYUpdateResult {
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

class TMarkovField {
private:
	// trees and Y
	std::vector<std::unique_ptr<TTree>> &_trees;
	TFieldStorage _Y;
	std::string _prefix;

	// fix values?
	bool _fix_Y = false;
	bool _fix_Z = false;

	/// Was Z initialized from children ?
	bool _z_initialized_from_children = false;
	// stuff for updating Y
	IndexArray _num_leaves_per_dim{};

	/// complete joint density of the markov random field
	std::vector<double> _complete_log_density;

	/// Output files
	coretools::TOutputFile _Y_trace_file;
	coretools::TOutputFile _joint_density_file;
	std::vector<coretools::TOutputFile> _Z_trace_files;

private:
	/// TODO: This will change from calculating directly from the trees to a per tree field and the
	/// sum log will just be the AND gate with error.
	void _calculate_log_prob_field(const IndexArray &index_in_leaves_space,
	                               std::array<coretools::TSumLogProbability, 2> &sum_log) const;

#ifdef USE_LOTUS
	static void _calc_lotus_LL(const IndexArray &index_in_leaves_space, std::array<double, 2> &prob,
	                           const TDataModel &data_model);
#endif

	template<bool IsSimulation, bool initYFromData>
	TYUpdateResult _update_Y(const IndexArray &index_in_leaves_space, const TDataModel &data_model,
	                         std::vector<size_t> &linear_indices_in_Y_space_to_insert) {
		auto current_state = _Y.is_one(index_in_leaves_space);
		// prepare log probabilities for the two possible states
		std::array<coretools::TSumLogProbability, 2> sum_log;

		// calculate probabilities in Markov random field
		if constexpr (!initYFromData) { _calculate_log_prob_field(index_in_leaves_space, sum_log); }
		std::array<coretools::TSumLogProbability, 2> sum_log_field = sum_log;

#ifdef USE_LOTUS
		std::array<double, 2> prob_lotus{1.0, 1.0};
#endif
		if constexpr (!IsSimulation) {
			// calculate log likelihood (lotus)
#ifdef USE_LOTUS
			_calc_lotus_LL(index_in_leaves_space, prob_lotus, data_model);
			for (size_t i = 0; i < 2; ++i) { sum_log[i].add(prob_lotus[i]); }
#endif
			// calculate log likelihood (simple error model)
#ifdef USE_SIMPLE_ERROR_MODEL
			std::array<double, 2> prob_simple{};
			_calc_simple_error_model_LL(index_in_leaves_space, prob_simple, data_model);
			for (size_t i = 0; i < 2; ++i) { sum_log[i].add(prob_simple[i]); }
#endif
			// calculate log likelihood mass spec data
			// TODO
		}

		// sample state
		const bool new_state = sample(sum_log);
		_set_new_Y<TFieldStorage>(linear_indices_in_Y_space_to_insert, current_state, new_state);

		// update Y accordingly
		TYUpdateResult result;
#ifdef USE_LOTUS
		result.prob_lotus_new_state = prob_lotus[new_state];
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
		if constexpr (!IsSimulation) {
			result.simple_model_disagrees =
			    _simple_error_model_disagrees(index_in_leaves_space, new_state, data_model);
		}
#endif

		_complete_log_density[omp_get_thread_num()] +=
		    sum_log_field[static_cast<size_t>(new_state)].getSum();

		return result;
	}

#ifdef USE_SIMPLE_ERROR_MODEL
	static void _calc_simple_error_model_LL(const IndexArray &multidim_index,
	                                        std::array<double, 2> &prob,
	                                        const TDataModel &data_model);
	[[nodiscard]] static bool _simple_error_model_disagrees(const IndexArray &multidim_index,
	                                                        bool new_state,
	                                                        const TDataModel &data_model);
#endif

	template<FieldStorage Field>
	void _set_new_Y(std::vector<size_t> &linear_indices_in_Y_space_to_insert,
	                IsOneResult<typename Field::TCell> &current_state, bool new_state) {
		bool cur_state = current_state.is_one;
		if (cur_state && !new_state) {
			current_state.cell->set_state(false);
		} else if (!cur_state && new_state) {
			if (current_state.in_container) {
				current_state.cell->set_state(true);
			} else {
				linear_indices_in_Y_space_to_insert.emplace_back(current_state.linear_index);
			}
		}
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
		std::vector<std::vector<size_t>> linear_indices_in_Y_space_to_insert(
		    ProgramOptions::NUMBER_OF_THREADS);
		TDataSweepAccumulator acc(ProgramOptions::NUMBER_OF_THREADS);

#pragma omp parallel for num_threads(ProgramOptions::NUMBER_OF_THREADS) default(none)              \
    shared(data_model, linear_indices_in_Y_space_to_insert, acc)
		for (size_t i = 0; i < _Y.total_size_of_container_space(); ++i) {
			auto leaf_index = _Y.get_multi_dimensional_index(i);
			auto result     = _update_Y<false, true>(
			    leaf_index, data_model, linear_indices_in_Y_space_to_insert[omp_get_thread_num()]);

			if constexpr (!IsSimulation) {
				acc.add(static_cast<size_t>(omp_get_thread_num()), result);
			}
		}

		_Y.insert_ones_in_container(linear_indices_in_Y_space_to_insert);
		if constexpr (!IsSimulation) { acc.commit(data_model); }
	}

	void _reset_log_joint_density() {
		_complete_log_density.clear();
		_complete_log_density.resize(ProgramOptions::NUMBER_OF_THREADS);
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
			line                 = {i, _Y.is_one(i).is_one};
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
		for (size_t i = 0; i < _Y.total_size_of_container_space(); ++i) {
			// A cell that is not a one now and was never counted a one carries no posterior. Which
			// cells are *stored* is a property of the backend -- the sparse field holds the cells
			// it was given, ones and zeros alike, and the dense one holds the whole container
			// space -- so leaving those cells out is what makes this file say the same thing under
			// either backend, and it drops no information: every column of such a row is the
			// default.

			auto storage = _Y.is_one(i);
			if (!storage.is_one && storage.cell->get_counter() == 0) { continue; }
			line                 = {storage.linear_index, storage.is_one};
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
	double _calculate_complete_joint_density();

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

		for (auto &_tree : _trees) {
			_tree->update_Z_and_nus_and_alphas_and_branch_lengths<IsSimulation, FixZ>(_Y);
		}
		if (_fix_Z) { return; }
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
	[[nodiscard]] const TFieldStorage &get_Y_matrix() const;

	// functions to perform stuff on Y after burnin / MCMC finished
	void burninHasFinished();
	void MCMCHasFinished();
	void oneBurninHasFinished();

	static size_t get_num_iterations_simulation() { return ProgramOptions::NUM_ITERATIONS; }
};
