//
// Created by madleina on 10.01.25.
//

#ifndef ACOL_TMARKOVFIELD_H
#define ACOL_TMARKOVFIELD_H

#include "TClique.h"
#include "TCurrentState.h"
#include "TTree.h"
#include "Types.h"
#include "coretools/Files/TOutputFile.h"
#include "coretools/Main/TError.h"
#include "coretools/devtools.h"
#include "omp.h"
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

//-----------------------------------
// TMarkovField
//-----------------------------------

class TLotus; // forward declaration

class TMarkovField {
private:
	// trees and Y
	std::vector<std::unique_ptr<TTree>> &_trees;
	TStorageYVector _Y;
	std::string _prefix;

	// stuff for updating Y
	size_t _K;
	size_t _num_outer_loops;
	std::vector<size_t> _num_leaves_per_dim_except_last;
	std::vector<TSheet> _sheets;
	TCurrentState _clique_last_dim;

	// fix values?
	bool _fix_Y = false;
	bool _fix_Z = false;

	// complete joint density of the markov random field
	std::vector<double> _complete_log_density;

	// output files
	coretools::TOutputFile _Y_trace_file;
	std::vector<coretools::TOutputFile> _Z_trace_files;
	coretools::TOutputFile _joint_density_file;

	// functions for updating Y
	void _update_sheets(bool first, const std::vector<size_t> &start_index_in_leaves_space,
	                    const std::vector<size_t> &previous_ix, size_t K_cur_sheet);
	void _fill_clique_along_last_dim(std::vector<size_t> start_index_in_leaves_space);
	void _calculate_log_prob_field(const std::vector<size_t> &index_in_leaves_space,
	                               std::array<coretools::TSumLogProbability, 2> &sum_log) const;
	bool _need_to_update_sheet(size_t sheet_ix, const std::vector<size_t> &start_index_in_leaves_space,
	                           const std::vector<size_t> &previous_ix) const;
	int _set_new_Y(bool new_state, const std::vector<size_t> &index_in_leaves_space,
	               std::vector<TStorageY> &linear_indices_in_Y_space_to_insert);
	void _update_counter_1_cliques(bool new_state, bool old_state, const std::vector<size_t> &index_in_leaves_space);

	void _simulate_Y();
	void _calc_lotus_LL(const std::vector<size_t> &index_in_leaves_space, size_t index_for_tmp_state,
	                    size_t leaf_index_last_dim, std::array<double, 2> &prob, TLotus &lotus);
	static void _prepare_lotus_LL(const std::vector<size_t> &start_index_in_leaves_space, size_t K_cur_sheet,
	                              TLotus &lotus);
	static void _update_cur_LL_lotus(TLotus &lotus, std::vector<coretools::TSumLogProbability> &new_LL);
	double _calculate_complete_joint_density();
	void _reset_log_joint_density() {
		_complete_log_density.clear();
		_complete_log_density.resize(NUMBER_OF_THREADS);
	}

	template<bool IsSimulation>
	std::pair<int, double> _update_Y(std::vector<size_t> index_in_leaves_space, size_t leaf_index_last_dim,
	                                 size_t index_for_tmp_state,
	                                 std::vector<TStorageY> &linear_indices_in_Y_space_to_insert, TLotus &lotus) {
		index_in_leaves_space.back() = leaf_index_last_dim;

		// prepare log probabilities for the two possible states
		std::array<coretools::TSumLogProbability, 2> sum_log;

		// calculate probabilities in Markov random field
		_calculate_log_prob_field(index_in_leaves_space, sum_log);
		std::array<coretools::TSumLogProbability, 2> sum_log_field = sum_log;

		// calculate log likelihood (lotus)
		std::array<double, 2> prob_lotus{1e-20, 1e-20};
		if constexpr (!IsSimulation) {
			_calc_lotus_LL(index_in_leaves_space, index_for_tmp_state, leaf_index_last_dim, prob_lotus, lotus);
			for (size_t i = 0; i < 2; ++i) { sum_log[i].add(prob_lotus[i]); }
		}

		// calculate log likelihood (virtual mass spec)...

		// sample state
		bool new_state = sample(sum_log);

		// update Y accordingly
		int diff_counter_1_in_last_dim =
		    _set_new_Y(new_state, index_in_leaves_space, linear_indices_in_Y_space_to_insert);
		double prob_new_state = prob_lotus[new_state];

		if (new_state) {
			_complete_log_density[omp_get_thread_num()] += sum_log_field[1].getSum();
		} else {
			_complete_log_density[omp_get_thread_num()] += sum_log_field[0].getSum();
		}

		return {diff_counter_1_in_last_dim, prob_new_state};
	}

	template<bool IsSimulation> void _update_all_Y(TLotus &lotus, size_t iteration) {
		_reset_log_joint_density();

		if (iteration == 0 && WRITE_Y_TRACE && !_Y_trace_file.isOpen() && !_fix_Y) {
			std::vector<size_t> Y_trace_header;
			for (size_t i = 0; i < _Y.total_size_of_container_space(); ++i) { Y_trace_header.push_back(i); }
			if constexpr (IsSimulation) {
				_Y_trace_file.open(_prefix + "_simulated_Y_trace.txt", Y_trace_header, "\t");
			} else {
				_Y_trace_file.open(_prefix + "_Y_trace.txt", Y_trace_header, "\t");
			}
		}

		if (_fix_Y) {
			if (_Y.empty()) {
				UERROR("Y is currently empty and fixed. Was Y read from a file ? "
				       "(--set_Y)");
			}
			return;
		}

		// loop over sheets in last dimension
		std::vector<coretools::TSumLogProbability> new_LL(NUMBER_OF_THREADS);
		std::vector<std::vector<TStorageY>> linear_indices_in_Y_space_to_insert(NUMBER_OF_THREADS);
		for (size_t k = 0; k < _num_outer_loops; ++k) {
			const size_t start_ix_in_leaves_last_dim = k * _K; // 0, _K, 2*_K, ...

			// loop over all dimensions except last (linearized)
			size_t num_inner_loops = coretools::containerProduct(_num_leaves_per_dim_except_last);
			std::vector<size_t> previous_ix;
			for (size_t i = 0; i < num_inner_loops; ++i) {
				// get multi-dimensional index from linear coordinate and set the start of the last dimension
				auto start_index_in_leaves_space   = coretools::getSubscripts(i, _num_leaves_per_dim_except_last);
				start_index_in_leaves_space.back() = start_ix_in_leaves_last_dim;
				// calculate size of current sheet (make sure not to overshoot)
				const size_t K_cur_sheet =
				    std::min(_K, _trees.back()->get_number_of_leaves() - start_ix_in_leaves_last_dim);
				// update sheet(s), if necessary
				_update_sheets(i == 0, start_index_in_leaves_space, previous_ix, K_cur_sheet);

				// fill clique along last dimension
				_fill_clique_along_last_dim(start_index_in_leaves_space);
				_prepare_lotus_LL(start_index_in_leaves_space, K_cur_sheet, lotus);

				// now loop along all leaves of the last dimension for updating (only K leaves for which we have
				// everything)
				const size_t end_ix_in_leaves_last_dim = start_ix_in_leaves_last_dim + K_cur_sheet;
				int diff_counter_1_in_last_dim         = 0;
#pragma omp parallel for num_threads(NUMBER_OF_THREADS) reduction(+ : diff_counter_1_in_last_dim)
				for (size_t j = start_ix_in_leaves_last_dim; j < end_ix_in_leaves_last_dim; ++j) {
					auto [diff, prob_new_state] =
					    _update_Y<IsSimulation>(start_index_in_leaves_space, j, j - start_ix_in_leaves_last_dim,
					                            linear_indices_in_Y_space_to_insert[omp_get_thread_num()], lotus);
					diff_counter_1_in_last_dim += diff;
					if constexpr (!IsSimulation) { new_LL[omp_get_thread_num()].add(prob_new_state); }
				}

				// insert new 1-valued indices into Y
				// Note: indices of where Y is one in _sheets is not accurate anymore, but we don't use them, so it's
				// ok
				_trees.back()
				    ->get_clique(start_index_in_leaves_space)
				    .update_counter_leaves_state_1(diff_counter_1_in_last_dim);

				previous_ix = start_index_in_leaves_space;
			}
		}

		_Y.insert_in_Y(linear_indices_in_Y_space_to_insert);
		// at the very end: sum the LL of all threads and store it in TLotus
		if constexpr (!IsSimulation) { _update_cur_LL_lotus(lotus, new_LL); }
		if (WRITE_Y_TRACE && (iteration % _Y.get_thinning_factor() == 0) && !_fix_Y) {
			_Y_trace_file.writeln(_Y.get_full_Y_binary_vector());
		}
	}

	void _read_Y_from_file(const std::string &filename);

	template<bool IsSimulation, bool FixZ> void _update_all_Z(size_t iteration) {
		if (iteration == 0 && WRITE_Z_TRACE && _Z_trace_files.empty() && !_fix_Z) {
			for (const auto &tree : _trees) {
				std::vector<size_t> Z_trace_header;
				for (size_t i = 0; i < tree->get_Z().total_size_of_container_space(); ++i) {
					Z_trace_header.push_back(i);
				}
				_Z_trace_files.emplace_back(_prefix + "_" + tree->get_tree_name() + "_Z_trace.txt", Z_trace_header,
				                            "\t");
			}
		}

		for (auto &_tree : _trees) { _tree->update_Z_and_mus_and_branch_lengths<IsSimulation, FixZ>(_Y); }
		if (_fix_Z) { return; }
		if (iteration % _Y.get_thinning_factor() == 0 && WRITE_Z_TRACE) {
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

		std::array<size_t, 2> line{};
		coretools::TOutputFile file(filename, header, "\t");
		double fraction;
		if constexpr (WriteFullY) {
			for (size_t i = 0; i < _Y.total_size_of_container_space(); ++i) {
				auto [found, position] = _Y.binary_search(i);
				if (found) {
					line = {i, _Y.is_one(position)};
				} else {
					line = {i, 0};
				}
				std::vector<size_t> leaf_index_of_Y = _Y.get_multi_dimensional_index(i);
				std::vector<std::string> node_names;
				for (size_t idx = 0; idx < leaf_index_of_Y.size(); ++idx) {
					size_t node_idx = _trees[idx]->get_node_index_from_leaf_index(leaf_index_of_Y[idx]);
					node_names.push_back(_trees[idx]->get_node_id(node_idx));
				};
				fraction = _Y.get_fraction_of_ones(i);
				file.writeln(line, node_names, fraction);
			}
		} else {
			for (size_t i = 0; i < _Y.size(); ++i) {
				line                                = {_Y[i].get_linear_index_in_Y_space(), _Y[i].is_one()};
				fraction                            = _Y.get_fraction_of_ones(_Y[i].get_linear_index_in_Y_space());
				std::vector<size_t> leaf_index_of_Y = _Y.get_multi_dimensional_index(i);
				std::vector<std::string> node_names;
				for (size_t idx = 0; idx < leaf_index_of_Y.size(); ++idx) {
					size_t node_idx = _trees[idx]->get_node_index_from_leaf_index(leaf_index_of_Y[idx]);
					node_names.push_back(_trees[idx]->get_node_id(node_idx));
				};
				file.writeln(line, node_names, fraction);
			}
		}
	};

public:
	TMarkovField(size_t n_iterations, std::vector<std::unique_ptr<TTree>> &Trees, std::string _prefix);
	~TMarkovField() = default;

	// updates
	void update(TLotus &lotus, size_t iteration);

	// simulation
	void simulate(TLotus &lotus);

	// get Y
	[[nodiscard]] const TStorageYVector &get_Y_vector() const;
	[[nodiscard]] const TStorageY &get_Y(size_t index) const;
	[[nodiscard]] size_t size_Y() const;

	// functions to perform stuff on Y after burnin / MCMC finished
	void burninHasFinished();
	void MCMCHasFinished();

	static size_t get_num_iterations_simulation() {
		return coretools::instances::parameters().get("num_iterations", 5000);
	}
};

#endif // ACOL_TMARKOVFIELD_H
