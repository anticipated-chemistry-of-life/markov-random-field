//
// Created by madleina on 10.01.25.
//

#include "TMarkovField.h"
#include "TClique.h"
#include "TCurrentState.h"
#include "TLotus.h"
#include "TStorageYVector.h"
#include "TTree.h"
#include "coretools/Main/TParameters.h"
#include <cstddef>
#include <vector>

TMarkovField::TMarkovField(size_t n_iterations, std::vector<std::unique_ptr<TTree>> &Trees)
    : _trees(Trees), _clique_last_dim(*_trees.back().get(), 1) {
	using namespace coretools::instances;

	// read K (sheet size for updating Y)
	_K = parameters().get("K", 100);

	// read: fix Y or Z?
	_fix_Y = !parameters().get("Y.update", true);
	if (_fix_Y) { logfile().list("Will fix Y during the MCMC."); }
	_fix_Z = !parameters().get("Z.update", true);
	if (_fix_Z) { logfile().list("Will fix Z during the MCMC."); }

	// number of outer loops = the number of times to repeat K such that all leaves of the last dimension are parsed
	_num_outer_loops = std::ceil((double)_trees.back()->get_number_of_leaves() / (double)_K);

	// set number of leaves per dimension (set the last dimension to one)
	_num_leaves_per_dim_except_last.resize(_trees.size(), 1);
	for (size_t i = 0; i < _trees.size() - 1; ++i) {
		_num_leaves_per_dim_except_last[i] = _trees[i]->get_number_of_leaves();
	}

	// create sheets: one per dimension except the last dimension
	_sheets.reserve(_trees.size() - 1);
	for (size_t i = 0; i < _trees.size() - 1; ++i) { _sheets.emplace_back(i, *_trees[i].get(), *_trees.back().get()); }

	// initialize Y
	auto num_leaves_per_dim   = _num_leaves_per_dim_except_last;
	num_leaves_per_dim.back() = _trees.back()->get_number_of_leaves();
	_Y.initialize(n_iterations, num_leaves_per_dim);
}

bool TMarkovField::_need_to_update_sheet(size_t sheet_ix, const std::vector<size_t> &start_index_in_leaves_space,
                                         const std::vector<size_t> &previous_ix) const {
	for (size_t j = 0; j < _sheets.size(); ++j) { // loop over all sheets
		if (j == sheet_ix) { continue; }          // ignore current sheet
		if (start_index_in_leaves_space[j] != previous_ix[j]) { return true; }
	}
	return false;
}

void TMarkovField::_update_sheets(bool first, const std::vector<size_t> &start_index_in_leaves_space,
                                  const std::vector<size_t> &previous_ix, size_t K_cur_sheet) {
	for (size_t j = 0; j < _sheets.size(); ++j) {
		if (first || _need_to_update_sheet(j, start_index_in_leaves_space, previous_ix)) {
			// first iteration or different index than before -> re-compute sheet
			_sheets[j].fill(start_index_in_leaves_space, K_cur_sheet, _Y);
		}
	}
}

void TMarkovField::_fill_clique_along_last_dim(const std::vector<size_t> &start_index_in_leaves_space) {
	// fill all Y along last dimension
	_clique_last_dim.fill_Y_along_last_dim(start_index_in_leaves_space, _trees.back()->get_number_of_leaves(), _Y);
	// fill all Z along last dimension
	_clique_last_dim.fill_Z_along_last_dim(start_index_in_leaves_space, _trees.back()->get_number_of_internal_nodes(),
	                                       _trees.back()->get_Z());
}

void TMarkovField::_calculate_log_prob_field(const std::vector<size_t> &index_in_leaves_space,
                                             std::array<coretools::TSumLogProbability, 2> &sum_log) const {
	for (size_t dim = 0; dim < _trees.size(); ++dim) {
		// get relevant clique
		const auto &clique         = _trees[dim]->get_clique(index_in_leaves_space);
		// translate index in leaves to the index in tree
		const size_t index_in_tree = _trees[dim]->get_node_index_from_leaf_index(index_in_leaves_space[dim]);
		// get leaf index in tree of last dimension
		const size_t leaf_index_in_tree_of_last_dim = index_in_leaves_space.back();
		// calculate P(parent | node = 0) and P(parent | node = 1)
		// Note: leaves can never be roots -> they always have a parent (no need to bother with stationary)
		if (dim == _trees.size() - 1) { // last dim -> use _clique_last_dim
			clique.calculate_log_prob_parent_to_node(
			    index_in_tree, _trees[dim]->get_binned_branch_length(index_in_tree), _trees[dim].get(),
			    leaf_index_in_tree_of_last_dim, _clique_last_dim, sum_log);
		} else { // use sheet
			clique.calculate_log_prob_parent_to_node(
			    index_in_tree, _trees[dim]->get_binned_branch_length(index_in_tree), _trees[dim].get(),
			    leaf_index_in_tree_of_last_dim, _sheets[dim], sum_log);
		}
	}
}

void TMarkovField::_update_counter_1_cliques(bool new_state, bool old_state,
                                             const std::vector<size_t> &index_in_leaves_space) {
	// update counter of leaves with value 1 for all dimensions except the last one
	// reason: we parallelize over the last dimension -> can not update the counter there, as this would result in race
	// condition
	for (size_t dim = 0; dim < _trees.size() - 1; ++dim) {
		_trees[dim]->get_clique(index_in_leaves_space).update_counter_leaves_state_1(new_state, old_state);
	}
}

void TMarkovField::_update_cur_LL_lotus(TLotus &lotus, std::vector<coretools::TSumLogProbability> &new_LL) {
	double sum_new_LL = 0.0;
	for (size_t i = 0; i < new_LL.size(); ++i) {
		// loop over all LL (stored per thread) and sum
		sum_new_LL += new_LL[i].getSum();
	}
	lotus.update_cur_LL(sum_new_LL);
}

int TMarkovField::_set_new_Y(bool new_state, const std::vector<size_t> &index_in_leaves_space,
                             std::vector<TStorageY> &linear_indices_in_Y_space_to_insert) {
	const size_t leaf_index_in_tree_of_last_dim = index_in_leaves_space.back();

	// get current state, exists and index in TStorageYVector
	// get it from _clique_last_dim, as this is re-computed for each update and is thus reliable with respect to indices
	// of where TStorageYVector is a one
	const auto [cur_state, exists_in_TStorageYVector, index_in_TStorageYVector] =
	    _clique_last_dim.get_state_exist_ix_TStorageYVector(leaf_index_in_tree_of_last_dim);

	if (cur_state && !new_state) { // 1 -> 0: set Y to zero
		_Y.set_to_zero(index_in_TStorageYVector);
	} else if (!cur_state && new_state) { // 0 -> 1
		if (exists_in_TStorageYVector) {
			_Y.set_to_one(index_in_TStorageYVector);
		} else { // does not yet exist -> remember linear index in Y to insert later
			const size_t linear_index_in_Y_space = _Y.get_linear_index_in_Y_space(index_in_leaves_space);
			linear_indices_in_Y_space_to_insert.emplace_back(linear_index_in_Y_space);
		}
	}

	// update values of 1-valued leaves in clique
	_update_counter_1_cliques(new_state, cur_state, index_in_leaves_space);
	int diff_counter_1_in_last_dim = (int)new_state - (int)cur_state;

	// update value in _sheets
	for (size_t dim = 0; dim < _trees.size() - 1; ++dim) {
		// translate leaf index to node index
		const size_t node_index_in_tree_of_dim =
		    _trees[dim]->get_node_index_from_leaf_index(index_in_leaves_space[dim]);
		_sheets[dim].set(node_index_in_tree_of_dim, leaf_index_in_tree_of_last_dim, new_state);
	}
	// Note: no need to update in _clique_last_dim, as this will anyway be overwritten for next Y

	return diff_counter_1_in_last_dim;
}

void TMarkovField::update(TLotus &lotus) {
	_update_all_Y<false>(lotus);
	_update_all_Z<false>();
}

void TMarkovField::_calc_lotus_LL(const std::vector<size_t> &index_in_leaves_space, size_t leaf_index_last_dim,
                                  std::array<double, 2> &prob, TLotus &lotus) {
	const bool cur_state = _clique_last_dim.get_Y(leaf_index_last_dim);
	lotus.calculate_LL_update_Y(index_in_leaves_space, cur_state, prob);
}

void TMarkovField::_prepare_lotus_LL(const std::vector<size_t> &start_index_in_leaves_space, size_t K_cur_sheet,
                                     TLotus &lotus) {
	lotus.fill_tmp_state_along_last_dim(start_index_in_leaves_space, K_cur_sheet);
}

void TMarkovField::simulate(TLotus &lotus) {
	// For simulation we always draw from the prior. (top-down)
	// 1. Draw branch len -> draw mus
	// 2. For every tree draw the root from those mus and the we BFS sample all the internal nodes based on the state of
	// the parent. At the end also draw the state of the leaves. For each tree we know which leaves and which internal
	// node are 0 or 1.
	// 2. draw Z since its unique per tree
	// 3. Draw Y (this function)
	// 4. Draw Lotus.
	// 5. From that simulated Lotus can we go back and infer the prior params
	//
	for (size_t tree_index = 0; tree_index < _trees.size(); ++tree_index) {
		auto &tree = _trees[tree_index];
		tree->simulate_Z(tree_index);
	}

	_simulate_Y();

	// for iteration in 1->max_iteration, (max_iteration should be passed from CLI)
	// we use tree.update_Z(); and then
	// update Y where likelihood of data is always one so it doesn't matter.
	size_t max_iteration = get_num_iterations_simulation();

	std::vector<size_t> Y_trace_header;
	for (size_t i = 0; i < _Y.total_size_of_container_space(); ++i) { Y_trace_header.push_back(i); }
	coretools::TOutputFile Y_trace_file("acol_simulated_Y_trace.txt", Y_trace_header, "\t");
	std::vector<coretools::TOutputFile> Z_trace_files;
	for (const auto &tree : _trees) {
		std::vector<size_t> Z_trace_header;
		for (size_t i = 0; i < tree->get_Z().total_size_of_container_space(); ++i) { Z_trace_header.push_back(i); }
		Z_trace_files.emplace_back("acol_simulated_Z_" + tree->get_tree_name() + "_trace.txt", Z_trace_header, "\t");
	}
	for (size_t iteration = 0; iteration < max_iteration; ++iteration) {
		_update_all_Y<true>(lotus);
		_update_all_Z<true>();
		_Y.add_to_counter(iteration);
		if (iteration % 10 == 0) {
			Y_trace_file.writeln(_Y.get_full_Y_binary_vector());
			for (size_t tree_idx = 0; tree_idx < _trees.size(); ++tree_idx) {
				const auto &tree = _trees[tree_idx];
				Z_trace_files[tree_idx].writeln(tree->get_Z().get_full_Z_binary_vector());
			}
		}
	}
	_write_Y_to_file<true>("acol_simulated_Y.txt");
	for (const auto &tree : _trees) {
		tree->write_Z_to_file<true>("acol_simulated_Z_" + tree->get_tree_name() + ".txt");
	}
}

void TMarkovField::_simulate_Y() {
	// to sample Y we need to know the state of the parent for each leaf that is represented in a Y entry.
	// We are going to iterate over all possible Y and sample givent the product of probabilities of the child given the
	// parent.
	// set number of leaves per dimension (set the last dimension to one)
	for (size_t linear_index_in_leaves_space = 0; linear_index_in_leaves_space < _Y.total_size_of_container_space();
	     ++linear_index_in_leaves_space) {
		std::vector<size_t> multidim_index_in_Y = _Y.get_multi_dimensional_index(linear_index_in_leaves_space);
		std::array<coretools::TSumLogProbability, 2> sum_log;
		for (size_t dim = 0; dim < _trees.size(); ++dim) {
			// get relevant clique
			const auto &clique = _trees[dim]->get_clique(multidim_index_in_Y);
			TCurrentState current_state(*_trees[dim].get(), clique.get_increment(), _trees[dim]->get_number_of_leaves(),
			                            _trees[dim]->get_number_of_internal_nodes());
			// translate index in leaves to the index in tree
			const size_t index_in_tree = _trees[dim]->get_node_index_from_leaf_index(multidim_index_in_Y[dim]);
			// calculate P(parent | node = 0) and P(parent | node = 1)
			// Note: leaves can never be roots -> they always have a parent (no need to bother with stationary)
			clique.calculate_log_prob_parent_to_node(index_in_tree,
			                                         _trees[dim]->get_binned_branch_length(index_in_tree),
			                                         _trees[dim].get(), 0, current_state, sum_log);
		}
		bool y_state = sample(sum_log);
		if (y_state) {
			_Y.insert_one(linear_index_in_leaves_space);
			for (auto &_tree : _trees) {
				// get relevant clique
				auto &clique = _tree->get_clique(multidim_index_in_Y);
				clique.update_counter_leaves_state_1(true, false);
			}
		}

		// TODO : refactor also Y sampling
	}
}

void TMarkovField::burninHasFinished() { _Y.reset_counts(); }

void TMarkovField::MCMCHasFinished() {
	// TODO: write function to write the posterior state of Y to file
}

const TStorageYVector &TMarkovField::get_Y_vector() const { return _Y; }
const TStorageY &TMarkovField::get_Y(size_t index_in_TStorageYVector) const { return _Y[index_in_TStorageYVector]; }
size_t TMarkovField::size_Y() const { return _Y.size(); }
