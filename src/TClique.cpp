//
// Created by madleina on 22.10.24.
//

#include "TClique.h"
#include "TCurrentState.h"
#include "TStorageYVector.h"
#include "TStorageZ.h"
#include "TStorageZVector.h"
#include "TTree.h"
#include "coretools/Math/TSumLog.h"
#include <cstddef>
#include <unordered_set>
#include <vector>

double TMatrices::_Delta;

TClique::TClique(const std::vector<size_t> &start_index_in_leaves_space, size_t variable_dimension, size_t n_nodes,
                 size_t increment) {
	_start_index_in_leaves_space = start_index_in_leaves_space;
	_variable_dimension          = variable_dimension;
	_n_nodes                     = n_nodes;
	_increment                   = increment;
}

TCurrentState TClique::create_current_state(const TStorageYVector &Y, TStorageZVector &Z, const TTree &tree) {
	TCurrentState current_state(tree, this->_increment, tree.get_number_of_leaves(),
	                            tree.get_number_of_internal_nodes());
	current_state.fill(_start_index_in_leaves_space, Y, Z);

	return current_state;
}

std::vector<TStorageZ>
TClique::update_Z(std::vector<double> &joint_prob_density, TCurrentState &current_state, TStorageZVector &Z,
                  const TTree *tree, TypeAlpha alpha, const TypeParamBinBranches *binned_branch_lengths,
                  const std::vector<size_t> &leaves_and_internal_nodes_without_roots_indices) const {
	std::vector<TStorageZ> linear_indices_in_Z_space_to_insert;

	const double stationary_0 = get_stationary_probability(false, alpha);

	for (const auto index_in_tree : tree->get_internal_nodes()) {
		// prepare log probabilities for the two possible states
		std::array<coretools::TSumLogProbability, 2> sum_log;
		const auto &node = tree->get_node(index_in_tree);
		if (node.isRoot()) { // calculate stationary
			_calculate_log_prob_root(stationary_0, sum_log);
		} else { // calculate P(node = 0 | parent) and P(node = 1 | parent)
			// note: for compatibility with update of Y, we need to pass leaf_index_in_tree_of_last_dim, but this
			// doesn't matter for this update (just pass 0)
			// Note: need to take oldValue of binned_branch_lengths because they are updated before the loop starts
			const auto bin_branch_len =
			    binned_branch_lengths->oldValue(leaves_and_internal_nodes_without_roots_indices[index_in_tree]);
			calculate_log_prob_parent_to_node(index_in_tree, bin_branch_len, tree, 0, current_state, sum_log);
		}

		// calculate P(child | node = 0) and P(child | node = 1) for all children of node
		_calculate_log_prob_node_to_children(index_in_tree, tree, current_state, sum_log, binned_branch_lengths,
		                                     leaves_and_internal_nodes_without_roots_indices);

		// sample new state and update Z accordingly
		const double log_prob_0 = sum_log[0].getSum();
		const double log_prob_1 = sum_log[1].getSum();
		bool new_state          = sample(log_prob_0, log_prob_1);

		if (new_state) {
			joint_prob_density[omp_get_thread_num()] += log_prob_1;
		} else {
			joint_prob_density[omp_get_thread_num()] += log_prob_0;
		}

		_update_current_state(Z, current_state, index_in_tree, new_state, linear_indices_in_Z_space_to_insert, tree);
	}

	return linear_indices_in_Z_space_to_insert;
}

void TClique::_update_current_state(TStorageZVector &Z, TCurrentState &current_state, size_t index_in_tree,
                                    bool new_state, std::vector<TStorageZ> &linear_indices_in_Z_space_to_insert,
                                    const TTree *tree) const {
	auto index_in_TStorageZVector = current_state.get_index_in_TStorageVector(index_in_tree);
	// 1 -> 0: can simply set Z to zero (element exists already)
	if (current_state.get(index_in_tree) && !new_state) { Z.set_to_zero(index_in_TStorageZVector); }
	// 0 -> 1: it depends if element already exists in Z (i.e. has been a 1 previously)
	if (!current_state.get(index_in_tree) && new_state) {
		if (current_state.exists_in_TStorageVector(index_in_tree)) { // it does exist -> set it to one
			Z.set_to_one(index_in_TStorageZVector);
		} else { // it does not exist -> remember linear index in Z to be inserted later!
			auto multidim_index_in_Z_space                 = _start_index_in_leaves_space;
			multidim_index_in_Z_space[_variable_dimension] = tree->get_index_within_internal_nodes(index_in_tree);
			size_t linear_index_in_Z_space                 = Z.get_linear_index_in_Z_space(multidim_index_in_Z_space);
			linear_indices_in_Z_space_to_insert.emplace_back(linear_index_in_Z_space);
		}
	}
	current_state.set(index_in_tree, new_state);
}

std::vector<TStorageZ>
TClique::initialize_Z_from_children(TCurrentState &current_state, TStorageZVector &Z, const TTree *tree,
                                    const TypeParamBinBranches *binned_branch_lengths,
                                    const std::vector<size_t> &leaves_and_internal_nodes_without_roots_indices) const {

	// initialise vector that will insert the Z not in parallel
	std::vector<TStorageZ> linear_indices_in_Z_space_to_insert;
	std::unordered_set<size_t> remaining;
	std::vector<int> was_initilized(tree->get_number_of_nodes(), false);
	for (size_t leaf_index : tree->get_leaf_nodes()) { was_initilized[leaf_index] = true; }

	// start with the leaves
	for (size_t internal_node_index : tree->get_internal_nodes()) { remaining.insert(internal_node_index); }

	// bottom-up update of Z
	while (!remaining.empty()) {
		for (auto it = remaining.begin(); it != remaining.end();) {
			const size_t node_index = *it;
			bool ready              = true;
			for (size_t child_index : tree->get_node(node_index).children_indices_in_tree()) {
				if (!was_initilized[child_index]) {
					ready = false;
					break;
				}
			}

			if (ready) {
				_set_Z_to_MLE(node_index, current_state, Z, tree, binned_branch_lengths,
				              leaves_and_internal_nodes_without_roots_indices, linear_indices_in_Z_space_to_insert);
				was_initilized[node_index] = true;
				it                         = remaining.erase(it);
			} else {
				++it;
			}
		}
	}

	return linear_indices_in_Z_space_to_insert;
}

void TClique::_set_Z_to_MLE(size_t node_index, TCurrentState &current_state, TStorageZVector &Z, const TTree *tree,
                            const TypeParamBinBranches *binned_branch_lengths,
                            const std::vector<size_t> &leaves_and_internal_nodes_without_roots_indices,
                            std::vector<TStorageZ> &linear_indices_in_Z_space_to_insert) const {
	std::array<coretools::TSumLogProbability, 2> sum_log;

	_calculate_log_prob_node_to_children(node_index, tree, current_state, sum_log, binned_branch_lengths,
	                                     leaves_and_internal_nodes_without_roots_indices);

	// sample new state and update Z accordingly
	const double log_prob_0 = sum_log[0].getSum();
	const double log_prob_1 = sum_log[1].getSum();

	bool new_state = log_prob_1 > log_prob_0;
	_update_current_state(Z, current_state, node_index, new_state, linear_indices_in_Z_space_to_insert, tree);
}

void TClique::_calculate_log_prob_root(double stationary_0, std::array<coretools::TSumLogProbability, 2> &sum_log) {
	sum_log[0].add(stationary_0);
	sum_log[1].add(1.0 - stationary_0);
}

void TClique::_calculate_log_prob_node_to_children(
    size_t index_in_tree, const TTree *tree, const TCurrentState &current_state,
    std::array<coretools::TSumLogProbability, 2> &sum_log, const TypeParamBinBranches *binned_branch_lengths,
    const std::vector<size_t> &leaves_and_internal_nodes_without_roots_indices) const {
	const auto &node = tree->get_node(index_in_tree);
	for (const auto &child_index : node.children_indices_in_tree()) {
		// Note: need to take old value of branch length because new values were proposed before the loop started
		auto bin_length = binned_branch_lengths->oldValue(leaves_and_internal_nodes_without_roots_indices[child_index]);
		const auto &matrix_for_bin = _cur_matrices[bin_length];
		const bool child_state     = current_state.get(child_index);
		for (size_t i = 0; i < 2; ++i) { // loop over possible values (0 or 1) of the node
			sum_log[i].add(matrix_for_bin(i, child_state));
		}
	}
}

void TClique::update_counter_leaves_state_1(bool new_state, bool old_state) {
	const int diff = (int)new_state - (int)old_state;
	update_counter_leaves_state_1(diff);
}
void TClique::update_counter_leaves_state_1(int difference) {
	assert((int)_counter_leaves_state_1 + difference >= 0);
	_counter_leaves_state_1 += difference;
}

bool sample(std::array<coretools::TSumLogProbability, 2> &sum_log) {
	const double log_Q = sum_log[1].getSum() - sum_log[0].getSum();
	return coretools::TAcceptOddsRatio::accept(log_Q);
}

bool sample(double log_prob_0, double log_prob_1) {
	const double log_Q = log_prob_1 - log_prob_0;
	return coretools::TAcceptOddsRatio::accept(log_Q);
}

size_t TClique::_get_parent_index(size_t index_in_tree, const TTree *tree) {
	const auto &node                  = tree->get_node(index_in_tree);
	const size_t parent_index_in_tree = node.parentIndex_in_tree();
	return parent_index_in_tree;
}
