//
// Created by madleina on 22.10.24.
//

#include "TClique.h"
#include "TStorageYVector.h"
#include "TStorageZ.h"
#include "TStorageZVector.h"
#include "TTree.h"
#include "coretools/Math/TSumLog.h"
#include <cstddef>
#include <utility>
#include <vector>

TClique::TClique(const std::vector<size_t> &start_index_in_leaves_space, size_t variable_dimension, size_t n_nodes,
                 size_t increment) {
	_start_index_in_leaves_space = start_index_in_leaves_space;
	_variable_dimension          = variable_dimension;
	_n_nodes                     = n_nodes;
	_increment                   = increment;
}

std::vector<TStorageZ> TClique::update_Z(const TStorageYVector &Y, TStorageZVector &Z, const TTree &tree) const {
	TCurrentState current_state(tree, this->_increment);
	current_state.fill(_start_index_in_leaves_space, Y, Z);
	std::vector<TStorageZ> linear_indices_in_Z_space_to_insert;

	const double stationary_0 = get_stationary_probability(false);

	// prepare log probabilities for the two possible states
	std::array<coretools::TSumLogProbability, 2> sum_log;

	for (const auto index_in_tree : tree.get_internal_nodes()) {
		const auto &node = tree.get_node(index_in_tree);
		if (node.isRoot()) { // calculate stationary
			_calculate_log_prob_root(stationary_0, sum_log);
		} else { // calculate P(node = 0 | parent) and P(node = 1 | parent)
			// note: for compatibility with update of Y, we need to pass leaf_index_in_tree_of_last_dim, but this
			// doesn't matter for this update (just pass 0)
			calculate_log_prob_parent_to_node(index_in_tree, tree, 0, current_state, sum_log);
		}

		// calculate P(child | node = 0) and P(child | node = 1) for all children of node
		_calculate_log_prob_node_to_children(index_in_tree, tree, current_state, sum_log);

		// sample new state and update Z accordingly
		bool new_state = sample(sum_log);
		_update_current_state(Z, current_state, index_in_tree, new_state, linear_indices_in_Z_space_to_insert, tree);
	}

	return linear_indices_in_Z_space_to_insert;
}

void TClique::_update_current_state(TStorageZVector &Z, const TCurrentState &current_state, size_t index_in_tree,
                                    bool new_state, std::vector<TStorageZ> &linear_indices_in_Z_space_to_insert,
                                    const TTree &tree) const {
	auto index_in_TStorageZVector = current_state.get_index_in_TStorageVector(index_in_tree);
	// 1 -> 0: can simply set Z to zero (element exists already)
	if (current_state.get(index_in_tree) && !new_state) { Z.set_to_zero(index_in_TStorageZVector); }
	// 0 -> 1: it depends if element already exists in Z (i.e. has been a 1 previously)
	if (!current_state.get(index_in_tree) && new_state) {
		if (current_state.exists_in_TStorageVector(index_in_tree)) { // it does exist -> set it to one
			Z.set_to_one(index_in_TStorageZVector);
		} else { // it does not exist -> remember linear index in Z to be inserted later!
			auto multidim_index_in_Z_space                 = _start_index_in_leaves_space;
			multidim_index_in_Z_space[_variable_dimension] = tree.get_index_within_internal_nodes(index_in_tree);
			size_t linear_index_in_Z_space                 = Z.get_linear_index_in_Z_space(multidim_index_in_Z_space);
			linear_indices_in_Z_space_to_insert.emplace_back(linear_index_in_Z_space);
		}
	}
}

void TClique::_calculate_log_prob_root(double stationary_0, std::array<coretools::TSumLogProbability, 2> &sum_log) {
	sum_log[0].add(stationary_0);
	sum_log[1].add(1.0 - stationary_0);
}

void TClique::_calculate_log_prob_node_to_children(size_t index_in_tree, const TTree &tree,
                                                   const TCurrentState &current_state,
                                                   std::array<coretools::TSumLogProbability, 2> &sum_log) const {
	const auto &node = tree.get_node(index_in_tree);
	for (const auto &child_index : node.children_indices_in_tree()) {
		auto bin_length            = tree.get_node(child_index).get_branch_length_bin();
		const auto &matrix_for_bin = _matrices[bin_length];
		for (size_t i = 0; i < 2; ++i) { // loop over possible values (0 or 1) of the node
			sum_log[i].add(matrix_for_bin(i, current_state.get(child_index)));
		}
	}
}

bool sample(std::array<coretools::TSumLogProbability, 2> &sum_log) {
	const double log_Q = sum_log[1].getSum() - sum_log[0].getSum();
	return coretools::TAcceptOddsRatio::accept(log_Q);
}

std::pair<size_t, TypeBinBranches> TClique::_get_parent_index_and_bin_length(size_t index_in_tree, const TTree &tree) {
	// calculates log P(node = 0 | parent) and log P(node = 1 | parent)
	const auto &node                  = tree.get_node(index_in_tree);
	auto bin_length                   = node.get_branch_length_bin();
	const size_t parent_index_in_tree = node.parentIndex_in_tree();

	return {parent_index_in_tree, bin_length};
}
