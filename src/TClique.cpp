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
#include <vector>

TClique::TClique(const std::vector<size_t> &start_index_in_leaves_space, size_t variable_dimension, size_t n_nodes,
                 size_t increment) {
	_start_index_in_leaves_space = start_index_in_leaves_space;
	_variable_dimension          = variable_dimension;
	_n_nodes                     = n_nodes;
	_increment                   = increment;
}

bool TClique::_compute_new_state(const TCurrentState &current_state, const TTree &tree, const TNode &node,
                                 coretools::TSumLogProbability &sum_log_0,
                                 coretools::TSumLogProbability &sum_log_1) const {
	for (const auto &child_index : node.children_indices_in_tree()) {
		auto bin_length            = tree.get_node(child_index).get_branch_length_bin();
		const auto &matrix_for_bin = _matrices[bin_length];
		double prob_0_to_child     = matrix_for_bin(0, current_state.get(child_index));
		double prob_1_to_child     = matrix_for_bin(1, current_state.get(child_index));
		sum_log_0.add(prob_0_to_child);
		sum_log_1.add(prob_1_to_child);
	}
	const double log_Q = sum_log_1.getSum() - sum_log_0.getSum();
	bool new_state     = coretools::TAcceptOddsRatio::accept(log_Q);
	return new_state;
};

std::vector<TStorageZ> TClique::update_Z(const TStorageYVector &Y, TStorageZVector &Z, const TTree &tree) const {
	TCurrentState current_state(tree, this->_increment);
	current_state.fill(_start_index_in_leaves_space, Y, Z);
	std::vector<TStorageZ> linear_indices_in_Z_space_to_insert;

	const double stationary_0 = this->get_stationary_probability(false);
	const double stationary_1 = this->get_stationary_probability(true);

	// if you are a root, the probability is just the stationary probability (according to our model)
	for (const auto index_in_tree : tree.get_internal_nodes()) {
		const auto &node = tree.get_node(index_in_tree);
		if (node.isRoot()) {
			auto sum_log   = TClique::_compute_roots(stationary_0, stationary_1);
			auto sum_log_0 = std::get<0>(sum_log);
			auto sum_log_1 = std::get<1>(sum_log);
			bool new_state = this->_compute_new_state(current_state, tree, node, sum_log_0, sum_log_1);
			_update_current_state(Z, current_state, index_in_tree, new_state, linear_indices_in_Z_space_to_insert,
			                      tree);
		} else if (node.isInternalNode()) {
			auto sum_log   = this->_compute_internal_nodes(node, current_state);
			auto sum_log_0 = std::get<0>(sum_log);
			auto sum_log_1 = std::get<1>(sum_log);
			bool new_state = this->_compute_new_state(current_state, tree, node, sum_log_0, sum_log_1);
			_update_current_state(Z, current_state, index_in_tree, new_state, linear_indices_in_Z_space_to_insert,
			                      tree);
		}
	}

	return linear_indices_in_Z_space_to_insert;
}

void TClique::_update_current_state(TStorageZVector &Z, const TCurrentState &current_state, size_t index_in_tree,
                                    bool new_state, std::vector<TStorageZ> &linear_indices_in_Z_space_to_insert,
                                    const TTree &tree) const {
	auto index_in_TStorageZVector = current_state.get_index_in_TStorageVector(index_in_tree);
	if (current_state.get(index_in_tree) && !new_state) { Z.set_to_zero(index_in_TStorageZVector); }
	if (!current_state.get(index_in_tree) && new_state) {
		if (current_state.exists_in_TStorageVector(index_in_tree)) {
			Z.set_to_zero(index_in_TStorageZVector);
		} else {
			auto multidim_index_in_Z_space                 = _start_index_in_leaves_space;
			multidim_index_in_Z_space[_variable_dimension] = tree.get_index_within_internal_nodes(index_in_tree);
			size_t linear_index_in_Z_space                 = Z.get_linear_index_in_Z_space(multidim_index_in_Z_space);
			linear_indices_in_Z_space_to_insert.emplace_back(linear_index_in_Z_space);
		}
	}
}

std::tuple<coretools::TSumLogProbability, coretools::TSumLogProbability>
TClique::_compute_internal_nodes(const TNode &node, const TCurrentState &current_state) const {
	coretools::TSumLogProbability sum_log_0;
	coretools::TSumLogProbability sum_log_1;
	auto bin_length            = node.get_branch_length_bin();
	const auto &matrix_for_bin = this->_matrices[bin_length];
	double prob_0_to_parent    = matrix_for_bin(current_state.get(node.parentIndex_in_tree()), 0);
	double prob_1_to_parent    = matrix_for_bin(current_state.get(node.parentIndex_in_tree()), 1);
	sum_log_0.add(prob_0_to_parent);
	sum_log_1.add(prob_1_to_parent);

	return {sum_log_0, sum_log_1};
}
