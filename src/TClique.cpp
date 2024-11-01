//
// Created by madleina on 22.10.24.
//

#include "TClique.h"
#include "TStorageYVector.h"
#include "TStorageZVector.h"
#include "TTree.h"
#include "coretools/Math/TAcceptOddsRation.h"
#include "coretools/Math/TSumLog.h"
#include <cstddef>
#include <vector>

void TClique::update_Z(const TStorageYVector &Y, const TStorageZVector &Z, const TTree &tree) {
	std::vector<bool> current_state(_n_nodes, false);
	for (size_t i = 0; i < _n_nodes; ++i) {
		if (tree.get_node(i).isLeaf()) {
			std::vector<size_t> multi_dim_index_in_Y  = _start_index;
			multi_dim_index_in_Y[_variable_dimension] = tree.get_index_within_leaves(i);
			const auto [found, index]                 = Y.binary_search(multi_dim_index_in_Y);
			current_state[i]                          = Y.is_one(multi_dim_index_in_Y);
		} else {
			std::vector<size_t> multi_dim_index_in_Z  = _start_index;
			multi_dim_index_in_Z[_variable_dimension] = tree.get_index_within_internal_nodes(i);
			current_state[i]                          = Z.is_one(multi_dim_index_in_Z);
		}
	}

	double stationary_0 = get_stationary_probability(false);
	double stationary_1 = get_stationary_probability(true);

	// if you are a root, the probability is just the stationary probability (according to our model)
	for (const auto &root_index : tree.get_root_nodes()) {

		const auto &node = tree.get_node(root_index);
		coretools::TSumLogProbability sum_log_0;
		coretools::TSumLogProbability sum_log_1;
		sum_log_0.add(stationary_0);
		sum_log_1.add(stationary_1);
		for (const auto &child_index : node.children()) {
			auto bin_length            = tree.get_node(child_index).get_branch_length_bin();
			const auto &matrix_for_bin = _matrices[bin_length];
			double prob_0_to_child     = matrix_for_bin(0, current_state[child_index]);
			double prob_1_to_child     = matrix_for_bin(1, current_state[child_index]);
			sum_log_0.add(prob_0_to_child);
			sum_log_1.add(prob_1_to_child);
		}

		const double log_Q = sum_log_1.getSum() - sum_log_0.getSum();
		bool new_state     = coretools::TAcceptOddsRatio::accept(log_Q);
		if (new_state != current_state[root_index]) { current_state[root_index] = new_state; }

		// TODO
	}
}

std::vector<bool> TClique::update_Z_test(const TStorageYVector &Y, const TStorageZVector &Z, const TTree &tree) {
	std::vector<bool> current_state(_n_nodes, false);
	for (size_t i = 0; i < _n_nodes; ++i) {
		if (tree.get_node(i).isLeaf()) {
			std::vector<size_t> multi_dim_index_in_Y  = _start_index;
			multi_dim_index_in_Y[_variable_dimension] = tree.get_index_within_leaves(i);
			const auto [found, index]                 = Y.binary_search(multi_dim_index_in_Y);
			current_state[i]                          = Y.is_one(multi_dim_index_in_Y);
		} else {
			std::vector<size_t> multi_dim_index_in_Z  = _start_index;
			multi_dim_index_in_Z[_variable_dimension] = tree.get_index_within_internal_nodes(i);
			current_state[i]                          = Z.is_one(multi_dim_index_in_Z);
		}
	}
	return current_state;
}
