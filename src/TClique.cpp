//
// Created by madleina on 22.10.24.
//

#include "TClique.h"
#include "TStorageYVector.h"
#include "TStorageZVector.h"
#include "TTree.h"
#include "coretools/Math/TSumLog.h"
#include "coretools/devtools.h"
#include "update_current_state.h"
#include <cstddef>
#include <iostream>
#include <vector>

std::vector<size_t> TClique::update_Z(const TStorageYVector &Y, TStorageZVector &Z, const TTree &tree) const {
	TCurrentState current_state(this->_start_index, this->_increment, Y, Z, tree);
	std::vector<size_t> Z_indices_not_to_update_in_parallel;

	double stationary_0 = this->get_stationary_probability(false);
	double stationary_1 = this->get_stationary_probability(true);

	// if you are a root, the probability is just the stationary probability (according to our model)
	for (const auto root_index : tree.get_root_nodes()) {
		OUT(root_index);
		OUT(tree.get_index_within_internal_nodes(root_index));
		const auto &node = tree.get_node(root_index);
		coretools::TSumLogProbability sum_log_0;
		coretools::TSumLogProbability sum_log_1;
		sum_log_0.add(stationary_0);
		sum_log_1.add(stationary_1);
		bool new_state  = this->_compute_new_state(current_state, tree, node, sum_log_0, sum_log_1);
		auto coordinate = current_state.get_coordinate_in_container(root_index);
		OUT(coordinate);
		if (1 == current_state.get(root_index) && 0 == new_state) {
			Z.set_to_zero(coordinate);
			current_state.set(root_index, new_state);
		}
		if (0 == current_state.get(root_index) && 1 == new_state && current_state.is_in_container(root_index)) {
			Z.set_to_zero(coordinate);
			current_state.set(root_index, new_state);
		} else if (0 == current_state.get(root_index) && 1 == new_state) {
			Z_indices_not_to_update_in_parallel.push_back(coordinate);
			current_state.set(root_index, new_state);
		}
	}

	for (const auto internal_index : tree.get_internal_nodes_without_roots()) {
		const auto &node = tree.get_node(internal_index);
		coretools::TSumLogProbability sum_log_0;
		coretools::TSumLogProbability sum_log_1;
		const auto &parent = tree.get_node(node.parentIndex());
		auto bin_length    = parent.get_branch_length_bin();
		OUT(bin_length);
		const auto &matrix_for_bin = this->_matrices[bin_length];
		matrix_for_bin.print();
		double prob_0_to_parent = matrix_for_bin(current_state.get(node.parentIndex()), 0);
		double prob_1_to_parent = matrix_for_bin(current_state.get(node.parentIndex()), 1);
		sum_log_0.add(prob_0_to_parent);
		sum_log_1.add(prob_1_to_parent);

		bool new_state  = this->_compute_new_state(current_state, tree, node, sum_log_0, sum_log_1);
		auto coordinate = current_state.get_coordinate_in_container(internal_index);
		if (1 == current_state.get(internal_index) && 0 == new_state) {
			Z.set_to_zero(coordinate);
			current_state.set(internal_index, new_state);
		}
		if (0 == current_state.get(internal_index) && 1 == new_state && current_state.is_in_container(internal_index)) {
			Z.set_to_zero(coordinate);
			current_state.set(internal_index, new_state);
		} else if (0 == current_state.get(internal_index) && 1 == new_state) {
			Z_indices_not_to_update_in_parallel.push_back(coordinate);
			current_state.set(internal_index, new_state);
		}
	}
	return Z_indices_not_to_update_in_parallel;
}
