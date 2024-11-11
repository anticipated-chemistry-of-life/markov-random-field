//
// Created by madleina on 22.10.24.
//

#include "TClique.h"
#include "TStorageYVector.h"
#include "TStorageZVector.h"
#include "TTree.h"
#include "coretools/Math/TAcceptOddsRation.h"
#include "coretools/Math/TSumLog.h"
#include "coretools/devtools.h"
#include "update_current_state.h"
#include <cstddef>
#include <filesystem>
#include <vector>

void TClique::update_Z(const TStorageYVector &Y, const TStorageZVector &Z, const TTree &tree) {
	TCurrentState current_state(this->_start_index, this->_increment, Y, Z, tree);

	double stationary_0 = get_stationary_probability(false);
	double stationary_1 = get_stationary_probability(true);

	// if you are a root, the probability is just the stationary probability (according to our model)
	for (const auto root_index : tree.get_root_nodes()) {
		const auto &node = tree.get_node(root_index);
		coretools::TSumLogProbability sum_log_0;
		coretools::TSumLogProbability sum_log_1;
		sum_log_0.add(stationary_0);
		sum_log_1.add(stationary_1);
		for (const auto &child_index : node.children()) {
			auto bin_length            = tree.get_node(child_index).get_branch_length_bin();
			const auto &matrix_for_bin = this->_matrices[bin_length];
			double prob_0_to_child     = matrix_for_bin(0, current_state.get(child_index));
			double prob_1_to_child     = matrix_for_bin(1, current_state.get(child_index));
			sum_log_0.add(prob_0_to_child);
			sum_log_1.add(prob_1_to_child);
		}

		const double log_Q = sum_log_1.getSum() - sum_log_0.getSum();
		bool new_state     = coretools::TAcceptOddsRatio::accept(log_Q);
		if (new_state != current_state.get(root_index)) { current_state.set(root_index, new_state); }
	}

	for (const auto internal_index : tree.get_internal_nodes_without_roots()) {
		const auto &node = tree.get_node(internal_index);
		coretools::TSumLogProbability sum_log_0;
		coretools::TSumLogProbability sum_log_1;
		const auto &parent         = tree.get_node(node.parentIndex());
		auto bin_length            = parent.get_branch_length_bin();
		const auto &matrix_for_bin = this->_matrices[bin_length];
		double prob_0_to_parent    = matrix_for_bin(current_state.get(node.parentIndex()), 0);
		double prob_1_to_parent    = matrix_for_bin(current_state.get(node.parentIndex()), 1);
		sum_log_0.add(prob_0_to_parent);
		sum_log_1.add(prob_1_to_parent);

		for (const auto &child_index : node.children()) {
			auto bin_length            = tree.get_node(child_index).get_branch_length_bin();
			const auto &matrix_for_bin = _matrices[bin_length];
			double prob_0_to_child     = matrix_for_bin(0, current_state.get(child_index));
			double prob_1_to_child     = matrix_for_bin(1, current_state.get(child_index));
			sum_log_0.add(prob_0_to_child);
			sum_log_1.add(prob_1_to_child);
		}
		const double log_Q = sum_log_1.getSum() - sum_log_0.getSum();
		bool new_state     = coretools::TAcceptOddsRatio::accept(log_Q);
		if (new_state != current_state.get(internal_index)) { current_state.set(internal_index, new_state); }
	}
}
