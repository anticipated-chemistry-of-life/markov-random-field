//
// Created by madleina on 22.10.24.
//

#include "TClique.h"
#include "TStorageYVector.h"
#include "TStorageZVector.h"
#include "TTree.h"
#include "coretools/Math/TSumLog.h"
#include "coretools/devtools.h"
#include "smart_binary_search.h"
#include <cstddef>
#include <vector>

std::vector<size_t> TClique::update_Z(const TStorageYVector &Y, TStorageZVector &Z, const TTree &tree) const {
	TCurrentState current_state(this->_start_index, this->_increment, Y, Z, tree);
	std::vector<size_t> linear_indices_in_Z_space_to_insert;

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
			TClique::_update_current_state(Z, current_state, index_in_tree, new_state,
			                               linear_indices_in_Z_space_to_insert);
		} else if (node.isInternalNode()) {
			auto sum_log   = this->_compute_internal_nodes(node, current_state);
			auto sum_log_0 = std::get<0>(sum_log);
			auto sum_log_1 = std::get<1>(sum_log);
			bool new_state = this->_compute_new_state(current_state, tree, node, sum_log_0, sum_log_1);
			TClique::_update_current_state(Z, current_state, index_in_tree, new_state,
			                               linear_indices_in_Z_space_to_insert);
		}
	}

	return linear_indices_in_Z_space_to_insert;
}
