//
// Created by madleina on 22.10.24.
//

#include "TClique.h"
#include "TStorageYVector.h"
#include "TStorageZVector.h"
#include "TTree.h"
#include <cstddef>
#include <vector>

void TClique::update_Z(const TStorageYVector &Y, const TStorageZVector &Z, const TTree &tree) {
	std::vector<bool> current_state(_n_nodes, false);
	for (size_t i = 0; i < _n_nodes; ++i) {
		if (tree.get_node(i).isLeaf()) {
			std::vector<size_t> multi_dim_index_in_Y  = _start_index;
			multi_dim_index_in_Y[_variable_dimension] = tree.get_index_within_leaves(i);
			current_state[i]                          = Y.is_one(multi_dim_index_in_Y);
		} else {
			std::vector<size_t> multi_dim_index_in_Z  = _start_index;
			multi_dim_index_in_Z[_variable_dimension] = tree.get_index_within_internal_nodes(i);
			current_state[i]                          = Z.is_one(multi_dim_index_in_Z);
		}
	}
	// TODO
}
