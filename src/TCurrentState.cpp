//
// Created by madleina on 10.01.25.
//

#include "TCurrentState.h"
#include "TTree.h"

//-----------------------------------
// TCurrentState
//-----------------------------------

TCurrentState::TCurrentState(const TTree &tree, size_t increment) : _increment(increment), _tree(tree) {}

void TCurrentState::fill(const std::vector<size_t> &start_index_in_leaves_space, const TStorageYVector &Y,
                         const TStorageZVector &Z) {
	_fill_Y(start_index_in_leaves_space, _tree.get_number_of_leaves(), Y);         // parse all Y (all leaves)
	_fill_Z(start_index_in_leaves_space, _tree.get_number_of_internal_nodes(), Z); // parse all Z (all internal nodes)
}

void TCurrentState::fill(const std::vector<size_t> &start_index_in_leaves_space, size_t num_nodes_to_parse,
                         const TStorageYVector &Y, const TStorageZVector &Z) {
	_fill_Y(start_index_in_leaves_space, num_nodes_to_parse, Y); // only parse num_nodes_to_parse in Y
	_fill_Z(start_index_in_leaves_space, num_nodes_to_parse, Z); // only parse num_nodes_to_parse in Z
}

void TCurrentState::_fill_Y(const std::vector<size_t> &start_index_in_leaves_space, size_t num_nodes_to_parse,
                            const TStorageYVector &Y) {
	std::tie(_current_state_Y, _exists_in_Y, _index_in_TStorageYVector) =
	    fill_current_state(Y, num_nodes_to_parse, start_index_in_leaves_space, _increment, Y.size());
}

void TCurrentState::_fill_Z(const std::vector<size_t> &start_index_in_leaves_space, size_t num_nodes_to_parse,
                            const TStorageZVector &Z) {
	std::tie(_current_state_Z, _exists_in_Z, _index_in_TStorageZVector) =
	    fill_current_state(Z, num_nodes_to_parse, start_index_in_leaves_space, _increment, Z.size());
}

bool TCurrentState::get(size_t index_in_tree) const {
	if (_tree.get_node(index_in_tree).isLeaf()) {
		return _current_state_Y[_tree.get_index_within_leaves(index_in_tree)];
	}
	return _current_state_Z[_tree.get_index_within_internal_nodes(index_in_tree)];
}

void TCurrentState::set(size_t index_in_tree, bool value) {
	if (_tree.get_node(index_in_tree).isLeaf()) {
		_current_state_Y[_tree.get_index_within_leaves(index_in_tree)] = value;
	} else {
		_current_state_Z[_tree.get_index_within_internal_nodes(index_in_tree)] = value;
	}
}

size_t TCurrentState::get_index_in_TStorageVector(size_t index_in_tree) const {
	if (_tree.get_node(index_in_tree).isLeaf()) {
		return _index_in_TStorageYVector[_tree.get_index_within_leaves(index_in_tree)];
	}
	return _index_in_TStorageZVector[_tree.get_index_within_internal_nodes(index_in_tree)];
}

bool TCurrentState::exists_in_TStorageVector(size_t index_in_tree) const {
	if (_tree.get_node(index_in_tree).isLeaf()) { return _exists_in_Y[_tree.get_index_within_leaves(index_in_tree)]; }
	return _exists_in_Z[_tree.get_index_within_internal_nodes(index_in_tree)];
}

//-----------------------------------
// TSheet
//-----------------------------------
TSheet::TSheet(const TTree &tree) : _tree(tree) {}

void TSheet::initialize(size_t dim_ix, size_t num_nodes_in_dim, const TTree &tree_last_dim) {
	_dim_ix = dim_ix;

	constexpr static size_t increment = 1; // always 1, since we move along the last dimension
	_cur_states.reserve(_tree.get_number_of_nodes());
	for (size_t i = 0; i < _tree.get_number_of_nodes(); ++i) { _cur_states.emplace_back(tree_last_dim, increment); }
}

void TSheet::fill(const std::vector<size_t> &start_index_in_leaves_space, size_t num_nodes_to_use_in_last_dim,
                  const TStorageYVector &Y, const std::vector<TStorageZVector> &Z_of_all_dimensions) {
	for (size_t i = 0; i < _tree.get_number_of_nodes(); ++i) { // loop over all nodes along current dimension
		if (_tree.isLeaf(i)) {                                 // use z of last dimension for filling
			_cur_states[i].fill(start_index_in_leaves_space, num_nodes_to_use_in_last_dim, Y,
			                    Z_of_all_dimensions.back());
		} else { // use z of your own dimension for filling
			_cur_states[i].fill(start_index_in_leaves_space, num_nodes_to_use_in_last_dim, Y,
			                    Z_of_all_dimensions[_dim_ix]);
		}
	}
}