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
	fill_Y(start_index_in_leaves_space, _tree.get_number_of_leaves(), Y);         // parse all Y (all leaves)
	fill_Z(start_index_in_leaves_space, _tree.get_number_of_internal_nodes(), Z); // parse all Z (all internal nodes)
}

void TCurrentState::fill_Y_along_last_dim(const std::vector<size_t> &start_index_in_leaves_space,
                                          size_t num_nodes_to_parse, const TStorageYVector &Y) {
	std::tie(_current_state_Y, _exists_in_Y, _index_in_TStorageYVector) =
	    fill_current_state<true>(Y, num_nodes_to_parse, start_index_in_leaves_space, _increment, Y.size());
}

void TCurrentState::fill_Z_along_last_dim(const std::vector<size_t> &start_index_in_leaves_space,
                                          size_t num_nodes_to_parse, const TStorageZVector &Z) {
	std::tie(_current_state_Z, _exists_in_Z, _index_in_TStorageZVector) =
	    fill_current_state<true>(Z, num_nodes_to_parse, start_index_in_leaves_space, _increment, Z.size());
}

void TCurrentState::fill_Y(const std::vector<size_t> &start_index_in_leaves_space, size_t num_nodes_to_parse,
                           const TStorageYVector &Y) {
	std::tie(_current_state_Y, _exists_in_Y, _index_in_TStorageYVector) =
	    fill_current_state(Y, num_nodes_to_parse, start_index_in_leaves_space, _increment, Y.size());
}

void TCurrentState::fill_Z(const std::vector<size_t> &start_index_in_leaves_space, size_t num_nodes_to_parse,
                           const TStorageZVector &Z) {
	std::tie(_current_state_Z, _exists_in_Z, _index_in_TStorageZVector) =
	    fill_current_state(Z, num_nodes_to_parse, start_index_in_leaves_space, _increment, Z.size());
}

bool TCurrentState::get(size_t index_in_tree) const { return get(index_in_tree, 0, 0); }

bool TCurrentState::get(size_t index_in_tree, size_t offset_leaves, size_t offset_internals) const {
	if (_tree.isLeaf(index_in_tree)) { return get_Y(_tree.get_index_within_leaves(index_in_tree) - offset_leaves); }
	return get_Z(_tree.get_index_within_internal_nodes(index_in_tree) - offset_internals);
}

bool TCurrentState::get_Y(size_t ix) const { return _current_state_Y[ix]; }

bool TCurrentState::get_Z(size_t ix) const { return _current_state_Z[ix]; }

void TCurrentState::set(size_t index_in_tree, bool value) {
	if (_tree.isLeaf(index_in_tree)) {
		_current_state_Y[_tree.get_index_within_leaves(index_in_tree)] = value;
	} else {
		_current_state_Z[_tree.get_index_within_internal_nodes(index_in_tree)] = value;
	}
}

void TCurrentState::set_Y(size_t index_in_leaves, bool value) { _current_state_Y[index_in_leaves] = value; }

size_t TCurrentState::get_index_in_TStorageVector(size_t index_in_tree) const {
	if (_tree.isLeaf(index_in_tree)) { return _index_in_TStorageYVector[_tree.get_index_within_leaves(index_in_tree)]; }
	return _index_in_TStorageZVector[_tree.get_index_within_internal_nodes(index_in_tree)];
}

bool TCurrentState::exists_in_TStorageVector(size_t index_in_tree) const {
	if (_tree.isLeaf(index_in_tree)) { return _exists_in_Y[_tree.get_index_within_leaves(index_in_tree)]; }
	return _exists_in_Z[_tree.get_index_within_internal_nodes(index_in_tree)];
}

auto TCurrentState::get_state_exist_ix_TStorageYVector(size_t index_in_leaves) const {
	const bool state  = _current_state_Y[index_in_leaves];
	const bool exists = _exists_in_Y[index_in_leaves];
	const size_t ix   = _index_in_TStorageYVector[index_in_leaves];
	return std::tuple<bool, size_t, size_t>(state, exists, ix);
}

//-----------------------------------
// TSheet
//-----------------------------------
TSheet::TSheet(size_t dim_ix, const TTree &tree, const TTree &tree_last_dim)
    : _dim_ix(dim_ix), _tree(tree), _tree_last_dim(tree_last_dim) {

	// create vector of current states that make up the sheet
	constexpr static size_t increment = 1; // always 1, since we move along the last dimension
	_cur_states.reserve(_tree.get_number_of_nodes());
	for (size_t i = 0; i < _tree.get_number_of_nodes(); ++i) { _cur_states.emplace_back(_tree_last_dim, increment); }
}

void TSheet::fill(std::vector<size_t> start_index_in_leaves_space, size_t K, const TStorageYVector &Y) {
	// start index and how many Y need to be parsed
	_start_ix_in_leaves_space_last_dim = start_index_in_leaves_space.back();

#pragma omp parallel for num_threads(NUMBER_OF_THREADS)
	for (size_t i = 0; i < _tree.get_number_of_nodes(); ++i) { // loop over all nodes along current dimension
		if (_tree.isLeaf(i)) { // only fill Y (ignore Z along last dimensions, not needed when updating)
			start_index_in_leaves_space[_dim_ix] = _tree.get_index_within_leaves(i);
			_cur_states[i].fill_Y_along_last_dim(start_index_in_leaves_space, K, Y);
		} else {
			// get start index: leaf space in all dimensions except _dim_ix, for which we use the internal node index
			start_index_in_leaves_space[_dim_ix] = _tree.get_index_within_internal_nodes(i);
			// fill Z. There are as many Z as there are leaves along the last dimension
			// use z of your own dimension for filling
			// Note: we do not need to fill Y here, as there are no Y when the node(i) is internal
			// Note: this will not fill nodes that are not part of Z, i.e. which are internal in the last dimension
			_cur_states[i].fill_Z_along_last_dim(start_index_in_leaves_space, K, _tree.get_Z());
		}
	}
}

bool TSheet::get(size_t node_index_in_tree_of_dim, size_t leaf_index_in_tree_of_last_dim) const {
	// calculate index in Y or Z: leaf index in last dimension, relative to start index
	const size_t ix = leaf_index_in_tree_of_last_dim - _start_ix_in_leaves_space_last_dim;

	if (_tree.isLeaf(node_index_in_tree_of_dim)) {
		// leaf in all dimensions -> return Y
		return _cur_states[node_index_in_tree_of_dim].get_Y(ix);
	}

	// is internal in _dim_ix but a leaf in last dim -> stored in Z
	return _cur_states[node_index_in_tree_of_dim].get_Z(ix);
}

void TSheet::set(size_t node_index_in_tree_of_dim, size_t leaf_index_in_tree_of_last_dim, bool value) {
	// calculate index in Y: leaf index in last dimension, relative to start index
	const size_t ix = leaf_index_in_tree_of_last_dim - _start_ix_in_leaves_space_last_dim;
	_cur_states[node_index_in_tree_of_dim].set_Y(ix, value);
}

auto TSheet::get_state_exist_ix_TStorageYVector(size_t node_index_in_tree_of_dim,
                                                size_t leaf_index_in_tree_of_last_dim) {
	// calculate index in Y: leaf index in last dimension, relative to start index
	const size_t ix = leaf_index_in_tree_of_last_dim - _start_ix_in_leaves_space_last_dim;
	return _cur_states[node_index_in_tree_of_dim].get_state_exist_ix_TStorageYVector(ix);
}