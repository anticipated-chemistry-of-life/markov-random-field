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
	if (_tree.isLeaf(index_in_tree)) {
		return _current_state_Y[_tree.get_index_within_leaves(index_in_tree) - offset_leaves];
	}
	return _current_state_Z[_tree.get_index_within_internal_nodes(index_in_tree) - offset_internals];
}

bool TCurrentState::get_Z(size_t index_in_tree, size_t offset_internals) const {
	return _current_state_Z[_tree.get_index_within_internal_nodes(index_in_tree) - offset_internals];
}

void TCurrentState::set(size_t index_in_tree, bool value) {
	if (_tree.isLeaf(index_in_tree)) {
		_current_state_Y[_tree.get_index_within_leaves(index_in_tree)] = value;
	} else {
		_current_state_Z[_tree.get_index_within_internal_nodes(index_in_tree)] = value;
	}
}

size_t TCurrentState::get_index_in_TStorageVector(size_t index_in_tree) const {
	if (_tree.isLeaf(index_in_tree)) { return _index_in_TStorageYVector[_tree.get_index_within_leaves(index_in_tree)]; }
	return _index_in_TStorageZVector[_tree.get_index_within_internal_nodes(index_in_tree)];
}

bool TCurrentState::exists_in_TStorageVector(size_t index_in_tree) const {
	if (_tree.isLeaf(index_in_tree)) { return _exists_in_Y[_tree.get_index_within_leaves(index_in_tree)]; }
	return _exists_in_Z[_tree.get_index_within_internal_nodes(index_in_tree)];
}

//-----------------------------------
// TSheet
//-----------------------------------
TSheet::TSheet(size_t dim_ix, const TTree &tree, const TTree &tree_last_dim)
    : _tree(tree), _tree_last_dim(tree_last_dim), _dim_ix(dim_ix) {

	// create vector of current states that make up the sheet
	constexpr static size_t increment = 1; // always 1, since we move along the last dimension
	_cur_states.reserve(_tree.get_number_of_nodes());
	for (size_t i = 0; i < _tree.get_number_of_nodes(); ++i) { _cur_states.emplace_back(_tree_last_dim, increment); }
}

void TSheet::_calculate_num_Y_Z(size_t start_ix_last_dim, size_t num_nodes_to_use_in_last_dim) {
	// TODO: we could calculate this once and forever, but would need to store this outside of class

	// calculate how many Y and Z there are along the last dimension before the start of the current sheet
	_offset_leaves_in_last_dim    = 0;
	_offset_internals_in_last_dim = 0;
	for (size_t i = 0; i < start_ix_last_dim; ++i) { // loop over 0 until start node of last dimension
		if (_tree_last_dim.isLeaf(i)) {
			++_offset_leaves_in_last_dim;
		} else {
			++_offset_internals_in_last_dim;
		}
	}

	// for the current sheet: calculate the number of Y and Z (along one slice of the sheet)
	_num_Y = 0;
	_num_Z = 0;

	// calculate end index (make sure not to overshoot)
	const size_t end = std::min(start_ix_last_dim + num_nodes_to_use_in_last_dim, _tree_last_dim.get_number_of_nodes());
	for (size_t i = start_ix_last_dim; i < end; ++i) { // loop over start until end node of last dimension
		if (_tree_last_dim.isLeaf(i)) {
			++_num_Y;
		} else {
			++_num_Z;
		}
	}
}

void TSheet::fill(std::vector<size_t> start_index_in_leaves_space, size_t num_nodes_to_use_in_last_dim,
                  const TStorageYVector &Y, const std::vector<TStorageZVector> &Z_of_all_dimensions) {
	// calculate how many Y and Z need to be parsed
	_calculate_num_Y_Z(start_index_in_leaves_space.back(), num_nodes_to_use_in_last_dim);

	for (size_t i = 0; i < _tree.get_number_of_nodes(); ++i) { // loop over all nodes along current dimension
		if (_tree.isLeaf(i)) {
			start_index_in_leaves_space[_dim_ix] = _tree.get_index_within_leaves(i);
			// use z of last dimension for filling
			_cur_states[i].fill_Y_along_last_dim(start_index_in_leaves_space, _num_Y, Y);
			_cur_states[i].fill_Z_along_last_dim(start_index_in_leaves_space, _num_Z, Z_of_all_dimensions.back());
		} else {
			// get start index: leaf space in all dimensions except _dim_ix, for which we use the internal node index
			start_index_in_leaves_space[_dim_ix] = _tree.get_index_within_internal_nodes(i);
			// fill Z. There are as many Z as there are leaves along the last dimension (= num_Y)
			// use z of your own dimension for filling
			// Note: we do not need to fill Y here, as there are no Y when the node(i) is internal
			// Note: this will not fill nodes that are not part of Z, i.e. which are internal in the last dimension
			_cur_states[i].fill_Z_along_last_dim(start_index_in_leaves_space, _num_Y, Z_of_all_dimensions[_dim_ix]);
		}
	}
}

bool TSheet::get(size_t index_in_tree_of_dim, size_t index_in_tree_of_last_dim) const {
	if (_tree.isLeaf(index_in_tree_of_dim)) {
		return _cur_states[index_in_tree_of_dim].get(index_in_tree_of_last_dim, _offset_leaves_in_last_dim,
		                                             _offset_internals_in_last_dim);
	}
	if (_tree_last_dim.isLeaf(index_in_tree_of_dim)) { // is internal in _dim_ix but a leaf in last dim -> stored in Z
		return _cur_states[index_in_tree_of_dim].get_Z(index_in_tree_of_last_dim, _offset_internals_in_last_dim);
	}
	// is internal in last dim -> not part of Z or Y -> doesn't matter, but should never be accessed
	DEVERROR("Tried to access element that is not part of Z nor Y: index_in_tree_of_dim = ", index_in_tree_of_dim,
	         ", index_in_tree_of_last_dim = ", index_in_tree_of_last_dim, ". Is this a bug?");
}