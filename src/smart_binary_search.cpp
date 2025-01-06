#include "smart_binary_search.h"
#include "TTree.h"

TCurrentState::TCurrentState(const std::vector<size_t> &multidim_index_of_first_in_container_space, size_t increment,
                             const TStorageYVector &Y, const TStorageZVector &Z, const TTree &tree)
    : _multidim_index_of_first_in_container_space(multidim_index_of_first_in_container_space), _increment(increment),
      _tree(tree) {
	if (increment == 1) {
		std::tie(_current_state_Y, _exists_in_Y, _index_in_TStorageYVector) =
		    fill_current_state_easy(Y, _multidim_index_of_first_in_container_space, tree.get_number_of_leaves());
		std::tie(_current_state_Z, _exists_in_Z, _index_in_TStorageZVector) = fill_current_state_easy(
		    Z, _multidim_index_of_first_in_container_space, tree.get_number_of_internal_nodes());
	} else {
		std::tie(_current_state_Y, _exists_in_Y, _index_in_TStorageYVector) = fill_current_state_hard(
		    Y, tree.get_number_of_leaves(), _multidim_index_of_first_in_container_space, _increment, Y.size());
		std::tie(_current_state_Z, _exists_in_Z, _index_in_TStorageZVector) = fill_current_state_hard(
		    Z, tree.get_number_of_internal_nodes(), _multidim_index_of_first_in_container_space, _increment, Z.size());
	}
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
