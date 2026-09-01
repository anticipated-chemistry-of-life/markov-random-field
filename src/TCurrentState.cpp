//
// Created by madleina on 10.01.25.
//

#include "TCurrentState.h"
#include "constants.h"
#include "tree/TPhylogeny.h"

//-----------------------------------
// TCurrentState
//-----------------------------------

TCurrentState::TCurrentState(const TPhylogeny &topology, size_t increment)
    : _increment(increment), _topology(topology) {}
TCurrentState::TCurrentState(const TPhylogeny &topology, size_t increment, size_t size_of_Y,
                             size_t size_of_Z)
    : _increment(increment), _topology(topology) {
	_current_state_Y.resize(size_of_Y, false);
	_exists_in_Y.resize(size_of_Y, false);
	_index_in_Y.resize(size_of_Y);

	_current_state_Z.resize(size_of_Z, false);
	_exists_in_Z.resize(size_of_Z, false);
	_index_in_Z.resize(size_of_Z);
}

void TCurrentState::fill(const IndexArray &start_index_in_leaves_space, const TFieldStorage &Y,
                         const TNodeStateStorage &Z) {
	fill_Y(start_index_in_leaves_space, _topology.n_leaves(),
	       Y); // parse all Y (all leaves)
	fill_Z(start_index_in_leaves_space, _topology.n_nodes(),
	       Z); // parse all Z (every node, leaves included)
}

void TCurrentState::fill_Y_along_last_dim(const IndexArray &start_index_in_leaves_space,
                                          size_t num_nodes_to_parse, const TFieldStorage &Y) {
	// along the last dimension -> increment is 1; the storage picks the line to walk from that
	// and its own shape (see TStorageYMatrix::fill_current_state).
	// _index_in_Y now holds the linear index in Y space of each parsed cell.
	Y.fill_current_state(start_index_in_leaves_space, num_nodes_to_parse, /*increment=*/1,
	                     _current_state_Y, _exists_in_Y, _index_in_Y);
}

void TCurrentState::fill_Z_along_last_dim(const IndexArray &start_index_in_leaves_space,
                                          size_t num_nodes_to_parse, const TNodeStateStorage &Z) {
	// along the last dimension -> increment is 1; the storage picks the line to walk from that
	// and its own shape (see TStorageZMatrix::fill_current_state).
	// _index_in_Z now holds the linear index in Z space of each parsed cell.
	Z.fill_current_state(start_index_in_leaves_space, num_nodes_to_parse, /*increment=*/1,
	                     _current_state_Z, _exists_in_Z, _index_in_Z);
}

void TCurrentState::fill_Y(const IndexArray &start_index_in_leaves_space, size_t num_nodes_to_parse,
                           const TFieldStorage &Y) {
	// The storage picks the line to walk from the increment and its own shape (see
	// TStorageYMatrix::fill_current_state).
	// _index_in_Y now holds the linear index in Y space of each parsed cell.
	Y.fill_current_state(start_index_in_leaves_space, num_nodes_to_parse, _increment,
	                     _current_state_Y, _exists_in_Y, _index_in_Y);
}

void TCurrentState::fill_Z(const IndexArray &start_index_in_leaves_space, size_t num_nodes_to_parse,
                           const TNodeStateStorage &Z) {
	// The storage picks the line to walk from the increment and its own shape (see
	// TStorageZMatrix::fill_current_state).
	// _index_in_Z now holds the linear index in Z space of each parsed cell.
	Z.fill_current_state(start_index_in_leaves_space, num_nodes_to_parse, _increment,
	                     _current_state_Z, _exists_in_Z, _index_in_Z);
}

bool TCurrentState::get(size_t index_in_tree) const {
	// A leaf's state is read from the field and an internal node's from the node state. The node
	// state now spans the leaves too, but nothing reads those rows yet (ADR-0005), so this still
	// routes a leaf to the field. Its row there is its node index either way (ADR-0004).
	if (_topology.is_leaf(index_in_tree)) { return get_Y(_topology.leaf_index(index_in_tree)); }
	return get_Z(index_in_tree);
}

bool TCurrentState::get_Y(size_t ix) const { return _current_state_Y[ix]; }

bool TCurrentState::get_Z(size_t ix) const { return _current_state_Z[ix]; }

void TCurrentState::set(size_t index_in_tree, bool value) {
	if (_topology.is_leaf(index_in_tree)) {
		_current_state_Y[_topology.leaf_index(index_in_tree)] = value;
	} else {
		_current_state_Z[index_in_tree] = value;
	}
}

void TCurrentState::set_Y(size_t index_in_leaves, bool value) {
	_current_state_Y[index_in_leaves] = value;
}

size_t TCurrentState::get_linear_index_in_storage(size_t index_in_tree) const {
	if (_topology.is_leaf(index_in_tree)) {
		return _index_in_Y[_topology.leaf_index(index_in_tree)];
	}
	return _index_in_Z[index_in_tree];
}

bool TCurrentState::exists_in_storage(size_t index_in_tree) const {
	if (_topology.is_leaf(index_in_tree)) {
		return _exists_in_Y[_topology.leaf_index(index_in_tree)];
	}
	return _exists_in_Z[index_in_tree];
}

TCellInY TCurrentState::get_state_exist_ix_in_Y(size_t index_in_leaves) const {
	return {static_cast<bool>(_current_state_Y[index_in_leaves]),
	        static_cast<bool>(_exists_in_Y[index_in_leaves]), _index_in_Y[index_in_leaves]};
}

//-----------------------------------
// TSheet
//-----------------------------------
TSheet::TSheet(size_t dim_ix, const TPhylogeny &topology, const TPhylogeny &topology_last_dim)
    : _dim_ix(dim_ix), _topology(topology), _topology_last_dim(topology_last_dim) {

	// create vector of current states that make up the sheet
	constexpr static size_t increment = 1; // always 1, since we move along the last dimension
	_cur_states.reserve(_topology.n_nodes());
	for (size_t i = 0; i < _topology.n_nodes(); ++i) {
		_cur_states.emplace_back(_topology_last_dim, increment);
	}
}

void TSheet::fill(const IndexArray &start_index_in_leaves_space, size_t K, const TFieldStorage &Y,
                  const TNodeStateStorage &Z) {
	// Worksharing fill: this runs on the team created in TMarkovField::_update_all_Y (all threads
	// call it), so we use `omp for`/`omp single` rather than spawning our own team. If ever called
	// outside a parallel region the constructs are orphaned and execute sequentially, which is also
	// correct.

	// start index and how many Y need to be parsed (scalar write -> single thread)
#pragma omp single
	{ _start_ix_in_leaves_space_last_dim = start_index_in_leaves_space.back(); }

#pragma omp for schedule(static)
	for (size_t i = 0; i < _topology.n_nodes();
	     ++i) { // loop over all nodes along current dimension
		IndexArray local_start_index_in_leaves_space =
		    start_index_in_leaves_space; // thread-local copy

		if (_topology.is_leaf(
		        i)) { // only fill Y (ignore Z along last dimensions, not needed when updating)
			local_start_index_in_leaves_space[_dim_ix] = _topology.leaf_index(i);
			_cur_states[i].fill_Y_along_last_dim(local_start_index_in_leaves_space, K, Y);
		} else {
			// leaf space in every dimension except _dim_ix, whose node state spans all nodes and
			// is therefore indexed by the node index itself
			local_start_index_in_leaves_space[_dim_ix] = i;
			// fill Z. There are as many Z as there are leaves along the last dimension
			// use z of your own dimension for filling
			// Note: we do not need to fill Y here, as there are no Y when the node(i) is internal
			// Note: this will not fill nodes that are not part of Z, i.e. which are internal in the
			// last dimension
			_cur_states[i].fill_Z_along_last_dim(local_start_index_in_leaves_space, K, Z);
		}
	}
}

bool TSheet::get(size_t node_index_in_tree_of_dim, size_t leaf_index_in_tree_of_last_dim) const {
	// calculate index in Y or Z: leaf index in last dimension, relative to start index
	const size_t ix = leaf_index_in_tree_of_last_dim - _start_ix_in_leaves_space_last_dim;

	if (_topology.is_leaf(node_index_in_tree_of_dim)) {
		// leaf in all dimensions -> return Y
		return _cur_states[node_index_in_tree_of_dim].get_Y(ix);
	}

	// is internal in _dim_ix but a leaf in last dim -> stored in Z
	return _cur_states[node_index_in_tree_of_dim].get_Z(ix);
}

void TSheet::set(size_t node_index_in_tree_of_dim, size_t leaf_index_in_tree_of_last_dim,
                 bool value) {
	// calculate index in Y: leaf index in last dimension, relative to start index
	const size_t ix = leaf_index_in_tree_of_last_dim - _start_ix_in_leaves_space_last_dim;
	_cur_states[node_index_in_tree_of_dim].set_Y(ix, value);
}

TCellInY TSheet::get_state_exist_ix_in_Y(size_t node_index_in_tree_of_dim,
                                         size_t leaf_index_in_tree_of_last_dim) const {
	// calculate index in Y: leaf index in last dimension, relative to start index
	const size_t ix = leaf_index_in_tree_of_last_dim - _start_ix_in_leaves_space_last_dim;
	return _cur_states[node_index_in_tree_of_dim].get_state_exist_ix_in_Y(ix);
}
