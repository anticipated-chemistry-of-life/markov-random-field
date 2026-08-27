//
// Created by madleina on 10.01.25.
//

#ifndef ACOL_TCURRENTSTATE_H
#define ACOL_TCURRENTSTATE_H

#include "constants.h"
#include "storages/y_storage/TStorageYMatrix.h"
#include "storages/z_storage/TStorageZMatrix.h"
#include <cstddef>
#include <tuple>
#include <vector>

class TPhylogeny; // only the topology is needed: no parameters, no cliques, no field

//-----------------------------------
// TCurrentState
//-----------------------------------

class TCurrentState {
private:
	// current state of Y
	std::vector<uint8_t> _current_state_Y;
	std::vector<uint8_t> _exists_in_Y;
	// Linear index in Y space of each parsed cell (the matrix maps it to (row, col)).
	std::vector<size_t> _index_in_TStorageYMatrix;

	// current state of Z
	std::vector<uint8_t> _current_state_Z;
	std::vector<uint8_t> _exists_in_Z;
	// Linear index in Z space of each parsed cell (Z is now a TSparseMatrix, like Y).
	std::vector<size_t> _index_in_TStorageZMatrix;

	// increment and topology
	size_t _increment;
	const TPhylogeny &_topology;

public:
	TCurrentState(const TPhylogeny &topology, size_t increment);
	TCurrentState(const TPhylogeny &topology, size_t increment, size_t size_of_Y, size_t size_of_Z);

	void fill(const IndexArray &start_index_in_leaves_space, const TStorageYMatrix &Y,
	          const TStorageZMatrix &Z);
	void fill_Y(const IndexArray &start_index_in_leaves_space, size_t num_nodes_to_parse,
	            const TStorageYMatrix &Y);
	void fill_Z(const IndexArray &start_index_in_leaves_space, size_t num_nodes_to_parse,
	            const TStorageZMatrix &Z);
	void fill_Y_along_last_dim(const IndexArray &start_index_in_leaves_space,
	                           size_t num_nodes_to_parse, const TStorageYMatrix &Y);
	void fill_Z_along_last_dim(const IndexArray &start_index_in_leaves_space,
	                           size_t num_nodes_to_parse, const TStorageZMatrix &Z);

	bool get(size_t index_in_tree) const;
	bool get(size_t index_in_tree, size_t offset_leaves, size_t offset_internals) const;
	bool get_Y(size_t ix) const;
	bool get_Z(size_t ix) const;
	void set(size_t index_in_tree, bool value);
	void set_Y(size_t index_in_leaves, bool value);
	size_t size_of_Y() const { return _current_state_Y.size(); }
	size_t size_of_Z() const { return _current_state_Z.size(); }

	size_t get_index_in_TStorageVector(size_t index_in_tree) const;
	bool exists_in_TStorageVector(size_t index_in_tree) const;

	std::tuple<bool, bool, size_t> get_state_exist_ix_TStorageYMatrix(size_t index_in_leaves) const;
};

//-----------------------------------
// TSheet
//-----------------------------------

class TSheet {
private:
	// dimension along which sheet runs (can be any except last dimension, as this one is covered by
	// K)
	size_t _dim_ix;
	// topology of the dimension corresponding to _dim_ix
	const TPhylogeny &_topology;
	// topology of the last dimension
	const TPhylogeny &_topology_last_dim;

	size_t _start_ix_in_leaves_space_last_dim = 0;

	// the sheet consists of multiple TCurrentStates (length of vector = number of nodes in _tree)
	// each _cur_states[i] is of size K
	std::vector<TCurrentState> _cur_states;

public:
	TSheet(size_t dim_ix, const TPhylogeny &topology, const TPhylogeny &topology_last_dim);
	~TSheet() = default;

	/// `Z` is the internal state of the dimension this sheet runs along, passed in the same way as
	/// the field rather than reached for through a tree.
	void fill(const IndexArray &start_index_in_leaves_space, size_t K, const TStorageYMatrix &Y,
	          const TStorageZMatrix &Z);

	bool get(size_t node_index_in_tree_of_dim, size_t leaf_index_in_tree_of_last_dim) const;
	void set(size_t node_index_in_tree_of_dim, size_t leaf_index_in_tree_of_last_dim, bool value);

	std::tuple<bool, bool, size_t>
	get_state_exist_ix_TStorageYMatrix(size_t node_index_in_tree_of_dim,
	                                   size_t leaf_index_in_tree_of_last_dim) const;
};

#endif // ACOL_TCURRENTSTATE_H
