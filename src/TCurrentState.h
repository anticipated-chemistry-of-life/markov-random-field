//
// Created by madleina on 10.01.25.
//

#ifndef ACOL_TCURRENTSTATE_H
#define ACOL_TCURRENTSTATE_H

#include "smart_binary_search.h"

//-----------------------------------
// TCurrentState
//-----------------------------------

class TCurrentState {
private:
	// current state of Y
	std::vector<bool> _current_state_Y;
	std::vector<bool> _exists_in_Y;
	std::vector<size_t> _index_in_TStorageYVector;

	// current state of Z
	std::vector<bool> _current_state_Z;
	std::vector<bool> _exists_in_Z;
	std::vector<size_t> _index_in_TStorageZVector;

	// increment and tree
	size_t _increment;
	const TTree &_tree;

public:
	TCurrentState(const TTree &tree, size_t increment);

	void fill(const std::vector<size_t> &start_index_in_leaves_space, const TStorageYVector &Y,
	          const TStorageZVector &Z);
	void fill_Y(const std::vector<size_t> &start_index_in_leaves_space, size_t num_nodes_to_parse,
	            const TStorageYVector &Y);
	void fill_Z(const std::vector<size_t> &start_index_in_leaves_space, size_t num_nodes_to_parse,
	            const TStorageZVector &Z);
	void fill_Y_along_last_dim(const std::vector<size_t> &start_index_in_leaves_space, size_t num_nodes_to_parse,
	                           const TStorageYVector &Y);
	void fill_Z_along_last_dim(const std::vector<size_t> &start_index_in_leaves_space, size_t num_nodes_to_parse,
	                           const TStorageZVector &Z);

	bool get(size_t index_in_tree) const;
	bool get(size_t index_in_tree, size_t offset_leaves, size_t offset_internals) const;
	bool get_Y(size_t ix) const;
	bool get_Z(size_t ix) const;
	void set(size_t index_in_tree, bool value);

	size_t get_index_in_TStorageVector(size_t index_in_tree) const;

	bool exists_in_TStorageVector(size_t index_in_tree) const;
};

//-----------------------------------
// TSheet
//-----------------------------------

class TSheet {
private:
	// dimension along which sheet runs (can be any except last dimension, as this one is covered by K)
	size_t _dim_ix;
	// tree of dimension corresponding to _dim_ix
	const TTree &_tree;
	// tree of last dimension
	const TTree &_tree_last_dim;

	size_t _start_ix_in_leaves_space_last_dim = 0;

	// the sheet consists of multiple TCurrentStates (length of vector = number of nodes in _tree)
	// each _cur_states[i] is of size K
	std::vector<TCurrentState> _cur_states;

public:
	TSheet(size_t dim_ix, const TTree &tree, const TTree &tree_last_dim);
	~TSheet() = default;

	void fill(std::vector<size_t> start_index_in_leaves_space, size_t K, const TStorageYVector &Y);

	bool get(size_t node_index_in_tree_of_dim, size_t leaf_index_in_tree_of_last_dim) const;
};

#endif // ACOL_TCURRENTSTATE_H
