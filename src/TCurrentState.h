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

	// functions
	void _fill_Y(const std::vector<size_t> &start_index_in_leaves_space, size_t num_nodes_to_parse,
	             const TStorageYVector &Y);
	void _fill_Z(const std::vector<size_t> &start_index_in_leaves_space, size_t num_nodes_to_parse,
	             const TStorageZVector &Z);

public:
	TCurrentState(const TTree &tree, size_t increment);

	void fill(const std::vector<size_t> &start_index_in_leaves_space, const TStorageYVector &Y,
	          const TStorageZVector &Z);
	void fill(const std::vector<size_t> &start_index_in_leaves_space, size_t num_nodes_to_parse,
	          const TStorageYVector &Y, const TStorageZVector &Z);

	bool get(size_t index_in_tree) const;
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
	const TTree &_tree;

	std::vector<TCurrentState> _cur_states;

public:
	TSheet(const TTree &_tree);
	~TSheet() = default;

	void initialize(size_t dim_ix, size_t num_nodes_in_dim, const TTree &tree_last_dim);
	void fill(const std::vector<size_t> &start_index_in_leaves_space, size_t index_end_in_last_dim,
	          const TStorageYVector &Y, const std::vector<TStorageZVector> &Z_of_all_dimensions);
};

#endif // ACOL_TCURRENTSTATE_H
