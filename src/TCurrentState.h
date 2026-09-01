//
// Created by madleina on 10.01.25.
//

#ifndef ACOL_TCURRENTSTATE_H
#define ACOL_TCURRENTSTATE_H

#include "constants.h"
#include "storages/storage_backend.h"
#include <cstddef>
#include <vector>

class TPhylogeny; // only the topology is needed: no parameters, no cliques, no field

/// One cell of the field as a clique holds it. The field update reads all three at once, because it
/// names the cell to its stream of uniforms before it knows what state the cell is about to take.
struct TCellInY {
	bool state          = false; ///< the state the storage holds for this cell
	bool is_stored      = false; ///< whether the storage holds the cell at all
	size_t linear_index = 0;     ///< where the cell sits, in field space
};

//-----------------------------------
// TCurrentState
//-----------------------------------

class TCurrentState {
private:
	// current state of Y
	std::vector<uint8_t> _current_state_Y;
	std::vector<uint8_t> _exists_in_Y;
	// Linear index in Y space of each parsed cell (the matrix maps it to (row, col)).
	std::vector<size_t> _index_in_Y;

	// current state of Z
	std::vector<uint8_t> _current_state_Z;
	std::vector<uint8_t> _exists_in_Z;
	// Linear index in Z space of each parsed cell (Z is now a TSparseMatrix, like Y).
	std::vector<size_t> _index_in_Z;

	// increment and topology
	size_t _increment;
	const TPhylogeny &_topology;

public:
	/// One increment serves both containers, and that is not an accident worth rediscovering: a
	/// clique's increment is the product of the leaf counts of the dimensions *after* its own
	/// (tree/clique/tree_clique.cpp), and a clique's own dimension is never one of those. So the
	/// only extent that ADR-0005 changed -- the node state's own dimension -- is never a factor in
	/// it, and the same stride walks the field and the node state alike.
	TCurrentState(const TPhylogeny &topology, size_t increment);
	TCurrentState(const TPhylogeny &topology, size_t increment, size_t size_of_Y, size_t size_of_Z);

	void fill(const IndexArray &start_index_in_leaves_space, const TFieldStorage &Y,
	          const TNodeStateStorage &Z);
	void fill_Y(const IndexArray &start_index_in_leaves_space, size_t num_nodes_to_parse,
	            const TFieldStorage &Y);
	void fill_Z(const IndexArray &start_index_in_leaves_space, size_t num_nodes_to_parse,
	            const TNodeStateStorage &Z);
	void fill_Y_along_last_dim(const IndexArray &start_index_in_leaves_space,
	                           size_t num_nodes_to_parse, const TFieldStorage &Y);
	void fill_Z_along_last_dim(const IndexArray &start_index_in_leaves_space,
	                           size_t num_nodes_to_parse, const TNodeStateStorage &Z);

	bool get(size_t index_in_tree) const;
	bool get_Y(size_t ix) const;
	bool get_Z(size_t ix) const;
	void set(size_t index_in_tree, bool value);
	void set_Y(size_t index_in_leaves, bool value);
	size_t size_of_Y() const { return _current_state_Y.size(); }
	size_t size_of_Z() const { return _current_state_Z.size(); }

	/// The linear index of a node's cell, and whether that cell is stored: in the field for a
	/// leaf, in the node state for anything else. Which of the two a node lands in is what
	/// these hide, so a caller updating one node needs no case distinction.
	size_t get_linear_index_in_storage(size_t index_in_tree) const;
	bool exists_in_storage(size_t index_in_tree) const;

	[[nodiscard]] TCellInY get_state_exist_ix_in_Y(size_t index_in_leaves) const;
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

	/// `Z` is the node state of the dimension this sheet runs along, passed in the same way as
	/// the field rather than reached for through a tree.
	void fill(const IndexArray &start_index_in_leaves_space, size_t K, const TFieldStorage &Y,
	          const TNodeStateStorage &Z);

	bool get(size_t node_index_in_tree_of_dim, size_t leaf_index_in_tree_of_last_dim) const;
	void set(size_t node_index_in_tree_of_dim, size_t leaf_index_in_tree_of_last_dim, bool value);

	[[nodiscard]] TCellInY get_state_exist_ix_in_Y(size_t node_index_in_tree_of_dim,
	                                               size_t leaf_index_in_tree_of_last_dim) const;
};

#endif // ACOL_TCURRENTSTATE_H
