//
// Created by madleina on 22.10.24.
//

#ifndef ACOL_TCLIQUE_H
#define ACOL_TCLIQUE_H

#include "Types.h"
#include "constants.h"
#include "coretools/Math/TSumLog.h"
#include "storages/storage_backend.h"
#include "tree/branch/TTransitionGrid.h"
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

class TPhylogeny;
class TTree;

/// The states a node-state walk has assigned so far, for the one clique it is walking.
///
/// The walk goes in post-order, so a parent is reached after its children and has to see the states
/// they were just given. The dense storage would show them, because every write lands in place. The
/// sparse one defers a write to a cell it does not yet hold until after the parallel region, and
/// would still answer with the old value -- so the two backends would part company inside a single
/// update. Carrying the states here is what keeps them on the same chain.
///
/// Seeded from the storages and indexed by node index. This is the whole of what the current-state
/// class was still needed for, once the field's own update stopped needing it: a vector scoped to
/// the loop that uses it, rather than a class with a fill protocol on both storages.
class TCliqueWalkStates {
private:
	std::vector<uint8_t> _state;

public:
	explicit TCliqueWalkStates(size_t n_nodes) : _state(n_nodes, 0) {}
	[[nodiscard]] bool get(size_t node) const { return _state[node] != 0; }
	void set(size_t node, bool state) { _state[node] = static_cast<uint8_t>(state); }
};

/** Class representing a clique in our model. A clique is defined as having a set of nodes that are
 * all leaves in all dimensions except one. Each clique has a transition grid, and the start index
 * of the nodes in the tree. The start index, the variable dimension, and the number of nodes are
 * needed to get the correct indices in our multidimensional space Y and Z.
 */
class TClique {
private:
	/// This clique's two-state process, discretised onto the tree's bin grid. Set once the
	/// parameters exist (TTree::guessInitialValues) and replaced wholesale whenever a proposal on
	/// alpha or nu is accepted; there is no mutable "try" copy.
	std::optional<TTransitionGrid> _transition_grid;

	// info about size and dimensionality of clique
	IndexArray _start_index_in_leaves_space;
	size_t _variable_dimension;
	size_t _n_nodes;
	size_t _increment;

	/// @brief Calculates the log probability of the root given the stationary distribution
	static void _calculate_log_prob_root(double stationary_0,
	                                     std::array<coretools::TSumLogProbability, 2> &sum_log);

	/// Give `node` its new state: in place where the storage already holds the cell, deferred to
	/// after the parallel region where it does not, because inserting reallocates a sparse row.
	void _write_new_state(TInternalStateStorage &Z, TCliqueWalkStates &walk, size_t node,
	                      bool new_state,
	                      std::vector<size_t> &linear_indices_in_Z_space_to_insert) const;

	/// @brief Calculates the log probability of a node to its children
	void _calculate_log_prob_node_to_children(
	    size_t index_in_tree, const TTree *tree, const TCliqueWalkStates &walk,
	    std::array<coretools::TSumLogProbability, 2> &sum_log) const;

	/// @brief Sets Z given the maximal likelihood given its children. This was created to avoid
	/// that Z is stuck in a state and cannot change.
	/// @param node_index The index of the internal node we want to set
	/// @param current_state The current state of the clique.
	/// @param Z the Z vector of that tree (i.e that clique)
	/// @param tree the tree of interest
	/// @param linear_indices_in_Z_space_to_insert Same as the variable name
	void _set_Z_to_MLE(size_t node_index, TCliqueWalkStates &walk, TInternalStateStorage &Z,
	                   const TTree *tree,
	                   std::vector<size_t> &linear_indices_in_Z_space_to_insert) const;

	static size_t _get_parent_index(size_t index_in_tree, const TTree *tree);

public:
	TClique(const IndexArray &start_index, size_t variable_dimension, size_t n_nodes,
	        size_t increment);
	~TClique() = default;

	/// @brief Install this clique's process. Called once the parameters exist, and again whenever a
	/// proposal on alpha or nu is accepted.
	void set_transition_grid(TTransitionGrid grid) { _transition_grid = std::move(grid); }

	/// @brief This clique's current process. Throws if the parameters have not been drawn yet,
	/// which used to read as a grid of zeros instead.
	[[nodiscard]] const TTransitionGrid &transition_grid() const {
		return _transition_grid.value();
	}

	/// @brief Update the Z dimension for this clique.
	/// @param Y The current state of the Y dimension.
	/// @param Z The current state of the Z dimension.
	/// @param tree The tree.
	/// Every state this clique starts from, read straight out of the storages. The caller keeps it
	/// for the whole of the clique's turn, because the parameter and branch-length moves that
	/// follow the node-state walk have to see the states that walk assigned.
	[[nodiscard]] TCliqueWalkStates read_states(const TFieldStorage &Y,
	                                            const TInternalStateStorage &Z,
	                                            const TPhylogeny &topology) const;

	std::vector<size_t> update_Z(std::vector<double> &joint_prob_density, TCliqueWalkStates &walk,
	                             TInternalStateStorage &Z, const TTree *tree) const;

	std::vector<size_t> initialize_Z_from_children(TCliqueWalkStates &walk,
	                                               TInternalStateStorage &Z,
	                                               const TTree *tree) const;

	/// The cell node `node` of this clique occupies, in the column the clique runs at. Setting the
	/// last dimension first and the variable one second is what makes this right for a clique along
	/// either dimension: when they are the same dimension, the node wins.
	[[nodiscard]] IndexArray cell_of(size_t node, size_t leaf_index_in_tree_of_last_dim) const {
		IndexArray cell           = _start_index_in_leaves_space;
		cell.back()               = leaf_index_in_tree_of_last_dim;
		cell[_variable_dimension] = node;
		return cell;
	}

	/// The state of `node` as the storages hold it. A leaf's state is the field's -- a leaf's index
	/// in leaf space is its node index (ADR-0004) -- and any other node's is the tree's own node
	/// state. This is the question the current-state class used to answer from a cache.
	[[nodiscard]] bool state_of(const TFieldStorage &Y, const TInternalStateStorage &Z,
	                            const TPhylogeny &topology, size_t node,
	                            size_t leaf_index_in_tree_of_last_dim) const;

	/// @brief Return the number of nodes in the clique
	/// @return Return the number of nodes in the clique
	[[nodiscard]] size_t get_number_of_nodes() const { return _n_nodes; }

	/// @brief Calculates the log probability of a node to its parent, under this clique's current
	/// process.
	/// The parent's state is the caller's to supply, because where it comes from differs: the
	/// field's update reads it from the node state, and the node state's own walk reads it from the
	/// states that walk has already assigned.
	void
	calculate_log_prob_parent_to_node(TypeBinnedBranchLengths binned_branch_length,
	                                  bool state_of_parent,
	                                  std::array<coretools::TSumLogProbability, 2> &sum_log) const {
		const auto &process = transition_grid();
		for (size_t i = 0; i < 2; ++i) { // loop over possible values (0 or 1) of the node
			sum_log[i].add(process.probability(binned_branch_length, state_of_parent, i));
		}
	}

	/// @return Returns the jump size of the clique
	[[nodiscard]] size_t get_increment() const { return _increment; }

	/// @return Returns the start index in the leaf space of that specific clique
	[[nodiscard]] const IndexArray &get_start_index_in_leaf_space() const {
		return _start_index_in_leaves_space;
	}

	/// @brief P(node | parent) under an explicitly given process, so a Metropolis proposal can ask
	/// the same question of the current grid and of its candidate.
	double calculate_prob_to_parent(size_t index_in_tree, const TTree *tree,
	                                TypeBinnedBranchLengths binned_branch_length,
	                                const TCliqueWalkStates &current_state,
	                                const TTransitionGrid &process) const {
		size_t parent_index = _get_parent_index(index_in_tree, tree);

		bool parent_state = current_state.get(parent_index);
		bool child_state  = current_state.get(index_in_tree);
		return process.probability(binned_branch_length, parent_state, child_state);
	}
};

bool sample(std::array<coretools::TSumLogProbability, 2> &sum_log);
bool sample(double log_prob_0, double log_prob_1);

#endif // ACOL_TCLIQUE_H
