//
// Created by madleina on 22.10.24.
//

#ifndef ACOL_TCLIQUE_H
#define ACOL_TCLIQUE_H

#include "Types.h"
#include "constants.h"
#include "coretools/Main/TError.h"
#include "coretools/Math/TSumLog.h"
#include "coretools/Types/probability.h"
#include "random/TCellUniforms.h"
#include "storages/storage_backend.h"
#include "tree/TPhylogeny.h"
#include "tree/branch/TTransitionGrid.h"
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

class TTree;
class TClique;

/// The states of one clique's nodes, read and written through a window.
///
/// The window runs over every node of the clique, leaves included. Nothing here reads the field.
/// A tree field is the leaf block of this very run of cells, so the walk below, and the alpha and
/// nu moves after it, read one tree alone. See ADR-0005.
///
/// The walk keeps no copy of the states it assigns. It reads them back from the window it wrote
/// them through. The window stays open across the moves that follow the walk, which read those
/// states too. ADR-0006 gives the argument for both.
class TCliqueStates {
private:
	/// This clique's cells of the node state: every node of the tree.
	TNodeStateStorage::TWindow _nodes;
	const TPhylogeny *_topology;

public:
	/// Opens the window over this clique's cells of the node state.
	TCliqueStates(TNodeStateStorage &Z, const TClique &clique, const TPhylogeny &topology);

	/// The state of `node`, from this tree's own node state.
	[[nodiscard]] bool get(size_t node) const { return _nodes.is_one(node); }

	/// Gives `node` its new state. A write of the state the cell already carries is dropped. The
	/// sparse window would otherwise buffer an insert for a cell it does not hold, and then throw
	/// that insert away.
	void set(size_t node, bool state) {
		// The walk assigns internal nodes. A leaf's state is drawn with the field and the other
		// tree's leaf, as one eight-state block, and not here (ADR-0005).
		DEBUG_ASSERT(!_topology->is_leaf(node));
		if (_nodes.is_one(node) != state) { _nodes.set_state(node, state); }
	}

	/// The linear index, in the node state's container space, of the cell `node` occupies. This is
	/// the cell's name to the stream of uniforms.
	[[nodiscard]] size_t linear_index_in_Z(size_t node) const { return _nodes.linear_index(node); }

	/// Hands out the cells the window could not write in place, as linear indices, and ends the
	/// window. The caller commits them after the parallel region. That is the only exit a window
	/// inside one may take (ADR-0006).
	[[nodiscard]] std::vector<size_t> take_buffered_inserts() {
		return _nodes.take_buffered_inserts();
	}
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

	/// @brief Calculates the log probability of a node to its children
	void _calculate_log_prob_node_to_children(
	    size_t index_in_tree, const TTree *tree, const TCliqueStates &states,
	    std::array<coretools::TSumLogProbability, 2> &sum_log) const;

	/// @brief Starts one internal node at the state its children make most likely. This is
	/// initialisation and not a sampler move: it runs once, before the chain's first update, and
	/// it takes the mode rather than a draw.
	/// @param node_index The index of the internal node we want to start
	/// @param states The states of this clique's nodes, read and written through its windows.
	/// @param tree the tree of interest
	void _initialize_node_from_children(size_t node_index, TCliqueStates &states,
	                                    const TTree *tree) const;

	static size_t _get_parent_index(size_t index_in_tree, const TTree *tree);

	/// The cell node `node` of this clique occupies, in the column the clique runs at. Setting the
	/// last dimension first and the variable one second is what makes this right for a clique along
	/// either dimension: when they are the same dimension, the node wins.
	[[nodiscard]] IndexArray cell_of(size_t node, size_t leaf_index_in_tree_of_last_dim) const {
		IndexArray cell           = _start_index_in_leaves_space;
		cell.back()               = leaf_index_in_tree_of_last_dim;
		cell[_variable_dimension] = node;
		return cell;
	}

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

	/// The window this clique's turn reads and writes through: the node state's column. The caller
	/// keeps it open for the whole of that turn, because the parameter and branch-length moves that
	/// follow the node-state walk have to see the states that walk assigned.
	[[nodiscard]] TCliqueStates open_states(TNodeStateStorage &Z,
	                                        const TPhylogeny &topology) const;

	/// The node state's column for this clique, on its own. The simulation draws every node from
	/// its parent and reads no leaf, so it needs no window on the field.
	[[nodiscard]] TNodeStateStorage::TWindow open_node_state_window(TNodeStateStorage &Z) const;

	/// The cell this clique's first node occupies. A window over the clique opens here and steps
	/// by the increment.
	[[nodiscard]] IndexArray first_cell() const {
		return cell_of(0, _start_index_in_leaves_space.back());
	}

	/// @brief Update the Z dimension for this clique.
	/// @param states The states of this clique's nodes, read and written through its windows.
	/// @param tree The tree.
	/// @param uniforms the node state's stream for this iteration. Each node draws the one uniform
	/// its own cell names, so the walk gives the same states whichever thread runs it.
	void update_Z(std::vector<double> &joint_prob_density, TCliqueStates &states, const TTree *tree,
	              const TCellUniforms &uniforms) const;

	void initialize_Z_from_children(TCliqueStates &states, const TTree *tree) const;

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

	/// @brief P(node | parent) under an explicitly given process, so a Metropolis proposal can ask
	/// the same question of the current grid and of its candidate.
	double calculate_prob_to_parent(size_t index_in_tree, const TTree *tree,
	                                TypeBinnedBranchLengths binned_branch_length,
	                                const TCliqueStates &states,
	                                const TTransitionGrid &process) const {
		size_t parent_index = _get_parent_index(index_in_tree, tree);

		bool parent_state = states.get(parent_index);
		bool child_state  = states.get(index_in_tree);
		return process.probability(binned_branch_length, parent_state, child_state);
	}
};

/// The two-state draw: state 1 with the probability its caller names, either as a probability or
/// as a pair of log probabilities.
///
/// Every caller supplies the uniform, and it comes from the cell being drawn (`TCellUniforms`)
/// rather than from a running generator. A caller that cannot name its cell therefore cannot draw.
/// ADR-0007 says why.
bool sample(coretools::Probability probability_of_one, double uniform);
bool sample(std::array<coretools::TSumLogProbability, 2> &sum_log, double uniform);
bool sample(double log_prob_0, double log_prob_1, double uniform);

#endif // ACOL_TCLIQUE_H
