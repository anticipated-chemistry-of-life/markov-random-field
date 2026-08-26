//
// Created by madleina on 22.10.24.
//

#ifndef ACOL_TCLIQUE_H
#define ACOL_TCLIQUE_H

#include "TCurrentState.h"
#include "Types.h"
#include "constants.h"
#include "coretools/Math/TSumLog.h"
#include "storages/y_storage/TStorageYMatrix.h"
#include "storages/z_storage/TStorageZMatrix.h"
#include "tree/branch/TTransitionGrid.h"
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

class TTree;

/** Class representing a clique in our model. A clique is defined as having a set of nodes that are
 * all leaves in all dimensions except one. Each clique has a transition grid, and the start index
 * of the nodes in the tree. The start index, the variable dimension, and the number of nodes are
 * needed to get the correct indices in our multidimensional space Y and Z.
 */
class TClique {
private:
	using TypeParamBinBranches = stattools::TParameter<SpecBinnedBranches, TTree>;

	/// This clique's two-state process, discretised onto the tree's bin grid. Set once the
	/// parameters exist (TTree::guessInitialValues) and replaced wholesale whenever a proposal on
	/// alpha or nu is accepted; there is no mutable "try" copy.
	std::optional<TTransitionGrid> _transition_grid;

	// info about size and dimensionality of clique
	IndexArray _start_index_in_leaves_space;
	size_t _variable_dimension;
	size_t _n_nodes;
	size_t _increment;

	// count the number of leaves with value 1 in a clique
	TypeCounter1 _counter_leaves_state_1 = 0;

	/// @brief Calculates the log probability of the root given the stationary distribution
	static void _calculate_log_prob_root(double stationary_0,
	                                     std::array<coretools::TSumLogProbability, 2> &sum_log);

	template<typename ContainerStates> // can either be TSheet or TCurrentStates
	inline bool _getState(const ContainerStates &states, size_t parent_index_in_tree,
	                      size_t leaf_index_in_tree_of_last_dim) const {
		if constexpr (std::is_same_v<ContainerStates, TSheet>) { // is a sheet
			return states.get(parent_index_in_tree, leaf_index_in_tree_of_last_dim);
		} else { // TCurrentState
			return states.get(parent_index_in_tree);
		}
	}

	void _update_current_state(TStorageZMatrix &Z, TCurrentState &current_state,
	                           size_t index_in_tree, bool new_state,
	                           std::vector<size_t> &linear_indices_in_Z_space_to_insert,
	                           const TTree *tree) const;

	/// @brief Calculates the log probability of a node to its children
	void _calculate_log_prob_node_to_children(
	    size_t index_in_tree, const TTree *tree, const TCurrentState &current_state,
	    std::array<coretools::TSumLogProbability, 2> &sum_log,
	    const TypeParamBinBranches *binned_branch_lengths,
	    const std::vector<size_t> &leaves_and_internal_nodes_without_roots_indices) const;

	/// @brief Sets Z given the maximal likelihood given its children. This was created to avoid
	/// that Z is stuck in a state and cannot change.
	/// @param node_index The index of the internal node we want to set
	/// @param current_state The current state of the clique.
	/// @param Z the Z vector of that tree (i.e that clique)
	/// @param tree the tree of interest
	/// @param binned_branch_lengths the vector of branch length
	/// @param leaves_and_internal_nodes_without_roots_indices Same as the variable name
	/// @param linear_indices_in_Z_space_to_insert Same as the variable name
	void _set_Z_to_MLE(size_t node_index, TCurrentState &current_state, TStorageZMatrix &Z,
	                   const TTree *tree, const TypeParamBinBranches *binned_branch_lengths,
	                   const std::vector<size_t> &leaves_and_internal_nodes_without_roots_indices,
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
	std::vector<size_t>
	update_Z(std::vector<double> &joint_prob_density, TCurrentState &current_state,
	         TStorageZMatrix &Z, const TTree *tree,
	         const TypeParamBinBranches *binned_branch_lengths,
	         const std::vector<size_t> &leaves_and_internal_nodes_without_roots_indices) const;

	std::vector<size_t> initialize_Z_from_children(
	    TCurrentState &current_state, TStorageZMatrix &Z, const TTree *tree,
	    const TypeParamBinBranches *binned_branch_lengths,
	    const std::vector<size_t> &leaves_and_internal_nodes_without_roots_indices) const;

	TCurrentState create_current_state(const TStorageYMatrix &Y, const TStorageZMatrix &Z,
	                                   const TTree &tree);

	/// @brief Return the number of nodes in the clique
	/// @return Return the number of nodes in the clique
	[[nodiscard]] size_t get_number_of_nodes() const { return _n_nodes; }

	/// @brief Calculates the log probability of a node to its parent, under this clique's current
	/// process.
	template<typename ContainerStates> // ContainerStates can either be TSheet or TCurrentStates
	void calculate_log_prob_parent_to_node(
	    size_t index_in_tree, TypeBinnedBranchLengths binned_branch_length, const TTree *tree,
	    size_t leaf_index_in_tree_of_last_dim, const ContainerStates &states,
	    std::array<coretools::TSumLogProbability, 2> &sum_log) const {
		const size_t parent_index_in_tree = _get_parent_index(index_in_tree, tree);
		const auto &process               = transition_grid();
		for (size_t i = 0; i < 2; ++i) { // loop over possible values (0 or 1) of the node
			const bool state_of_parent =
			    _getState(states, parent_index_in_tree, leaf_index_in_tree_of_last_dim);
			sum_log[i].add(process.probability(binned_branch_length, state_of_parent, i));
		}
	}

	/// Updates the counter of the clique. This is used by the collapser.
	void update_counter_leaves_state_1(bool new_state, bool old_state);
	void update_counter_leaves_state_1(int difference);
	[[nodiscard]] size_t get_counter_leaves_state_1() const { return _counter_leaves_state_1; }

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
	                                const TCurrentState &current_state,
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
