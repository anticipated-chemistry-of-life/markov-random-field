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
#include "tree/clique/TCliqueView.h"
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

class TTree;

/** Class representing a clique in our model. A clique is defined as having a set of nodes that are
 * all leaves in all dimensions except one. Each clique has a transition grid and its own
 * multidimensional index. Turning that index into the cell a node occupies is TCliqueView's, and
 * nothing here does it.
 */
class TClique {
private:
	/// This clique's two-state process, discretised onto the tree's bin grid. Set once the
	/// parameters exist (TTree::guessInitialValues) and replaced wholesale whenever a proposal on
	/// alpha or nu is accepted; there is no mutable "try" copy.
	std::optional<TTransitionGrid> _transition_grid;

	// info about size and dimensionality of clique
	IndexArray _start_index_in_leaves_space;

	/// @brief Calculates the log probability of the root given the stationary distribution
	static void _calculate_log_prob_root(double stationary_0,
	                                     std::array<coretools::TSumLogProbability, 2> &sum_log);

	/// @brief Calculates the log probability of a node to its children
	void _calculate_log_prob_node_to_children(
	    size_t index_in_tree, const TTree *tree, const TNodeStateCliqueView &states,
	    std::array<coretools::TSumLogProbability, 2> &sum_log) const;

	/// @brief Starts one internal node at the state its children make most likely. This is
	/// initialisation and not a sampler move: it runs once, before the chain's first update, and
	/// it takes the mode rather than a draw.
	/// @param node_index The index of the internal node we want to start
	/// @param states This clique's cells of the node state, read and written through its view.
	/// @param tree the tree of interest
	void _initialize_node_from_children(size_t node_index, TNodeStateCliqueView &states,
	                                    const TTree *tree) const;

	static size_t _get_parent_index(size_t index_in_tree, const TTree *tree);

public:
	explicit TClique(const IndexArray &start_index);
	~TClique() = default;

	/// @brief Install this clique's process. Called once the parameters exist, and again whenever a
	/// proposal on alpha or nu is accepted.
	void set_transition_grid(TTransitionGrid grid) { _transition_grid = std::move(grid); }

	/// @brief This clique's current process. Throws if the parameters have not been drawn yet,
	/// which used to read as a grid of zeros instead.
	[[nodiscard]] const TTransitionGrid &transition_grid() const {
		return _transition_grid.value();
	}

	/// This clique's own multidimensional index: a leaf in every dimension but its tree's own,
	/// which carries a 0. Setting that dimension to a node index gives that node's cell, which is
	/// what TCliqueView does and what nothing else does.
	[[nodiscard]] const IndexArray &clique_index() const { return _start_index_in_leaves_space; }

	/// @brief Update the Z dimension for this clique.
	/// @param states This clique's cells of the node state, read and written through its view.
	/// @param tree The tree.
	/// @param uniforms the node state's stream for this iteration. Each node draws the one uniform
	/// its own cell names, so the walk gives the same states whichever thread runs it.
	///
	/// The walk keeps no running density. It used to add each drawn node's own log probability,
	/// which scored that node against its parent *and* against every child, so each internal edge
	/// counted twice. The joint density is a question about the configuration the walk leaves
	/// behind, and tree/node_state_density.h answers it there.
	void update_Z(TNodeStateCliqueView &states, const TTree *tree,
	              const TCellUniforms &uniforms) const;

	void initialize_Z_from_children(TNodeStateCliqueView &states, const TTree *tree) const;

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

	/// @brief P(node | parent) under an explicitly given process, so a Metropolis proposal can ask
	/// the same question of the current grid and of its candidate.
	double calculate_prob_to_parent(size_t index_in_tree, const TTree *tree,
	                                TypeBinnedBranchLengths binned_branch_length,
	                                const TNodeStateCliqueView &states,
	                                const TTransitionGrid &process) const {
		size_t parent_index = _get_parent_index(index_in_tree, tree);

		bool parent_state = states.is_one(parent_index);
		bool child_state  = states.is_one(index_in_tree);
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
