//
// The clique side of a tree: which grid a clique carries, where its cells are, and the walk that
// gives its nodes their states.
//

#include "../TTree.h"
#include "constants.h"
#include "coretools/Math/TSumLog.h"
#include "coretools/algorithms.h"
#include "random/two_state_draw.h"

#include <array>
#include <cstddef>

IndexArray TTree::_clique_index(size_t c) const {
	// The 0 this tree's own dimension carries comes from the 1 _dimension_cliques holds there.
	// Setting that dimension to a node index gives that node's cell, which TCliqueView does and
	// nothing else does.
	return coretools::getSubscriptsAsArray(c, _dimension_cliques);
}

const TTransitionGrid &
TTree::transition_grid_of_cell(const IndexArray &index_in_leaves_space) const {
	size_t ix_clique = 0;
	size_t stride    = 1;

	for (size_t i = 0; i < _dimension_cliques.size(); ++i) {
		const size_t idx = (i == _dimension) ? 0 : index_in_leaves_space[i];
		ix_clique += idx * stride;
		stride *= _dimension_cliques[i];
	}

	return transition_grid(ix_clique);
}

void TTree::_initialize_cliques(const IndexArray &num_leaves_per_tree,
                                const std::vector<std::unique_ptr<TTree>> &all_trees) {
	// clique of a tree: runs along that dimension
	// the cliques of a tree are can only contain leaves in all trees except the one we are working
	// on.
	_dimension_cliques             = num_leaves_per_tree;
	_dimension_cliques[_dimension] = 1;

	// we then caclulate how many cliques we will have in total for that tree. Which is the product
	// of the number of leaves in each tree except the one we are working on (that is why we set it
	// to 1 before).
	const size_t clique_count = coretools::containerProduct(_dimension_cliques);

	// One grid slot per clique, in clique order, and all of them empty. A grid needs alpha and nu,
	// which stattools has not drawn yet. TTree::guessInitialValues installs them.
	_transition_grids.resize(clique_count);

	for (size_t i = 0; i < clique_count; ++i) {
		const IndexArray clique_index = _clique_index(i);

		// build clique name from leaf names in all other dimensions
		std::string name;
		for (size_t d = 0; d < all_trees.size(); ++d) {
			if (d == _dimension) continue;
			size_t node_idx = all_trees[d]->get_node_index_from_leaf_index(clique_index[d]);
			if (!name.empty()) name += "_";
			name += all_trees[d]->get_node_id(node_idx);
		}
		_clique_names.push_back(name);
	}
}

/// Gives every internal node of clique `c` a new state, one node at a time.
///
/// Each node draws the one uniform its own cell names, so the walk gives the same states whichever
/// thread runs it.
///
/// The walk keeps no running density. It used to add each drawn node's own log probability, which
/// scored that node against its parent *and* against every child, so each internal edge counted
/// twice. The joint density is a question about the configuration the walk leaves behind, and
/// tree/node_state_density.h answers it there.
void TTree::_update_Z_of_clique(size_t c, TNodeStateCliqueView &states,
                                const TCellUniforms &uniforms) const {
	const TTransitionGrid &process = transition_grid(c);
	const double stationary_0      = process.stationary(false);

	for (const auto index_in_tree : get_internal_nodes()) {
		// prepare log probabilities for the two possible states
		std::array<coretools::TSumLogProbability, 2> sum_log;
		if (_topology().is_root(index_in_tree)) { // calculate stationary
			_log_prob_root(stationary_0, sum_log);
		} else { // calculate P(node = 0 | parent) and P(node = 1 | parent)
			// Note: the *previous* bin, because branch lengths are proposed before the loop starts.
			// The parent comes after this node in post-order, so the walk has not touched it and
			// its state is still the one this update started from.
			const auto bin_branch_len  = get_previous_binned_branch_length(index_in_tree);
			const bool state_of_parent = states.is_one(_topology().parent_of(index_in_tree));
			for (size_t i = 0; i < 2; ++i) { // loop over possible values (0 or 1) of the node
				sum_log[i].add(process.probability(bin_branch_len, state_of_parent, i));
			}
		}

		// calculate P(child | node = 0) and P(child | node = 1) for all children of node
		_log_prob_node_to_children(index_in_tree, process, states, sum_log);

		// sample new state and update Z accordingly. The cell decides which uniform it draws, so
		// the state this node gets does not depend on which thread walked this clique.
		const double log_prob_0 = sum_log[0].getSum();
		const double log_prob_1 = sum_log[1].getSum();
		const bool new_state    = two_state_draw::sample(
		    log_prob_0, log_prob_1, uniforms.at(states.linear_index(index_in_tree)));

		// The view writes the cell the node state already holds in place. It defers the one the
		// node state cannot take. See ADR-0006. A later read on this view sees the new state
		// either way, which is what the parent of this node needs.
		states.set_state(index_in_tree, new_state);
	}
}

void TTree::_initialize_clique_from_children(size_t c, TNodeStateCliqueView &states) const {
	// Bottom-up start of Z, as one forward walk. The internal nodes are stored as the non-root
	// block in post-order followed by the roots (ADR-0004), so every node's children are already
	// done by the time it comes up -- leaves before all of them, and each parent after its own
	// children.
	const TTransitionGrid &process = transition_grid(c);
	for (const size_t node_index : get_internal_nodes()) {
		_initialize_node_from_children(node_index, process, states);
	}
}

/// Starts one internal node at the state its children make most likely. This is initialisation and
/// not a sampler move: it runs once, before the chain's first update, and it takes the mode rather
/// than a draw.
void TTree::_initialize_node_from_children(size_t node_index, const TTransitionGrid &process,
                                           TNodeStateCliqueView &states) const {
	std::array<coretools::TSumLogProbability, 2> sum_log;

	_log_prob_node_to_children(node_index, process, states, sum_log);

	const double log_prob_0 = sum_log[0].getSum();
	const double log_prob_1 = sum_log[1].getSum();

	// The mode, not a draw: this is where the chain starts, and the first update moves it.
	const bool most_likely_state = log_prob_1 > log_prob_0;
	states.set_state(node_index, most_likely_state);
}

/// The log probability of the root under the stationary distribution.
void TTree::_log_prob_root(double stationary_0,
                           std::array<coretools::TSumLogProbability, 2> &sum_log) {
	sum_log[0].add(stationary_0);
	sum_log[1].add(1.0 - stationary_0);
}

/// The log probability of a node to its children, under this clique's process.
void TTree::_log_prob_node_to_children(
    size_t index_in_tree, const TTransitionGrid &process, const TNodeStateCliqueView &states,
    std::array<coretools::TSumLogProbability, 2> &sum_log) const {
	for (const auto &child_index : children_of(index_in_tree)) {
		// Note: the *previous* bin, because new values were proposed before the loop started.
		// Children come before their parent in post-order, so this reads a state the walk has
		// already assigned -- which the view shows even where it could not write it in place.
		auto bin_length        = get_previous_binned_branch_length(child_index);
		const bool child_state = states.is_one(child_index);
		for (size_t i = 0; i < 2; ++i) { // loop over possible values (0 or 1) of the node
			sum_log[i].add(process.probability(bin_length, i, child_state));
		}
	}
}
