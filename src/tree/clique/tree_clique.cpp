//
// The clique side of a tree: which grid a clique carries, where its cells are, and the bottom-up
// start that gives its nodes their first states.
//

#include "../TTree.h"
#include "constants.h"
#include "coretools/Math/TSumLog.h"
#include "coretools/algorithms.h"
#include "tree/node_state_walk.h"

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

	// The same child terms the node-state walk adds, and from the same place. The start reads them
	// alone: it has no parent term, because it takes the state the children make most likely.
	node_state_walk::add_log_prob_of_children(_topology(), process, _previous_bins(), states,
	                                          node_index, sum_log);

	const double log_prob_0 = sum_log[0].getSum();
	const double log_prob_1 = sum_log[1].getSum();

	// The mode, not a draw: this is where the chain starts, and the first update moves it.
	const bool most_likely_state = log_prob_1 > log_prob_0;
	states.set_state(node_index, most_likely_state);
}
