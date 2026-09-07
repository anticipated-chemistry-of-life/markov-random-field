#include "../TTree.h"
#include "constants.h"
#include "coretools/algorithms.h"
#include <utility>

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
	_n_cliques = coretools::containerProduct(_dimension_cliques);

	// initialize cliques
	for (size_t i = 0; i < _n_cliques; ++i) {
		// get start index of each clique in leaves space
		auto start_index_in_leaves_space = coretools::getSubscriptsAsArray(i, _dimension_cliques);

		// build clique name from leaf names in all other dimensions
		std::string name;
		for (size_t d = 0; d < all_trees.size(); ++d) {
			if (d == _dimension) continue;
			size_t node_idx =
			    all_trees[d]->get_node_index_from_leaf_index(start_index_in_leaves_space[d]);
			if (!name.empty()) name += "_";
			name += all_trees[d]->get_node_id(node_idx);
		}
		_clique_names.push_back(name);
	}
}

void TTree::_simulation_prepare_cliques(size_t c) {
	_transition_grid_per_clique[c] = TTransitionGrid(_alpha_c->value(c), _nu_c[c], _grid());
};
