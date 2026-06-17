#include "../TTree.h"

std::vector<TClique> &TTree::get_cliques() { return _cliques; }

const TClique &TTree::get_clique(const std::vector<size_t> &index_in_leaves_space) const {
	size_t ix_clique = 0;
	size_t stride    = 1;

	for (size_t i = 0; i < _dimension_cliques.size(); ++i) {
		const size_t idx = (i == _dimension) ? 0 : index_in_leaves_space[i];
		ix_clique += idx * stride;
		stride *= _dimension_cliques[i];
	}

	return _cliques[ix_clique];
}

TClique &TTree::get_clique(const std::vector<size_t> &index_in_leaves_space) {
	size_t ix_clique = 0;
	size_t stride    = 1;

	for (size_t i = 0; i < _dimension_cliques.size(); ++i) {
		const size_t idx = (i == _dimension) ? 0 : index_in_leaves_space[i];
		ix_clique += idx * stride;
		stride *= _dimension_cliques[i];
	}

	return _cliques[ix_clique];
}

void TTree::_initialize_cliques(const std::vector<size_t> &num_leaves_per_tree,
                                const std::vector<std::unique_ptr<TTree>> &all_trees) {
	// clique of a tree: runs along that dimension
	// the cliques of a tree are can only contain leaves in all trees except the one we are working
	// on.
	_dimension_cliques             = num_leaves_per_tree;
	_dimension_cliques[_dimension] = 1;

	// we then caclulate how many cliques we will have in total for that tree. Which is the product
	// of the number of leaves in each tree except the one we are working on (that is why we set it
	// to 1 before).
	const size_t n_cliques = coretools::containerProduct(_dimension_cliques);

	// calculate increment: product of the number of leaves of all subsequent dimensions
	size_t increment = 1;
	for (size_t i = _dimension + 1; i < all_trees.size(); ++i) {
		increment *= all_trees[i]->get_number_of_leaves();
	}

	// initialize cliques
	for (size_t i = 0; i < n_cliques; ++i) {
		// get start index of each clique in leaves space
		std::vector<size_t> start_index_in_leaves_space =
		    coretools::getSubscripts(i, _dimension_cliques);
		_cliques.emplace_back(start_index_in_leaves_space, _dimension, _nodes.size(), increment);
		_cliques.back().initialize(_delta, _number_of_bins);

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

void TTree::_simulation_prepare_cliques(size_t c, TClique &clique) const {
	clique.initialize(this->get_delta(), this->get_number_of_bins());
	clique.set_lambda(_alpha_c->value(c), _nu_c[c]);
};
