//
// Created by madleina on 22.01.25.
//

#ifndef ACOL_TCOLLAPSER_H
#define ACOL_TCOLLAPSER_H

#include "TTree.h"
#include <string>
#include <vector>

class TCollapser {
private:
	// trees
	const std::vector<std::unique_ptr<TTree>> &_trees;

	// dimensions
	std::vector<size_t> _dimensions_to_keep;
	std::vector<size_t> _dimensions_to_collapse;

	template<bool ModifyCounter>
	bool _x_is_one_in_clique(size_t dim_along_which_clique_runs, const std::vector<size_t> &index_in_leaves) const {
		// count the number of leaves with state = 1 in current clique (corresponds to old_state)
		const auto c = _trees[dim_along_which_clique_runs]->get_clique(index_in_leaves).get_counter_leaves_state_1();

		// new_state is always zero once we get here
		// 0 and 0 -> c does not change
		if constexpr (!ModifyCounter) { return c > 0; }

		// new_state is zero and old_state is one
		// -> depends on the states of the other Y's in the clique
		// -> if old_state was the only 1 in the clique -> new_state will not be one anymore
		return ((int)c - 1) > 0;
	}

	template<bool ModifyCounter> bool _x_is_one(std::vector<size_t> index_in_leaves) const {
		assert(do_collapse());

		// define dimension along which we take the cliques: always the last dimension to collapse
		const size_t dim_along_which_clique_runs     = _dimensions_to_collapse.back();
		index_in_leaves[dim_along_which_clique_runs] = 0;

		if (_dimensions_to_collapse.size() == 1) {
			// only need to consider a single clique
			return _x_is_one_in_clique<ModifyCounter>(dim_along_which_clique_runs, index_in_leaves);
		}

		// define sub-space of dimensions to get the starting positions of all cliques to consider
		// the last dimension to collapse is used as a clique -> no need to consider that one -> take size() - 1
		std::vector<size_t> dimensions_clique_starts(_dimensions_to_collapse.size() - 1);
		for (size_t i = 0; i < _dimensions_to_collapse.size() - 1; ++i) {
			const size_t tree_index     = _dimensions_to_collapse[i];
			dimensions_clique_starts[i] = _trees[tree_index]->get_number_of_leaves();
		}

		const size_t num_loops = coretools::containerProduct(dimensions_clique_starts);
		for (size_t i = 0; i < num_loops; ++i) { // loop over linear index of clique start
			const auto ix = coretools::getSubscripts(i, dimensions_clique_starts);
			// set clique starting index
			for (size_t j = 0; j < _dimensions_to_collapse.size() - 1; ++j) {
				const size_t tree_index     = _dimensions_to_collapse[j];
				index_in_leaves[tree_index] = ix[j];
			}
			// if at least one clique has a one: x is one
			if (_x_is_one_in_clique<ModifyCounter>(dim_along_which_clique_runs, index_in_leaves)) { return true; }
		}

		// we checked all the cliques -> never a one -> return false
		return false;
	}

public:
	TCollapser(const std::vector<std::unique_ptr<TTree>> &trees);
	~TCollapser() = default;

	std::vector<size_t> initialize(const std::vector<std::string> &dimension_names_to_keep, std::string_view data_name);

	// getters
	bool x_is_one(const std::vector<size_t> &index_in_leaves, bool old_state) const;
	bool x_is_one(const std::vector<size_t> &index_in_leaves) const;

	std::vector<size_t> collapse(const std::vector<size_t> &index_in_full_space) const;
	std::vector<size_t> expand(const std::vector<size_t> &index_in_collapsed_space) const;
	bool do_collapse() const;
	size_t num_dim_to_keep() const;
	size_t num_dim_to_collapse() const;
	size_t dim_to_keep(size_t i) const;
	size_t dim_to_collapse(size_t i) const;
};

#endif // ACOL_TCOLLAPSER_H
