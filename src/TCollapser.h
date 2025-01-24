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
	const std::vector<TTree> &_trees;

	// dimensions
	std::vector<size_t> _dimensions_to_keep;
	std::vector<size_t> _dimensions_to_collapse;

	bool _x_is_one(size_t dim_along_which_clique_runs, const std::vector<size_t> &index_in_leaves,
	               bool old_state) const;

public:
	TCollapser(const std::vector<TTree> &trees);
	~TCollapser() = default;

	std::vector<size_t> initialize(const std::vector<std::string> &dimension_names_to_keep, std::string_view data_name);

	// getters
	bool x_is_one(std::vector<size_t> index_in_leaves, bool new_state, bool old_state) const;
	bool x_is_one(std::vector<size_t> index_in_leaves) const;

	std::vector<size_t> collapse(const std::vector<size_t> &index_in_full_space) const;
	std::vector<size_t> expand(const std::vector<size_t> &index_in_collapsed_space) const;
	bool do_collapse() const;
	size_t num_dim_to_keep() const;
	size_t num_dim_to_collapse() const;
	size_t dim_to_keep(size_t i) const;
	size_t dim_to_collapse(size_t i) const;
};

#endif // ACOL_TCOLLAPSER_H
