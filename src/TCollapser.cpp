//
// Created by madleina on 22.01.25.
//

#include "TCollapser.h"
#include "coretools/Main/TLog.h"

TCollapser::TCollapser(const std::vector<TTree> &trees) : _trees(trees) {}

std::vector<size_t> TCollapser::initialize(const std::vector<std::string> &dimension_names_to_keep,
                                           std::string_view data_name) {
	using namespace coretools::instances;
	// all dimensions that are not present in header will be collapsed
	if (dimension_names_to_keep.size() > _trees.size()) {
		UERROR(data_name, " can not have more dimensions than there are trees (", dimension_names_to_keep.size(),
		       " vs ", _trees.size(), ")");
	}

	std::vector<size_t> len_per_dimension;
	for (const auto &tree_name : dimension_names_to_keep) {
		bool found = false;
		for (size_t j = 0; j < _trees.size(); ++j) {
			if (_trees[j].get_tree_name() == tree_name) {
				// already found before -> duplicated!
				if (found) { UERROR("Duplicate column name '", tree_name, "' in ", data_name, " file!"); }
				// else: remember this dimension -> we will not collapse it
				_dimensions_to_keep.push_back(j);
				len_per_dimension.push_back(_trees[j].get_number_of_leaves());
				found = true;
			}
		}
		if (!found) {
			UERROR("Could not find tree with name '", tree_name, "' in trees (required by ", data_name, ").");
		}
	}

	if (_dimensions_to_keep.empty()) { UERROR("No dimensions in ", data_name, " file are kept!"); }
	if (_dimensions_to_keep.back() != _trees.size() - 1) {
		UERROR("Last dimension of trees and ", data_name, " must be identical (", _dimensions_to_keep.back(), " vs ",
		       _trees.size() - 1, ")!");
	}

	// find dimensions to collapse
	for (size_t i = 0; i < _trees.size(); ++i) {
		if (std::find(_dimensions_to_keep.begin(), _dimensions_to_keep.end(), i) == _dimensions_to_keep.end()) {
			// not found in keep -> collapse
			_dimensions_to_collapse.emplace_back(i);
		}
	}

	// report to logfile
	logfile().startIndent("Will keep the following dimensions for ", data_name, ": ", _dimensions_to_keep, ":");
	for (const auto i : _dimensions_to_keep) { logfile().list(_trees[i].get_tree_name()); }
	logfile().endIndent();
	if (_dimensions_to_collapse.empty()) {
		logfile().list("Will not collapse any dimensions for ", data_name, ".");
	} else {
		logfile().startIndent("Will collapse the following dimensions for ", data_name, ": ", _dimensions_to_collapse,
		                      ":");
		for (const auto i : _dimensions_to_collapse) { logfile().list(_trees[i].get_tree_name()); }
		logfile().endIndent();
	}

	return len_per_dimension;
}

bool TCollapser::_x_is_one(size_t dim_along_which_clique_runs, const std::vector<size_t> &index_in_leaves,
                           bool old_state) const {
	// count the number of leaves with state = 1 in current clique (corresponds to old_state)
	const auto c = _trees[dim_along_which_clique_runs].get_clique(index_in_leaves).get_counter_leaves_state_1();

	// new_state is always zero once we get here
	// 0 and 0 -> c does not change
	if (!old_state) { return c > 0; }

	// new_state is zero and old_state is one
	// -> depends on the states of the other Y's in the clique
	// -> if old_state was the only 1 in the clique -> new_state will not be one anymore
	return ((int)c - 1) > 0;
}

bool TCollapser::x_is_one(std::vector<size_t> index_in_leaves, bool new_state, bool old_state) const {
	// we return x_is_one for the new_state
	if (!do_collapse()) { return new_state; }
	if (new_state) { return true; }

	// define dimension along which we take the cliques: always the last dimension to collapse
	const size_t dim_along_which_clique_runs     = _dimensions_to_collapse.back();
	index_in_leaves[dim_along_which_clique_runs] = 0;

	if (_dimensions_to_collapse.size() == 1) {
		// only need to consider a single clique
		return _x_is_one(dim_along_which_clique_runs, index_in_leaves, old_state);
	}

	// define sub-space of dimensions to get the starting positions of all cliques to consider
	// the last dimension to collapse is used as a clique -> no need to consider that one -> take size() - 1
	std::vector<size_t> dimensions_clique_starts(_dimensions_to_collapse.size() - 1);
	for (size_t i = 0; i < _dimensions_to_collapse.size() - 1; ++i) {
		const size_t tree_index     = _dimensions_to_collapse[i];
		dimensions_clique_starts[i] = _trees[tree_index].get_number_of_leaves();
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
		if (_x_is_one(dim_along_which_clique_runs, index_in_leaves, old_state)) { return true; }
	}

	// we checked all the cliques -> never a one -> return false
	return false;
}

bool TCollapser::x_is_one(std::vector<size_t> index_in_leaves, bool state) const {
	return x_is_one(index_in_leaves, state, state);
}

std::vector<size_t> TCollapser::collapse(const std::vector<size_t> &index_in_full_space) const {
	std::vector<size_t> index_in_collapsed_space(_dimensions_to_keep.size());
	for (size_t i = 0; i < _dimensions_to_keep.size(); ++i) {
		index_in_collapsed_space[i] = index_in_full_space[_dimensions_to_keep[i]];
	}

	return index_in_collapsed_space;
}

std::vector<size_t> TCollapser::expand(const std::vector<size_t> &index_in_collapsed_space) const {
	std::vector<size_t> index_in_leaves(_trees.size(), 0);
	for (size_t i = 0; i < _dimensions_to_keep.size(); ++i) {
		index_in_leaves[_dimensions_to_keep[i]] = index_in_collapsed_space[i];
	}

	return index_in_leaves;
}

bool TCollapser::do_collapse() const { return !_dimensions_to_collapse.empty(); }
size_t TCollapser::num_dim_to_keep() const { return _dimensions_to_keep.size(); }
size_t TCollapser::num_dim_to_collapse() const { return _dimensions_to_collapse.size(); }
size_t TCollapser::dim_to_keep(size_t i) const { return _dimensions_to_keep[i]; }
size_t TCollapser::dim_to_collapse(size_t i) const { return _dimensions_to_collapse[i]; }