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

bool TCollapser::_ix_is_one_no_collapse(bool state) { return state; }

std::pair<bool, bool> TCollapser::x_is_one_new_old(std::vector<size_t> index_in_leaves, bool new_state,
                                                   bool old_state) const {
	if (!do_collapse()) { return {new_state, old_state}; }


}

bool TCollapser::x_is_one(std::vector<size_t> index_in_leaves, bool state) const {}

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