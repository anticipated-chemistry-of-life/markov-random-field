#pragma once

#include "TStorageYVector.h"
#include "TTree.h"
#include <cstddef>
#include <string>
#include <vector>

class TLinks {
private:
	// trees should be a const ref because we don't want to change the trees and don't want to copy them
	const std::vector<TTree> &_trees;
	TStorageYVector _data_Y;

	// private functions
	size_t _tree_index_of_node(const std::string &node_id) const {
		for (size_t i = 0; i < _trees.size(); ++i) {
			if (_trees[i].in_tree(node_id)) { return i; }
		}
		UERROR("Node '", node_id, "' is not in any of the trees !");
	}

	static std::vector<size_t> _get_dimensions_Y_space(const std::vector<TTree> &trees) {
		std::vector<size_t> dimensions_Y_space(trees.size());
		for (size_t i = 0; i < trees.size(); ++i) { dimensions_Y_space[i] = trees[i].get_number_of_leaves(); }
		return dimensions_Y_space;
	}

public:
	// we add the trees when we construct the object.
	explicit TLinks(const std::vector<TTree> &trees, size_t n_iterations)
	    : _trees(trees), _data_Y(n_iterations, _get_dimensions_Y_space(trees)) {}
	~TLinks() = default;

	bool in_any_tree(const std::string &node_id) const {
		bool result = false;
		for (const auto &tree : _trees) {
			if (tree.in_tree(node_id)) { result = true; }
		};
		return result;
	}

	bool is_leaf_in_any_tree(const std::string &node_id) const {
		bool result = false;
		for (const auto &tree : _trees) {
			if (tree.in_tree(node_id) && tree.get_node(tree.get_node_index(node_id)).isLeaf()) { result = true; }
		};
		return result;
	}

	void load_from_file(const std::string &filename);

	[[nodiscard]] const TStorageYVector &get_TStorageYVector() const { return _data_Y; }
};
