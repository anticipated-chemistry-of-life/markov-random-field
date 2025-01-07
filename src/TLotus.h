#pragma once

#include "TStorageYVector.h"
#include "TTree.h"
#include <cstddef>
#include <string>
#include <vector>

// For this class, we need to enforce that the first tree will always be the one of the species and
// the second tree will always be the one of the molecules. The rest of the trees, we won't care.
class TLotus {
private:
	// trees should be a const ref because we don't want to change the trees and don't want to copy them
	const std::vector<TTree> &_trees;
	TStorageYVector _data_L;
	std::unordered_map<std::string, size_t> _species_counter;
	std::unordered_map<std::string, size_t> _molecules_counter;

	// private functions
	static std::vector<size_t> _get_dimensions_Lotus_space(const std::vector<TTree> &trees) {
		std::vector<size_t> dimensions_Y_space(trees.size());
		for (size_t i = 0; i < trees.size(); ++i) { dimensions_Y_space[i] = trees[i].get_number_of_leaves(); }

		// return only the first two dimensions i.e. the species and the molecules
		return {dimensions_Y_space[0], dimensions_Y_space[1]};
	}

	size_t _get_tree_index_of_node(const std::string &node_id) const {
		for (size_t i = 0; i < _trees.size(); ++i) {
			const auto &tree = _trees[i];
			if (tree.in_tree(node_id)) {
				if (tree.get_node(tree.get_node_index(node_id)).isLeaf()) {
					return i;
				} else {
					UERROR("Node '", node_id,
					       "' is not a leaf ! So far, we have defined our model to only accept leaves.");
				}
			}
		};
		UERROR("Node '", node_id, "' doesn't exist in any of the provided trees !");
	}

public:
	// we add the trees when we construct the object.
	explicit TLotus(const std::vector<TTree> &trees) : _trees(trees), _data_L(0, _get_dimensions_Lotus_space(trees)) {}

	// default destructor
	~TLotus() = default;

	void load_from_file(const std::string &filename);

	[[nodiscard]] const TStorageYVector &get_TStorageYVector() const { return _data_L; }

	// To construct the _data_X_of_Lotus we will need to have as input a vector of TStorageYVector but this time with as
	// many dimensions as there are trees.
};
