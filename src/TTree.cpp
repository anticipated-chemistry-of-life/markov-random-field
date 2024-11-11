//
// Created by Marco Visani on 26.06.23.
//

#include "TTree.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/algorithms.h"
#include "coretools/devtools.h"
#include <cstddef>
#include <string>
#include <vector>

// Tree destructor implementation
TTree::~TTree() = default;

void TTree::_initialize_grid_branch_lengths(size_t number_of_branches) {
	// read a, b and K from command-line
	_a               = coretools::instances::parameters().get("a", coretools::Probability(0.0));
	double default_b = std::min(1.0, 1.0 / (double)number_of_branches * 10);
	_b               = coretools::instances::parameters().get("b", coretools::Probability(default_b));
	_number_of_bins  = coretools::instances::parameters().get("K", 100);
	if (_number_of_bins >= (size_t)std::numeric_limits<TypeBinBranches>::max()) {
		UERROR("More bins (", _number_of_bins, ") required than type allows (",
		       std::numeric_limits<TypeBinBranches>::max(), ")! Please decrease K or change type of bins.");
	}

	// calculate Delta
	_delta = ((double)_b - (double)_a) / (double)_number_of_bins;
}

void TTree::_bin_branch_lengths(std::vector<double> &branch_lengths) {
	_initialize_grid_branch_lengths(branch_lengths.size());
	// normalize such that they sum to one
	coretools::normalize(branch_lengths);

	std::vector<double> grid(_number_of_bins);
	for (size_t k = 0; k < _number_of_bins; ++k) { grid[k] = (_a + _delta * (k + 1)); }

	for (size_t i = 0; i < branch_lengths.size(); ++i) {
		// find bin
		auto it = std::lower_bound(grid.begin(), grid.end(), branch_lengths[i]);
		if (it == grid.end()) {
			// last bin
			_nodes[i].set_bin_branch_length_to_parent(grid.size() - 1);
		} else {
			_nodes[i].set_bin_branch_length_to_parent(std::distance(grid.begin(), it));
		}
	}
};

void TTree::load_from_file(const std::string &filename) {
	coretools::instances::logfile().listFlush("Reading tree from file '", filename, "' ...");
	coretools::TInputFile file(filename, coretools::FileType::Header);

	if (file.numCols() != 3) {
		UERROR("File '", filename, "' is expected to have 3 columns, but has ", file.numCols(), " !");
	}

	std::vector<double> branch_lengths;

	// read each line of the file
	for (; !file.empty(); file.popFront()) {
		// we read line by line the edge list from the file
		// if the child node (column 0) is not in the tree, we add the node,
		// if the parent is not in the tree, we add the node and set the parent as root
		// as long as it has no parents.
		std::string child  = std::string(file.get(0));
		std::string parent = std::string(file.get(1));
		auto branch_length = file.get<double>(2);

		if (!in_tree(child) && !in_tree(parent)) {
			// we add the parent
			_nodes.emplace_back(parent, -1, true);
			branch_lengths.push_back(0.0);
			_node_map[parent] = _nodes.size() - 1;

			// since we just added the parent, we know that it is at the last element of the vector
			size_t parent_index = _nodes.size() - 1;
			_nodes.emplace_back(child, parent_index, false);
			branch_lengths.push_back(branch_length);
			_node_map[child] = _nodes.size() - 1;
			_nodes[parent_index].addChild(_nodes.size() - 1);
		} else if (!in_tree(child) && in_tree(parent)) {
			size_t parent_index = get_node_index(parent);
			// we add the child node to the tree
			_nodes.emplace_back(child, parent_index, false);
			branch_lengths.push_back(branch_length);
			_node_map[child] = _nodes.size() - 1;
			_nodes[parent_index].addChild(_nodes.size() - 1);
		} else if (in_tree(child) && !in_tree(parent)) {
			// if the child node was in the tree but not the parent
			// that means that the child was a root and is now becoming
			// a child of the new parent which will be added to the tree
			size_t child_index = get_node_index(child);
			_nodes[child_index].set_is_root(false);
			branch_lengths[child_index] = branch_length;
			_nodes[child_index].set_parent_index(_nodes.size());
			_nodes.emplace_back(parent, -1, true);
			branch_lengths.push_back(0.0);
			_node_map[parent] = _nodes.size() - 1;
			_nodes[_nodes.size() - 1].addChild(child_index);
		} else {
			// if both nodes were already in the tree that means that the
			// child would in theroy have two parents which is not allowed
			// so we throw an error
			UERROR("Node: '", child, "' has already a parent in the tree. Adding an other parent is not allowed !");
		}
	}

	// identify roots and leaves
	this->_leafIndices.resize(_nodes.size(), -1);
	this->_rootIndices.resize(_nodes.size(), -1);
	this->_internalIndices.resize(_nodes.size(), -1);
	this->_internalIndicesWithoutRoots.resize(_nodes.size(), -1);
	for (auto it = _nodes.begin(); it != _nodes.end(); ++it) {
		if (it->isLeaf()) {
			this->_leafIndices[it - _nodes.begin()] = _leaves.size();
			this->_leaves.push_back(it - _nodes.begin());
		} else {
			this->_internalIndices[it - _nodes.begin()] = _internal_nodes.size();
			this->_internal_nodes.push_back(it - _nodes.begin());
			if (it->isRoot()) {
				_rootIndices[it - _nodes.begin()] = _roots.size();
				_roots.push_back(it - _nodes.begin());
			} else if (it->isInternalNode()) {
				_internalIndicesWithoutRoots[it - _nodes.begin()] = _internal_nodes_without_roots.size();
				_internal_nodes_without_roots.push_back(it - _nodes.begin());
			}
		}
	}

	// binned branches
	_bin_branch_lengths(branch_lengths);

	coretools::instances::logfile().done();
	coretools::instances::logfile().conclude("Read ", _nodes.size(), " nodes of which ", _roots.size(),
	                                         " are roots and ", _leaves.size(), " are leaves and ",
	                                         _nodes.size() - _roots.size() - _leaves.size(), " are internal nodes.");
}

const TNode &TTree::get_node(const std::string &Id) const {
	auto it = _node_map.find(Id);
	if (it == _node_map.end()) { UERROR("Node '", Id, "' does not exist!"); }
	return _nodes[it->second]; // Retrieve node from vector using the index
}

const TNode &TTree::get_node(size_t index) const { return _nodes[index]; }

size_t TTree::get_node_index(const std::string &Id) const {
	auto it = _node_map.find(Id);
	if (it == _node_map.end()) { UERROR("Node '", Id, "' does not exist!"); }
	return it->second; // Return the index from the map
}

void TTree::initialize_cliques(const std::vector<TTree> &all_trees) {

	// we initialize the number of leaves we have in each tree
	std::vector<size_t> num_leaves_per_tree(all_trees.size());

	// we get the number of leaves for each tree
	for (size_t i = 0; i < all_trees.size(); ++i) { num_leaves_per_tree[i] = all_trees[i].get_number_of_leaves(); }

	// the cliques of a tree are can only contain leaves in all trees except the one we are working on.
	num_leaves_per_tree[_dimension] = 1;

	// we then caclulate how many cliques we will have in total for that tree. Which is the product of the number of
	// leaves in each tree except the one we are working on (that is why we set it to 1 before).
	size_t n_cliques = coretools::containerProduct(num_leaves_per_tree);

	size_t increment = 1;
	for (size_t i = _dimension + 1; i < all_trees.size(); ++i) { increment *= all_trees[i].size(); }
	for (size_t i = 0; i < n_cliques; ++i) {
		std::vector<size_t> multi_dim_index = coretools::getSubscripts(i, num_leaves_per_tree);
		_cliques.emplace_back(multi_dim_index, _dimension, _nodes.size(), increment);
	}
}
