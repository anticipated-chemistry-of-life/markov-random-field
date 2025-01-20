//
// Created by Marco Visani on 26.06.23.
//

#include "TTree.h"
#include "TClique.h"
#include "TStorageZ.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/algorithms.h"

#include "omp.h"
#include <cstddef>
#include <string>
#include <vector>

TTree::TTree(size_t dimension, const std::string &filename, const std::string &tree_name) {
	_dimension = dimension;
	_load_from_file(filename, tree_name);
}

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
	for (size_t k = 0; k < _number_of_bins; ++k) { grid[k] = (_a + _delta * (k + 1.0)); }

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

void TTree::_load_from_file(const std::string &filename, const std::string &tree_name) {
	coretools::instances::logfile().listFlush("Reading tree from file '", filename, "' ...");
	coretools::TInputFile file(filename, coretools::FileType::Header);
	this->_tree_name = tree_name;
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
			_nodes[parent_index].addChild_index_in_tree(_nodes.size() - 1);
		} else if (!in_tree(child) && in_tree(parent)) {
			size_t parent_index = get_node_index(parent);
			// we add the child node to the tree
			_nodes.emplace_back(child, parent_index, false);
			branch_lengths.push_back(branch_length);
			_node_map[child] = _nodes.size() - 1;
			_nodes[parent_index].addChild_index_in_tree(_nodes.size() - 1);
		} else if (in_tree(child) && !in_tree(parent)) {
			// if the child node was in the tree but not the parent
			// that means that the child was a root and is now becoming
			// a child of the new parent which will be added to the tree
			size_t child_index = get_node_index(child);

			if (!_nodes[child_index].isRoot()) {
				// if the child was not a root and the parent was not in the tree
				// we throw an error because the child already has a parent
				UERROR("Node: '", child, "' has already a parent in the tree. Adding an other parent is not allowed !");
			}
			_nodes[child_index].set_is_root(false);
			branch_lengths[child_index] = branch_length;
			_nodes[child_index].set_parent_index_in_tree(_nodes.size());
			_nodes.emplace_back(parent, -1, true);
			branch_lengths.push_back(0.0);
			_node_map[parent] = _nodes.size() - 1;
			_nodes[_nodes.size() - 1].addChild_index_in_tree(child_index);
		} else {
			size_t child_index  = get_node_index(child);
			size_t parent_index = get_node_index(parent);
			auto node           = _nodes[child_index];
			auto node_parent    = _nodes[parent_index];
			if (!node.isRoot()) {
				// if the child was not a root and the parent was already in the tree
				// we throw an error because the child already has a parent
				UERROR("Node: '", child, "' has already a parent in the tree. Adding an other parent is not allowed !");

			} else {
				// if the child was a root and the parent was already in the tree
				// we set the parent as the parent of the child
				size_t parent_index = get_node_index(parent);
				_nodes[child_index].set_is_root(false);
				branch_lengths[child_index] = branch_length;
				_nodes[child_index].set_parent_index_in_tree(parent_index);
				_nodes[parent_index].addChild_index_in_tree(child_index);
			}
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
	                                         _internal_nodes_without_roots.size(), " are internal nodes.");
}

const TNode &TTree::get_node(const std::string &Id) const {
	auto it = _node_map.find(Id);
	if (it == _node_map.end()) { UERROR("Node '", Id, "' does not exist!"); }
	return _nodes[it->second]; // Retrieve node from vector using the index
}

const TNode &TTree::get_node(size_t index) const { return _nodes[index]; }
bool TTree::isLeaf(size_t index) const { return _nodes[index].isLeaf(); }

size_t TTree::get_node_index(const std::string &Id) const {
	auto it = _node_map.find(Id);
	if (it == _node_map.end()) { UERROR("Node '", Id, "' does not exist in the tree !"); }
	return it->second; // Return the index from the map
}

void TTree::initialize_cliques_and_Z(const std::vector<TTree> &all_trees) {

	// we initialize the number of leaves we have in each tree
	std::vector<size_t> num_leaves_per_tree(all_trees.size());
	for (size_t i = 0; i < all_trees.size(); ++i) { num_leaves_per_tree[i] = all_trees[i].get_number_of_leaves(); }

	_initialize_Z(num_leaves_per_tree);
	_initialize_cliques(num_leaves_per_tree, all_trees);
}

void TTree::_initialize_Z(std::vector<size_t> num_leaves_per_tree) {
	num_leaves_per_tree[_dimension] = this->get_number_of_internal_nodes();

	_Z.initialize_dimensions(num_leaves_per_tree);
}

void TTree::_initialize_cliques(const std::vector<size_t> &num_leaves_per_tree, const std::vector<TTree> &all_trees) {
	// clique of a tree: runs along that dimension
	// the cliques of a tree are can only contain leaves in all trees except the one we are working on.
	_dimension_cliques             = num_leaves_per_tree;
	_dimension_cliques[_dimension] = 1;

	// we then caclulate how many cliques we will have in total for that tree. Which is the product of the number of
	// leaves in each tree except the one we are working on (that is why we set it to 1 before).
	const size_t n_cliques = coretools::containerProduct(_dimension_cliques);

	// calculate increment: product of the number of leaves of all subsequent dimensions
	size_t increment = 1;
	for (size_t i = _dimension + 1; i < all_trees.size(); ++i) { increment *= all_trees[i].get_number_of_leaves(); }

	// initialize cliques
	for (size_t i = 0; i < n_cliques; ++i) {
		// get start index of each clique in leaves space
		std::vector<size_t> start_index_in_leaves_space = coretools::getSubscripts(i, _dimension_cliques);
		_cliques.emplace_back(start_index_in_leaves_space, _dimension, _nodes.size(), increment);
	}
}

void TTree::update_Z(const TStorageYVector &Y) {
	std::vector<std::vector<TStorageZ>> indices_to_insert(this->_cliques.size());

#pragma omp parallel for num_threads(NUMBER_OF_THREADS) schedule(static)
	for (size_t i = 0; i < _cliques.size(); ++i) { indices_to_insert[i] = _cliques[i].update_Z(Y, _Z, *this); }
	_Z.insert_in_Z(indices_to_insert);
}

const TStorageZVector &TTree::get_Z() const { return _Z; };
std::vector<TClique> &TTree::get_cliques() { return this->_cliques; }
const TClique &TTree::get_clique(std::vector<size_t> index_in_leaves_space) const {
	index_in_leaves_space[_dimension] = 0; // set to start index
	const size_t ix_clique            = coretools::getLinearIndex(index_in_leaves_space, _dimension_cliques);
	return _cliques[ix_clique];
}

TClique &TTree::get_clique(std::vector<size_t> index_in_leaves_space) {
	index_in_leaves_space[_dimension] = 0; // set to start index
	const size_t ix_clique            = coretools::getLinearIndex(index_in_leaves_space, _dimension_cliques);
	return _cliques[ix_clique];
}
