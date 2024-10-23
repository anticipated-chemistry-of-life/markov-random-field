//
// Created by Marco Visani on 26.06.23.
//

#include "TTree.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Main/TLog.h"
#include <iostream>
#include <string>

// Tree destructor implementation
TTree::~TTree() = default;

void TTree::load_from_file(const std::string &filename) {
	coretools::instances::logfile().listFlush("Reading tree from file '", filename, "' ...");
	coretools::TInputFile file(filename, coretools::FileType::Header);

	if (file.numCols() != 3) {
		UERROR("File '", filename, "' is expected to have 3 columns, but has ", file.numCols(), " !");
	}

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
			_nodes.emplace_back(parent, 0.0, -1, true);
			_node_map[parent] = _nodes.size() - 1;

			// since we just added the parent, we know that it is at the last element of the vector
			size_t parent_index = _nodes.size() - 1;
			_nodes.emplace_back(child, branch_length, parent_index, false);
			_node_map[child] = _nodes.size() - 1;
			_nodes[parent_index].addChild(_nodes.size() - 1);
		} else if (!in_tree(child) && in_tree(parent)) {
			size_t parent_index = get_node_index(parent);
			// we add the child node to the tree
			_nodes.emplace_back(child, branch_length, parent_index, false);
			_node_map[child] = _nodes.size() - 1;
			_nodes[parent_index].addChild(_nodes.size() - 1);
		} else if (in_tree(child) && !in_tree(parent)) {
			// if the child node was in the tree but not the parent
			// that means that the child was a root and is now becoming
			// a child of the new parent which will be added to the tree
			size_t child_index = get_node_index(child);
			_nodes[child_index].set_is_root(false);
			_nodes[child_index].set_branch_length_to_parent(branch_length);
			_nodes[child_index].set_parent_index(_nodes.size());
			_nodes.emplace_back(parent, 0.0, -1, true);
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
	_leafIndices.resize(_nodes.size(), -1);
	_rootIndices.resize(_nodes.size(), -1);
	for (auto it = _nodes.begin(); it != _nodes.end(); ++it) {
		if (it->isLeaf()) {
			_leafIndices[it - _nodes.begin()] = _leaves.size();
			_leaves.push_back(it - _nodes.begin());
		} else if (it->isRoot()) {
			_rootIndices[it - _nodes.begin()] = _roots.size();
			_roots.push_back(it - _nodes.begin());
		}
	}

	coretools::instances::logfile().done();
	coretools::instances::logfile().conclude("Read ", _nodes.size(), " nodes of which ", _roots.size(),
	                                         " are roots and ", _leaves.size(), " are leaves and ",
	                                         _nodes.size() - _roots.size() - _leaves.size(), " are internal nodes.");
}

TNode TTree::get_node(const std::string &Id) const {
	auto it = _node_map.find(Id);
	if (it == _node_map.end()) { UERROR("Node '", Id, "' does not exist!"); }
	return _nodes[it->second]; // Retrieve node from vector using the index
}

size_t TTree::get_node_index(const std::string &Id) const {
	auto it = _node_map.find(Id);
	if (it == _node_map.end()) { UERROR("Node '", Id, "' does not exist!"); }
	return it->second; // Return the index from the map
}
