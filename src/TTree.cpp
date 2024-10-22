//
// Created by Marco Visani on 26.06.23.
//

#include "TTree.h"
#include "coretools/Files/TInputFile.h"
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
		std::string child    = std::string(file.get(0));
		std::string parent   = std::string(file.get(1));
		double branch_length = file.get<double>(2);

		if (!in_tree(child) && !in_tree(parent)) {
			// we add the parent
			_nodes.emplace_back(parent, -1, -1, true);

			// since we just added the parent, we know that it is at the last element of the vector
			int parent_index = _nodes.size() - 1;
			_nodes.emplace_back(child, branch_length, parent_index, false);
			_nodes[parent_index].addChild(_nodes.size() - 1);
		} else if (!in_tree(child) && in_tree(parent)) {
			TNode parent_node  = get_node(parent);
			auto node_iterator = std::find(_nodes.begin(), _nodes.end(), parent_node.id());
			int parent_index   = node_iterator - _nodes.begin();
			// we add the child node to the tree
			_nodes.emplace_back(child, branch_length, parent_index, false);
			_nodes[parent_index].addChild(_nodes.size() - 1);
		} else if (in_tree(child) && !in_tree(parent)) {
			// if the child node was in the tree but not the parent
			// that means that the child was a root and is now becoming
			// a child of the new parent which will be added to the tree
			TNode child_node   = get_node(child);
			auto node_iterator = std::find(_nodes.begin(), _nodes.end(), child_node.id());
			int child_index    = node_iterator - _nodes.begin();
			_nodes[child_index].set_is_root(false);
			_nodes[child_index].set_branch_length_to_parent(branch_length);
			_nodes[child_index].set_parent_index(_nodes.size());
			_nodes.emplace_back(parent, -1, -1, true);
			_nodes[_nodes.size() - 1].addChild(child_index);
		} else {
			// if both nodes are in the tree we just need to add the child to the parent relationship
			// and set the different requirements
			TNode child_node     = get_node(child);
			TNode parent_node    = get_node(parent);
			auto child_iterator  = std::find(_nodes.begin(), _nodes.end(), child_node.id());
			auto parent_iterator = std::find(_nodes.begin(), _nodes.end(), parent_node.id());
			int child_index      = child_iterator - _nodes.begin();
			int parent_index     = parent_iterator - _nodes.begin();
			_nodes[child_index].set_is_root(false);
			_nodes[child_index].set_branch_length_to_parent(branch_length);
			_nodes[child_index].set_parent_index(parent_index);
			_nodes[parent_index].addChild(child_index);
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
	                                         " are roots and ", _leaves.size(), " are leaves.");
}

TNode TTree::get_node(std::string &Id) {
	auto node_iterator = std::find(_nodes.begin(), _nodes.end(), Id);
	if (node_iterator == _nodes.end()) { UERROR("Node '", Id, "' does not exist !"); }
	int node_index = node_iterator - _nodes.begin();
	return _nodes[node_index];
}

bool TTree::in_tree(const std::string &node) {
	auto node_iterator = std::find(_nodes.begin(), _nodes.end(), node);
	if (node_iterator == _nodes.end()) { return false; }
	return true;
}
