#include "../TTree.h"

void TTree::_load_from_file(const std::string &filename, const std::string &tree_name) {
	coretools::instances::logfile().listFlush("Reading tree from file '", filename, "' ...");
	coretools::TInputFile file(filename, coretools::FileType::Header);
	this->_tree_name = tree_name;
	if (file.numCols() != 3) {
		throw coretools::TUserError("File '", filename, "' is expected to have 3 columns, but has ", file.numCols(),
		                            " !");
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
		if (child == parent) {
			throw coretools::TUserError("Node '", child, "' can not be parent of itself ! Got ", child,
			                            "for the child and ", parent, " for the parent.");
		}
		auto branch_length = file.get<double>(2);
		if (branch_length <= 0.0) {
			throw coretools::TUserError("You can't have a negative branch length or equal to 0.0 !");
		}

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
				throw coretools::TUserError(
				    "Node: '", child, "' has already a parent in the tree. Adding an other parent is not allowed !");
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
				throw coretools::TUserError(
				    "Node: '", child, "' has already a parent in the tree. Adding an other parent is not allowed !");

			} else {
				// if the child was a root and the parent was already in the tree
				// we set the parent as the parent of the child
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
	this->_leaves_and_internal_nodes_without_roots_indices.resize(_nodes.size(), -1);
	for (auto it = _nodes.begin(); it != _nodes.end(); ++it) {
		if (it->isLeaf()) {
			this->_leafIndices[it - _nodes.begin()] = _leaves.size();
			this->_leaves.push_back(it - _nodes.begin());

			this->_leaves_and_internal_nodes_without_roots_indices[it - _nodes.begin()] =
			    _leaves_and_internal_nodes_without_roots.size();
			this->_leaves_and_internal_nodes_without_roots.push_back(it - _nodes.begin());
		} else {
			this->_internalIndices[it - _nodes.begin()] = _internal_nodes.size();
			this->_internal_nodes.push_back(it - _nodes.begin());
			if (it->isRoot()) {
				_rootIndices[it - _nodes.begin()] = _roots.size();
				_roots.push_back(it - _nodes.begin());
			} else if (it->isInternalNode()) {
				_internalIndicesWithoutRoots[it - _nodes.begin()] = _internal_nodes_without_roots.size();
				_internal_nodes_without_roots.push_back(it - _nodes.begin());

				this->_leaves_and_internal_nodes_without_roots_indices[it - _nodes.begin()] =
				    _leaves_and_internal_nodes_without_roots.size();
				this->_leaves_and_internal_nodes_without_roots.push_back(it - _nodes.begin());
			}
		}
	}

	// binned branches
	_bin_branch_lengths_from_tree(branch_lengths);

	coretools::instances::logfile().done();
	coretools::instances::logfile().conclude("Read ", _nodes.size(), " nodes of which ", _roots.size(),
	                                         " are roots and ", _leaves.size(), " are leaves and ",
	                                         _internal_nodes_without_roots.size(), " are internal nodes.");
}
