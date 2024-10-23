//
// Created by Marco Visani on 26.06.23.
//

#ifndef METABOLITE_INFERENCE_TREE_H
#define METABOLITE_INFERENCE_TREE_H
#include "TBranchLengths.h"
#include <cstddef>
#include <string>
#include <vector>

class TNode {
private:
	std::string _id;              // unique identifier for the node
	size_t _parentIndex;          // pointer to parent node
	double _branchLengthToParent; // branch length to parent node
	// size_t _discretised_branch_length{}; // discretised branch length
	std::vector<size_t> _children; // vector to child nodes indices
	bool _is_root;

public:
	// Constructor, initializes a TNode with id, _branchLengthToParent and parent node
	TNode(std::string IdString, double BranchLengthToParent, size_t Parent, bool is_root)
	    : _id(std::move(IdString)), _parentIndex(Parent), _branchLengthToParent(BranchLengthToParent),
	      _is_root(is_root) {};

	// Function to add child
	void addChild(size_t ChildIndex) { _children.push_back(ChildIndex); };

	// Method to check if the node is a leaf (has no _children)
	[[nodiscard]] const std::string &id() const { return _id; };
	[[nodiscard]] size_t parentIndex() const { return _parentIndex; };
	[[nodiscard]] size_t numChildren() const { return _children.size(); };
	[[nodiscard]] const std::vector<size_t> &children() const { return _children; };
	[[nodiscard]] bool isLeaf() const { return _children.empty(); };
	[[nodiscard]] bool isRoot() const { return _is_root; };
	[[nodiscard]] double branchLengthToParent() const { return _branchLengthToParent; };
	void set_is_root(bool is_root) { _is_root = is_root; }
	void set_branch_length_to_parent(double branch_length) { _branchLengthToParent = branch_length; }
	void set_parent_index(size_t parent_index) { _parentIndex = parent_index; }

	bool operator==(const std::string &Id) const { return _id == Id; };
};

class TTree {
private:
	std::vector<TNode> _nodes;                         // a map to store nodes with their ids
	std::unordered_map<std::string, size_t> _node_map; // for fast access to nodes
	std::vector<size_t> _leaves;
	std::vector<size_t> _roots;
	std::vector<size_t> _leafIndices;
	std::vector<size_t> _rootIndices;

public:
	// Destructor, deletes all nodes
	~TTree();

	// Get node by its id
	TNode get_node(const std::string &Id) const;

	// Get the node index by its id
	size_t get_node_index(const std::string &Id) const;

	// Method to load the tree from a file
	void load_from_file(const std::string &filename);

	// Find the number of roots in the tree
	[[nodiscard]] size_t count_roots() const { return _roots.size(); }

	// Method to get all the leaves of the tree
	[[nodiscard]] const std::vector<size_t> &get_leaf_nodes() const { return _leaves; }

	// Method to get all the root nodes
	[[nodiscard]] const std::vector<size_t> &get_root_nodes() const { return _roots; }

	bool in_tree(const std::string &node_id) const { return _node_map.find(node_id) != _node_map.end(); };

	std::vector<double> get_all_branch_lengths() const {
		std::vector<double> branch_lengths;
		branch_lengths.resize(_nodes.size());
		for (size_t i = 0; i < _nodes.size(); i++) { branch_lengths[i] = _nodes[i].branchLengthToParent(); }
		return branch_lengths;
	}
};
#endif // METABOLITE_INFERENCE_TREE_H
