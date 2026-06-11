#pragma once
#include <cstddef>
#include <string>
#include <vector>

class TNode {
private:
	std::string _id;                               // unique identifier for the node
	size_t _parent_index_in_tree;                  // pointer to parent node
	std::vector<size_t> _children_indices_in_tree; // vector to child nodes indices
	bool _is_root;

public:
	// Constructor, initializes a TNode with id, _branchLengthToParent and parent node
	TNode(std::string IdString, size_t Parent_index_in_tree, bool is_root)
	    : _id(std::move(IdString)), _parent_index_in_tree(Parent_index_in_tree),
	      _is_root(is_root) {};

	// Function to add child
	void addChild_index_in_tree(size_t ChildIndex) {
		_children_indices_in_tree.push_back(ChildIndex);
	};

	// Method to check if the node is a leaf (has no _children)
	[[nodiscard]] const std::string &get_id() const { return _id; };
	[[nodiscard]] size_t parent_index_in_tree() const { return _parent_index_in_tree; };
	[[nodiscard]] size_t num_children() const { return _children_indices_in_tree.size(); };
	[[nodiscard]] const std::vector<size_t> &children_indices_in_tree() const {
		return _children_indices_in_tree;
	};
	[[nodiscard]] inline bool is_leaf() const { return _children_indices_in_tree.empty(); };
	[[nodiscard]] bool is_root() const { return _is_root; };
	[[nodiscard]] bool is_internal_node() const {
		return !_is_root && !_children_indices_in_tree.empty();
	};
	void set_is_root(bool is_root) { _is_root = is_root; }
	void set_parent_index_in_tree(size_t parent_index_in_tree) {
		_parent_index_in_tree = parent_index_in_tree;
	}

	bool operator==(const std::string &Id) const { return _id == Id; };
};
