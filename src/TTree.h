//
// Created by Marco Visani on 26.06.23.
//

#ifndef METABOLITE_INFERENCE_TREE_H
#define METABOLITE_INFERENCE_TREE_H
#include "TClique.h"
#include "Types.h"
#include "coretools/Types/probability.h"
#include <cstddef>
#include <string>
#include <vector>

class TNode {
private:
	std::string _id;                         // unique identifier for the node
	size_t _parentIndex;                     // pointer to parent node
	TypeBinBranches _binned_branch_length{}; // discretised branch length
	std::vector<size_t> _children;           // vector to child nodes indices
	bool _is_root;

public:
	// Constructor, initializes a TNode with id, _branchLengthToParent and parent node
	TNode(std::string IdString, size_t Parent, bool is_root)
	    : _id(std::move(IdString)), _parentIndex(Parent), _is_root(is_root) {};

	// Function to add child
	void addChild(size_t ChildIndex) { _children.push_back(ChildIndex); };

	// Method to check if the node is a leaf (has no _children)
	[[nodiscard]] const std::string &id() const { return _id; };
	[[nodiscard]] size_t parentIndex() const { return _parentIndex; };
	[[nodiscard]] size_t numChildren() const { return _children.size(); };
	[[nodiscard]] const std::vector<size_t> &children() const { return _children; };
	[[nodiscard]] bool isLeaf() const { return _children.empty(); };
	[[nodiscard]] bool isRoot() const { return _is_root; };
	[[nodiscard]] TypeBinBranches get_branch_length_bin() const { return _binned_branch_length; };
	void set_is_root(bool is_root) { _is_root = is_root; }
	void set_bin_branch_length_to_parent(TypeBinBranches branch_length) { _binned_branch_length = branch_length; }
	void set_parent_index(size_t parent_index) { _parentIndex = parent_index; }

	bool operator==(const std::string &Id) const { return _id == Id; };
};

class TTree {
private:
	std::vector<TNode> _nodes;                         // a map to store nodes with their ids
	std::unordered_map<std::string, size_t> _node_map; // for fast access to nodes
	std::vector<size_t> _leaves;
	std::vector<size_t> _roots;
	std::vector<size_t> _internal_nodes;
	std::vector<size_t> _leafIndices;
	std::vector<size_t> _rootIndices;
	std::vector<size_t> _internalIndices;

	// For binning branch lengths
	coretools::Probability _a;
	coretools::Probability _b;
	double _delta          = 0.0;
	size_t _number_of_bins = 0;

	// cliques
	std::vector<TClique> _cliques;

	// dimension of the tree
	size_t _dimension;

	void _bin_branch_lengths(std::vector<double> &branch_lengths);
	void _initialize_grid_branch_lengths(size_t number_of_branches);

public:
	explicit TTree(size_t dimension) { _dimension = dimension; };
	// Destructor, deletes all nodes
	~TTree();

	// Get node by its id
	const TNode &get_node(const std::string &Id) const;
	const TNode &get_node(size_t index) const;

	// Get the node index by its id
	size_t get_node_index(const std::string &Id) const;

	// Method to load the tree from a file
	void load_from_file(const std::string &filename);

	// Find the number of roots in the tree
	[[nodiscard]] size_t count_roots() const { return _roots.size(); }

	// Method to get all the leaves of the tree
	[[nodiscard]] const std::vector<size_t> &get_leaf_nodes() const { return _leaves; }
	size_t get_number_of_leaves() const { return _leaves.size(); }
	size_t get_index_within_leaves(size_t node_index) const { return _leafIndices[node_index]; }
	size_t get_index_within_internal_nodes(size_t node_index) const { return _internalIndices[node_index]; }
	// Method to get all the root nodes
	[[nodiscard]] const std::vector<size_t> &get_root_nodes() const { return _roots; }

	bool in_tree(const std::string &node_id) const { return _node_map.find(node_id) != _node_map.end(); };

	std::vector<TypeBinBranches> get_all_binned_branch_lengths() const {
		std::vector<TypeBinBranches> branch_lengths;
		branch_lengths.resize(_nodes.size());
		for (size_t i = 0; i < _nodes.size(); i++) { branch_lengths[i] = _nodes[i].get_branch_length_bin(); }
		return branch_lengths;
	}

	void initialize_cliques(const std::vector<TTree> &trees);
};
#endif // METABOLITE_INFERENCE_TREE_H
