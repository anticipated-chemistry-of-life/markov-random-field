//
// Created by Marco Visani on 26.06.23.
//

#ifndef METABOLITE_INFERENCE_TREE_H
#define METABOLITE_INFERENCE_TREE_H

#include "Types.h"
#include "coretools/Types/probability.h"
#include <cstddef>
#include <string>
#include <vector>
class TClique;

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
	[[nodiscard]] bool isInternalNode() const { return !_is_root && !_children.empty(); };
	[[nodiscard]] TypeBinBranches get_branch_length_bin() const { return _binned_branch_length; };
	void set_is_root(bool is_root) { _is_root = is_root; }
	void set_bin_branch_length_to_parent(TypeBinBranches branch_length) { _binned_branch_length = branch_length; }
	void set_parent_index(size_t parent_index) { _parentIndex = parent_index; }

	bool operator==(const std::string &Id) const { return _id == Id; };
};

class TTree {
public:
	explicit TTree(size_t dimension);
	TTree()  = default;
	// Destructor, deletes all nodes
	~TTree() = default;

	size_t size() const { return _nodes.size(); };

	/** Get node by its id
	 * @param Id: the id of the node
	 * @return a reference to the node with the given id
	 */
	const TNode &get_node(const std::string &Id) const;

	/** Get node by its index
	 * @param Id: the index of the node wihtin the tree
	 * @return a reference to the Node with the given id
	 */
	const TNode &get_node(size_t index) const;

	/** Get the index of a node by its id
	 * @param Id: the id of the node
	 * @return the index of the node with the given id
	 */
	size_t get_node_index(const std::string &Id) const;

	/** Load a tree from a file
	 * @param filename: the name of the file to load the tree from. This should contain three columns: child, parent,
	 * branch_length.
	 * @return the loaded tree
	 */
	void load_from_file(const std::string &filename);

	/** Gives the number of roots within the tree
	 * @return the number of roots
	 */
	[[nodiscard]] size_t count_roots() const { return _roots.size(); }

	/** Method to get all the leaves of the tree.
	 * @return Returns a vector of length equal to the number of leaves in the tree. Each element of the vector is the
	 * index of the leaf node within the tree.
	 */
	[[nodiscard]] const std::vector<size_t> &get_leaf_nodes() const { return _leaves; }

	/** @return the number of leaves in the tree
	 */
	size_t get_number_of_leaves() const { return _leaves.size(); }
	size_t get_number_of_nodes() const { return _nodes.size(); }
	size_t get_number_of_internal_nodes() const { return _internal_nodes.size(); }

	/** @param node_index: the index of the node within the tree
	 * @return The index of the node within the leaves vector (which is smaller than the total number of nodes in the
	 * tree). If the node is not a leaf, the function will return -1.
	 */
	size_t get_index_within_leaves(size_t node_index) const { return _leafIndices[node_index]; }

	/** @param node_index: the index of the node within the tree
	 * @return The index of the node within the internal nodes vector (which is smaller than the total number of nodes
	 * in the tree). If the node is not an internal node, the function will return -1.
	 */
	size_t get_index_within_internal_nodes(size_t node_index) const { return _internalIndices[node_index]; }

	/** @return The root nodes of the tree
	 */
	[[nodiscard]] const std::vector<size_t> &get_root_nodes() const { return _roots; }
	const std::vector<size_t> &get_internal_nodes() const { return _internal_nodes; }
	const std::vector<size_t> &get_internal_nodes_without_roots() const { return _internal_nodes_without_roots; }

	/** Checks whether a node is in the tree
	 * @param node_id: the id of the node
	 * @return true if the node is in the tree, false otherwise
	 */
	bool in_tree(const std::string &node_id) const { return _node_map.find(node_id) != _node_map.end(); };

	std::vector<TypeBinBranches> get_all_binned_branch_lengths() const {
		std::vector<TypeBinBranches> branch_lengths;
		branch_lengths.resize(_nodes.size());
		for (size_t i = 0; i < _nodes.size(); i++) { branch_lengths[i] = _nodes[i].get_branch_length_bin(); }
		return branch_lengths;
	}

	void initialize_cliques(const std::vector<TTree> &trees);

	coretools::Probability get_a() const { return _a; }
	coretools::Probability get_b() const { return _b; }
	double get_delta() const { return _delta; }
	size_t get_number_of_bins() const { return _number_of_bins; }
	const std::vector<TClique> &get_cliques() const { return _cliques; }

private:
	std::vector<TNode> _nodes;                         // a map to store nodes with their ids
	std::unordered_map<std::string, size_t> _node_map; // for fast access to nodes
	std::vector<size_t> _leaves;
	std::vector<size_t> _roots;
	std::vector<size_t> _internal_nodes;
	std::vector<size_t> _internal_nodes_without_roots;
	std::vector<size_t> _leafIndices;
	std::vector<size_t> _rootIndices;
	std::vector<size_t> _internalIndices;
	std::vector<size_t> _internalIndicesWithoutRoots;

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
};
#endif // METABOLITE_INFERENCE_TREE_H
