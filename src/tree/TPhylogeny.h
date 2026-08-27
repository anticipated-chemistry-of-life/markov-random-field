#ifndef ACOL_TPHYLOGENY_H
#define ACOL_TPHYLOGENY_H

#include <cstddef>
#include <limits>
#include <span>
#include <string>
#include <unordered_map>
#include <vector>

/// One edge exactly as a tree file states it: a child, its parent, and the length of the branch
/// between them. This is what build_phylogeny validates, and it is the reason a phylogeny can be
/// built in a test without writing a file.
struct TEdge {
	std::string child;
	std::string parent;
	double branch_length;
};

class TPhylogeny;
struct TPhylogenyBuilder;

/// Turn an edge list into a validated phylogeny. The only way to build one.
///
/// Throws coretools::TUserError if the edges do not describe a forest: a node claimed as a child
/// twice, a node that is its own parent, a non-positive branch length, a cycle, a childless root,
/// no root at all, no leaf at all, or a node unreachable from any root.
[[nodiscard]] TPhylogeny build_phylogeny(const std::vector<TEdge> &edges);

/// Parse a tree file into edges and build a phylogeny from them. File shape is validated here,
/// structure in build_phylogeny.
[[nodiscard]] TPhylogeny read_phylogeny(const std::string &filename);

/// A tree's topology and nothing else: which nodes exist, how they are related, and how long each
/// branch is.
///
/// It carries no parameters, no cliques and no field, and in particular no stattools and no
/// command-line options -- which is exactly what lets it be built inside a unit test. TTree adds
/// the parameters and the MCMC moves over them.
///
/// Nodes are addressed by index, and stored in canonical order: all leaves, then all internal
/// non-root nodes in post-order, then all roots, with file order preserved inside the leaf block
/// and the root block. So a node's index alone says what kind of node it is, and a child's index
/// is always below its parent's. See ADR-0004.
class TPhylogeny {
public:
	/// The parent of a root. Never a valid node index.
	static constexpr size_t NO_PARENT    = std::numeric_limits<size_t>::max();
	/// Returned when a node has no index in the space asked for -- a leaf has no internal-node
	/// index, a root has no branch index.
	static constexpr size_t NOT_IN_SPACE = std::numeric_limits<size_t>::max();

private:
	std::vector<std::string> _ids;
	std::vector<size_t> _parents;        // NO_PARENT for a root
	std::vector<size_t> _child_offsets;  // CSR, size n_nodes + 1
	std::vector<size_t> _children;       // CSR, flat
	std::vector<double> _branch_lengths; // one per branch, in branch order
	std::unordered_map<std::string, size_t> _index_by_id;

	// Category lists, each in ascending node order.
	std::vector<size_t> _leaves;
	std::vector<size_t> _roots;
	std::vector<size_t> _internal_nodes; // roots included
	std::vector<size_t> _internal_nodes_without_roots;
	std::vector<size_t> _branches; // every non-root node: the node its branch hangs below

	// Maps back into those lists, each of size n_nodes, NOT_IN_SPACE where inapplicable.
	std::vector<size_t> _leaf_index;
	std::vector<size_t> _internal_index;
	std::vector<size_t> _branch_index;

	TPhylogeny() = default;
	/// The only thing that may fill one in. Both build_phylogeny and read_phylogeny go through it,
	/// which is what lets read_phylogeny stream a file rather than materialise its edges first.
	friend struct TPhylogenyBuilder;

	/// Out of line so that this header stays free of the error machinery.
	[[noreturn]] void _throw_no_branch(size_t node) const;

public:
	[[nodiscard]] size_t n_nodes() const { return _ids.size(); }
	[[nodiscard]] size_t n_leaves() const { return _leaves.size(); }
	[[nodiscard]] size_t n_roots() const { return _roots.size(); }
	[[nodiscard]] size_t n_branches() const { return _branches.size(); }
	[[nodiscard]] size_t n_internal_nodes() const { return _internal_nodes.size(); }

	[[nodiscard]] const std::string &id_of(size_t node) const { return _ids[node]; }
	[[nodiscard]] bool contains(const std::string &id) const { return _index_by_id.count(id) > 0; }
	/// Throws coretools::TUserError if no node carries this id.
	[[nodiscard]] size_t index_of(const std::string &id) const;

	[[nodiscard]] size_t parent_of(size_t node) const { return _parents[node]; }
	[[nodiscard]] std::span<const size_t> children_of(size_t node) const {
		return {_children.data() + _child_offsets[node],
		        _child_offsets[node + 1] - _child_offsets[node]};
	}

	[[nodiscard]] bool is_leaf(size_t node) const {
		return _child_offsets[node + 1] == _child_offsets[node];
	}
	[[nodiscard]] bool is_root(size_t node) const { return _parents[node] == NO_PARENT; }

	[[nodiscard]] size_t leaf_index(size_t node) const { return _leaf_index[node]; }
	[[nodiscard]] size_t internal_index(size_t node) const { return _internal_index[node]; }
	[[nodiscard]] size_t branch_index(size_t node) const { return _branch_index[node]; }

	/// The branch below `node`, exactly as the file stated it. Throws a dev error for a root,
	/// which has no branch; ask is_root first.
	[[nodiscard]] double branch_length_of(size_t node) const {
		if (_branch_index[node] == NOT_IN_SPACE) { _throw_no_branch(node); }
		return _branch_lengths[_branch_index[node]];
	}
	[[nodiscard]] const std::vector<double> &branch_lengths() const { return _branch_lengths; }

	[[nodiscard]] const std::vector<size_t> &leaves() const { return _leaves; }
	[[nodiscard]] const std::vector<size_t> &roots() const { return _roots; }
	[[nodiscard]] const std::vector<size_t> &internal_nodes() const { return _internal_nodes; }
	[[nodiscard]] const std::vector<size_t> &internal_nodes_without_roots() const {
		return _internal_nodes_without_roots;
	}
	[[nodiscard]] const std::vector<size_t> &branches() const { return _branches; }
};

#endif // ACOL_TPHYLOGENY_H
