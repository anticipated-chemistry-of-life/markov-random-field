#ifndef ACOL_TPHYLOGENY_H
#define ACOL_TPHYLOGENY_H

#include <cstddef>
#include <limits>
#include <ranges>
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
///
/// That is why every question about which category a node belongs to, and about its index within
/// that category, is arithmetic here rather than a lookup: the whole of it derives from three
/// numbers -- the node count, the leaf count and the root count. The arithmetic stays behind named
/// methods rather than being written out at call sites, so that the layout has exactly one home.
class TPhylogeny {
public:
	/// The parent of a root. Never a valid node index.
	static constexpr size_t NO_PARENT = std::numeric_limits<size_t>::max();

private:
	std::vector<std::string> _ids;
	std::vector<size_t> _parents;        // NO_PARENT for a root
	std::vector<size_t> _child_offsets;  // CSR, size n_nodes + 1
	std::vector<size_t> _children;       // CSR, flat
	std::vector<double> _branch_lengths; // one per non-root node, indexed by that node
	std::unordered_map<std::string, size_t> _index_by_id;

	// The layout, in full: with the node count these decide every category and every index space.
	size_t _n_leaves = 0;
	size_t _n_roots  = 0;

	TPhylogeny() = default;
	/// The only thing that may fill one in. Both build_phylogeny and read_phylogeny go through it,
	/// which is what lets read_phylogeny stream a file rather than materialise its edges first.
	friend struct TPhylogenyBuilder;

	/// Out of line so that this header stays free of the error machinery.
	[[noreturn]] void _throw_no_branch(size_t node) const;

public:
	[[nodiscard]] size_t n_nodes() const { return _ids.size(); }
	[[nodiscard]] size_t n_leaves() const { return _n_leaves; }
	[[nodiscard]] size_t n_roots() const { return _n_roots; }
	/// Every node except a root has exactly one branch, and the non-root nodes are the first two
	/// blocks -- so this is also the index one past the last branch, and the first root.
	[[nodiscard]] size_t n_branches() const { return n_nodes() - _n_roots; }
	/// Roots included, as ever: an internal node is any node with a child.
	[[nodiscard]] size_t n_internal_nodes() const { return n_nodes() - _n_leaves; }

	[[nodiscard]] const std::string &id_of(size_t node) const { return _ids[node]; }
	[[nodiscard]] bool contains(const std::string &id) const { return _index_by_id.count(id) > 0; }
	/// Throws coretools::TUserError if no node carries this id.
	[[nodiscard]] size_t index_of(const std::string &id) const;

	[[nodiscard]] size_t parent_of(size_t node) const { return _parents[node]; }
	[[nodiscard]] std::span<const size_t> children_of(size_t node) const {
		return {_children.data() + _child_offsets[node],
		        _child_offsets[node + 1] - _child_offsets[node]};
	}

	/// Leaves are the first block, so this is a comparison and not a look at the node's children.
	/// It is the most frequently asked question in the sweep.
	[[nodiscard]] bool is_leaf(size_t node) const { return node < _n_leaves; }
	/// Roots are the last block.
	[[nodiscard]] bool is_root(size_t node) const { return node >= n_branches(); }

	/// A leaf's index in leaf space is its node index. Meaningless for a non-leaf.
	[[nodiscard]] size_t leaf_index(size_t node) const { return node; }
	/// An internal node's index in internal-node space. Meaningless for a leaf.
	[[nodiscard]] size_t internal_index(size_t node) const { return node - _n_leaves; }
	/// A branch is identified by the node it hangs below, and that node's index is its index in
	/// branch space. Meaningless for a root, which has no branch.
	[[nodiscard]] size_t branch_index(size_t node) const { return node; }

	/// The branch below `node`, exactly as the file stated it. Throws a dev error for a root,
	/// which has no branch; ask is_root first.
	[[nodiscard]] double branch_length_of(size_t node) const {
		if (is_root(node)) { _throw_no_branch(node); }
		return _branch_lengths[node];
	}
	[[nodiscard]] const std::vector<double> &branch_lengths() const { return _branch_lengths; }

	// Each category is a contiguous block, so a category is a range of node indices rather than a
	// list of them.
	[[nodiscard]] auto leaves() const { return std::views::iota(size_t{0}, _n_leaves); }
	[[nodiscard]] auto roots() const { return std::views::iota(n_branches(), n_nodes()); }
	[[nodiscard]] auto internal_nodes() const { return std::views::iota(_n_leaves, n_nodes()); }
	[[nodiscard]] auto internal_nodes_without_roots() const {
		return std::views::iota(_n_leaves, n_branches());
	}
	[[nodiscard]] auto branches() const { return std::views::iota(size_t{0}, n_branches()); }
};

#endif // ACOL_TPHYLOGENY_H
