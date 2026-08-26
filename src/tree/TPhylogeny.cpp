#include "tree/TPhylogeny.h"

#include "coretools/Files/TInputFile.h"
#include "coretools/Main/TError.h"

#include <cstddef>
#include <queue>
#include <string>
#include <vector>

size_t TPhylogeny::index_of(const std::string &id) const {
	auto it = _index_by_id.find(id);
	if (it == _index_by_id.end()) {
		throw coretools::TUserError("Node '", id, "' does not exist in the tree !");
	}
	return it->second;
}

namespace {

/// Every node reaches a node with no parent within n_nodes steps, or it sits on a cycle.
///
/// Each node has at most one parent -- the builder rejects a second one -- so a component either
/// hangs below a root or contains exactly one cycle. Walking up is therefore enough, and it can
/// name the node it got stuck on, which "there is no root" cannot.
void reject_cycles(const std::vector<size_t> &parents, const std::vector<std::string> &ids) {
	const size_t n = parents.size();
	for (size_t start = 0; start < n; ++start) {
		size_t node = start;
		for (size_t steps = 0; steps <= n; ++steps) {
			if (parents[node] == TPhylogeny::NO_PARENT) { break; }
			node = parents[node];
			if (steps == n) {
				throw coretools::TUserError("Node '", ids[start],
				                            "' sits on a cycle: following its parents never "
				                            "reaches a root.");
			}
		}
	}
}

} // namespace

TPhylogeny build_phylogeny(const std::vector<TEdge> &edges) {
	TPhylogeny tree;

	// Children are collected per node and flattened into CSR at the end: the incremental walk below
	// does not know a node's degree until every edge has been seen.
	std::vector<std::vector<size_t>> children_of;
	std::vector<double>
	    length_by_node; // parallel to nodes; a root's entry stays 0.0 and is dropped

	auto add_node = [&](const std::string &id, size_t parent, double length) {
		tree._ids.push_back(id);
		tree._parents.push_back(parent);
		children_of.emplace_back();
		length_by_node.push_back(length);
		tree._index_by_id[id] = tree._ids.size() - 1;
		return tree._ids.size() - 1;
	};

	for (const auto &edge : edges) {
		if (edge.child == edge.parent) {
			throw coretools::TUserError("Node '", edge.child,
			                            "' can not be parent of itself ! Got ", edge.child,
			                            "for the child and ", edge.parent, " for the parent.");
		}
		if (edge.branch_length <= 0.0) {
			throw coretools::TUserError(
			    "You can't have a negative branch length or equal to 0.0 !");
		}

		const bool has_child  = tree.contains(edge.child);
		const bool has_parent = tree.contains(edge.parent);

		if (!has_child && !has_parent) {
			const size_t parent_ix = add_node(edge.parent, TPhylogeny::NO_PARENT, 0.0);
			const size_t child_ix  = add_node(edge.child, parent_ix, edge.branch_length);
			children_of[parent_ix].push_back(child_ix);
		} else if (!has_child) {
			const size_t parent_ix = tree._index_by_id.at(edge.parent);
			const size_t child_ix  = add_node(edge.child, parent_ix, edge.branch_length);
			children_of[parent_ix].push_back(child_ix);
		} else {
			// The child already exists, so it may only be adopted if nothing has claimed it yet.
			const size_t child_ix = tree._index_by_id.at(edge.child);
			if (tree._parents[child_ix] != TPhylogeny::NO_PARENT) {
				throw coretools::TUserError(
				    "Node: '", edge.child,
				    "' has already a parent in the tree. Adding an other parent is not allowed !");
			}
			const size_t parent_ix = has_parent ? tree._index_by_id.at(edge.parent)
			                                    : add_node(edge.parent, TPhylogeny::NO_PARENT, 0.0);
			tree._parents[child_ix]  = parent_ix;
			length_by_node[child_ix] = edge.branch_length;
			children_of[parent_ix].push_back(child_ix);
		}
	}

	const size_t n = tree._ids.size();

	// Flatten children into CSR so that is_leaf and children_of answer without a pointer chase.
	tree._child_offsets.resize(n + 1, 0);
	for (size_t i = 0; i < n; ++i) {
		tree._child_offsets[i + 1] = tree._child_offsets[i] + children_of[i].size();
	}
	tree._children.reserve(tree._child_offsets[n]);
	for (const auto &kids : children_of) {
		tree._children.insert(tree._children.end(), kids.begin(), kids.end());
	}

	reject_cycles(tree._parents, tree._ids);

	// Classify. The order here is the order the category lists carry, and it is load-bearing:
	// branch index is an index into _branches, which downstream code indexes branch lengths by.
	tree._leaf_index.assign(n, TPhylogeny::NOT_IN_SPACE);
	tree._internal_index.assign(n, TPhylogeny::NOT_IN_SPACE);
	tree._branch_index.assign(n, TPhylogeny::NOT_IN_SPACE);
	for (size_t i = 0; i < n; ++i) {
		if (tree.is_leaf(i)) {
			tree._leaf_index[i] = tree._leaves.size();
			tree._leaves.push_back(i);
		} else {
			tree._internal_index[i] = tree._internal_nodes.size();
			tree._internal_nodes.push_back(i);
			if (tree.is_root(i)) {
				tree._roots.push_back(i);
			} else {
				tree._internal_nodes_without_roots.push_back(i);
			}
		}
		if (!tree.is_root(i)) {
			tree._branch_index[i] = tree._branches.size();
			tree._branches.push_back(i);
		}
	}

	if (tree._roots.empty()) {
		throw coretools::TUserError("The tree has no root: every node has a parent.");
	}
	if (tree._leaves.empty()) {
		throw coretools::TUserError("The tree has no leaf: every node has a child.");
	}
	for (const size_t root : tree._roots) {
		if (tree.is_leaf(root)) {
			throw coretools::TUserError(
			    "Node '", tree._ids[root],
			    "' is both a root and a leaf: it has no parent and no children. A root must have "
			    "at least one child.");
		}
	}

	// Reachability. Cycles are already rejected and every node has at most one parent, so this
	// cannot fail today -- it is asserted because it is a post-condition callers rely on, not
	// because a known input reaches it.
	std::vector<bool> seen(n, false);
	std::queue<size_t> pending;
	for (const size_t root : tree._roots) {
		seen[root] = true;
		pending.push(root);
	}
	size_t reached = tree._roots.size();
	while (!pending.empty()) {
		const size_t node = pending.front();
		pending.pop();
		for (const size_t child : tree.children_of(node)) {
			if (!seen[child]) {
				seen[child] = true;
				++reached;
				pending.push(child);
			}
		}
	}
	if (reached != n) {
		for (size_t i = 0; i < n; ++i) {
			if (!seen[i]) {
				throw coretools::TUserError("Node '", tree._ids[i],
				                            "' is not reachable from any root.");
			}
		}
	}

	// One length per branch, in branch order: a root has no branch, so its entry is dropped here
	// rather than carried around as a zero.
	tree._branch_lengths.reserve(tree._branches.size());
	for (const size_t node : tree._branches) {
		tree._branch_lengths.push_back(length_by_node[node]);
	}

	return tree;
}

TPhylogeny read_phylogeny(const std::string &filename) {
	coretools::TInputFile file(filename, coretools::FileType::Header);
	if (file.numCols() != 3) {
		throw coretools::TUserError("File '", filename, "' is expected to have 3 columns, but has ",
		                            file.numCols(), " !");
	}

	std::vector<TEdge> edges;
	for (; !file.empty(); file.popFront()) {
		edges.push_back({std::string(file.get(0)), std::string(file.get(1)), file.get<double>(2)});
	}

	return build_phylogeny(edges);
}
