#include "tree/TPhylogeny.h"

#include "coretools/Files/TInputFile.h"
#include "coretools/Main/TError.h"

#include <cstddef>
#include <queue>
#include <string>
#include <utility>
#include <vector>

size_t TPhylogeny::index_of(const std::string &id) const {
	auto it = _index_by_id.find(id);
	if (it == _index_by_id.end()) {
		throw coretools::TUserError("Node '", id, "' does not exist in the tree !");
	}
	return it->second;
}

void TPhylogeny::_throw_no_branch(size_t node) const {
	throw coretools::TDevError("Node '", _ids[node],
	                           "' is a root: it has no branch, and therefore no branch length.");
}

/// Everything that can be done while edges arrive one at a time, split from the whole-tree
/// post-processing that needs them all. That split is what lets read_phylogeny stream a file
/// instead of materialising its edges: a taxonomy of a few million rows would otherwise sit in
/// memory in full, in string form, next to the tree it is being turned into.
struct TPhylogenyBuilder {
	TPhylogeny tree;

	// Children are collected per node and flattened into CSR by finish(): an incremental walk does
	// not know a node's degree until every edge has been seen.
	std::vector<std::vector<size_t>> children_of;
	// Parallel to nodes; a root's entry stays 0.0 and is dropped by finish().
	std::vector<double> length_by_node;

	size_t add_node(const std::string &id, size_t parent, double length) {
		tree._ids.push_back(id);
		tree._parents.push_back(parent);
		children_of.emplace_back();
		length_by_node.push_back(length);
		tree._index_by_id[id] = tree._ids.size() - 1;
		return tree._ids.size() - 1;
	}

	// While building, the CSR is not there yet, so these answer from what is.
	[[nodiscard]] bool is_leaf(size_t node) const { return children_of[node].empty(); }
	[[nodiscard]] bool is_root(size_t node) const {
		return tree._parents[node] == TPhylogeny::NO_PARENT;
	}

	void add_edge(const std::string &child, const std::string &parent, double branch_length);
	TPhylogeny finish();

private:
	void _reject_a_non_forest() const;
	[[nodiscard]] std::vector<size_t> _canonical_order() const;
	void _permute(const std::vector<size_t> &new_of_old);
	void _assemble();
};

void TPhylogenyBuilder::add_edge(const std::string &child, const std::string &parent,
                                 double branch_length) {
	if (child == parent) {
		throw coretools::TUserError("Node '", child, "' can not be parent of itself ! Got ", child,
		                            "for the child and ", parent, " for the parent.");
	}
	if (branch_length <= 0.0) {
		throw coretools::TUserError("You can't have a negative branch length or equal to 0.0 !");
	}

	const bool has_child  = tree.contains(child);
	const bool has_parent = tree.contains(parent);

	if (!has_child && !has_parent) {
		const size_t parent_ix = add_node(parent, TPhylogeny::NO_PARENT, 0.0);
		const size_t child_ix  = add_node(child, parent_ix, branch_length);
		children_of[parent_ix].push_back(child_ix);
	} else if (!has_child) {
		const size_t parent_ix = tree._index_by_id.at(parent);
		const size_t child_ix  = add_node(child, parent_ix, branch_length);
		children_of[parent_ix].push_back(child_ix);
	} else {
		// The child already exists, so it may only be adopted if nothing has claimed it yet.
		const size_t child_ix = tree._index_by_id.at(child);
		if (tree._parents[child_ix] != TPhylogeny::NO_PARENT) {
			throw coretools::TUserError(
			    "Node: '", child,
			    "' has already a parent in the tree. Adding an other parent is not allowed !");
		}
		const size_t parent_ix   = has_parent ? tree._index_by_id.at(parent)
		                                      : add_node(parent, TPhylogeny::NO_PARENT, 0.0);
		tree._parents[child_ix]  = parent_ix;
		length_by_node[child_ix] = branch_length;
		children_of[parent_ix].push_back(child_ix);
	}
}

/// The forest post-conditions, checked before anything is derived from the topology. They come
/// first because the canonical ordering below is a walk down from the roots: on a tree that is not
/// a forest that walk does not reach every node, and would leave nodes without an index rather
/// than reporting what is actually wrong with the file.
void TPhylogenyBuilder::_reject_a_non_forest() const {
	const size_t n = tree._ids.size();

	std::vector<size_t> roots;
	bool has_a_leaf = false;
	for (size_t i = 0; i < n; ++i) {
		if (is_root(i)) { roots.push_back(i); }
		has_a_leaf = has_a_leaf || is_leaf(i);
	}

	if (roots.empty()) {
		throw coretools::TUserError("The tree has no root: every node has a parent.");
	}
	if (!has_a_leaf) {
		throw coretools::TUserError("The tree has no leaf: every node has a child.");
	}
	for (const size_t root : roots) {
		if (is_leaf(root)) {
			throw coretools::TUserError(
			    "Node '", tree._ids[root],
			    "' is both a root and a leaf: it has no parent and no children. A root must have "
			    "at least one child.");
		}
	}

	// Reachability, which doubles as the cycle check. Every node has at most one parent, so
	// walking up from an unreachable node can only end on a cycle: had it ended on a node without
	// a parent, that node is a root and the walk was a path down from one. So there is nothing a
	// separate cycle pass would catch that this does not, and this is O(n) rather than O(n*depth).
	std::vector<bool> seen(n, false);
	std::queue<size_t> pending;
	for (const size_t root : roots) {
		seen[root] = true;
		pending.push(root);
	}
	size_t reached = roots.size();
	while (!pending.empty()) {
		const size_t node = pending.front();
		pending.pop();
		for (const size_t child : children_of[node]) {
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
				                            "' sits on or below a cycle: following its parents "
				                            "never reaches a root.");
			}
		}
	}
}

/// Where each node ends up: `new_of_old[old_index]` is the node's index in the canonical order --
/// leaves, then internal non-root nodes in post-order, then roots. See ADR-0004 for why.
///
/// File order survives inside the leaf and root blocks because both are filled by an ascending
/// walk over the nodes as built, and inside the middle block because each node's children are
/// visited in the order the file introduced them.
std::vector<size_t> TPhylogenyBuilder::_canonical_order() const {
	const size_t n = tree._ids.size();
	std::vector<size_t> new_of_old(n, TPhylogeny::NOT_IN_SPACE);

	size_t n_roots   = 0;
	size_t next_leaf = 0;
	for (size_t i = 0; i < n; ++i) {
		if (is_leaf(i)) {
			new_of_old[i] = next_leaf++;
		} else if (is_root(i)) {
			++n_roots;
		}
	}
	size_t next_internal = next_leaf;
	size_t next_root     = n - n_roots;

	// Post-order, iteratively: a taxonomy is deep enough that a recursive walk is a stack overflow
	// waiting to happen. The second element is how many of the node's children have been pushed.
	std::vector<std::pair<size_t, size_t>> stack;
	for (size_t root = 0; root < n; ++root) {
		if (!is_root(root)) { continue; }
		stack.emplace_back(root, 0);
		while (!stack.empty()) {
			auto &[node, next_child] = stack.back();
			if (next_child < children_of[node].size()) {
				const size_t child = children_of[node][next_child];
				++next_child;
				// A leaf already has its index and no descendants, so it never goes on the stack.
				if (!is_leaf(child)) { stack.emplace_back(child, 0); }
				continue; // emplace_back may have moved the vector: re-take the reference.
			}
			// Every child is done, so this is the post-visit.
			new_of_old[node] = is_root(node) ? next_root++ : next_internal++;
			stack.pop_back();
		}
	}

	// Total by construction: every leaf was numbered by the scan above, and _reject_a_non_forest
	// has already established that every other node is reachable from a root, which is exactly
	// what the walk visits. Checked anyway, because a permutation with a hole in it would show up
	// far from here as a node silently overwriting another.
	for (size_t i = 0; i < n; ++i) {
		if (new_of_old[i] == TPhylogeny::NOT_IN_SPACE) {
			throw coretools::TDevError("Node '", tree._ids[i], "' was left out of the node order.");
		}
	}
	return new_of_old;
}

/// Move every node to where `_canonical_order` put it. This is the one place the ordering is
/// applied: everything above builds in file order and knows nothing about it.
void TPhylogenyBuilder::_permute(const std::vector<size_t> &new_of_old) {
	const size_t n = tree._ids.size();

	std::vector<std::string> ids(n);
	std::vector<size_t> parents(n);
	std::vector<std::vector<size_t>> kids(n);
	std::vector<double> lengths(n);
	for (size_t old = 0; old < n; ++old) {
		const size_t at = new_of_old[old];
		ids[at]         = std::move(tree._ids[old]);
		parents[at]     = is_root(old) ? TPhylogeny::NO_PARENT : new_of_old[tree._parents[old]];
		kids[at]        = std::move(children_of[old]);
		lengths[at]     = length_by_node[old];
		for (size_t &child : kids[at]) { child = new_of_old[child]; }
	}

	// The ids themselves did not change, so the map is renumbered rather than rebuilt: rehashing
	// every node name is the expensive part, and none of the keys moved.
	for (auto &[id, index] : tree._index_by_id) { index = new_of_old[index]; }

	tree._ids      = std::move(ids);
	tree._parents  = std::move(parents);
	children_of    = std::move(kids);
	length_by_node = std::move(lengths);
}

/// Derive everything a TPhylogeny answers from the node order it now has.
void TPhylogenyBuilder::_assemble() {
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

	// Classify. Ascending node order, so each category list comes out as the contiguous block the
	// canonical ordering put it in, and each map back into one is arithmetic on the node index --
	// which is what the next ticket deletes these vectors in favour of.
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

	// One length per branch, in branch order: a root has no branch, so its entry is dropped here
	// rather than carried around as a zero.
	tree._branch_lengths.reserve(tree._branches.size());
	for (const size_t node : tree._branches) {
		tree._branch_lengths.push_back(length_by_node[node]);
	}
}

TPhylogeny TPhylogenyBuilder::finish() {
	_reject_a_non_forest();
	_permute(_canonical_order());
	_assemble();
	return std::move(tree);
}

TPhylogeny build_phylogeny(const std::vector<TEdge> &edges) {
	TPhylogenyBuilder builder;
	for (const auto &edge : edges) {
		builder.add_edge(edge.child, edge.parent, edge.branch_length);
	}
	return builder.finish();
}

TPhylogeny read_phylogeny(const std::string &filename) {
	coretools::TInputFile file(filename, coretools::FileType::Header);
	if (file.numCols() != 3) {
		throw coretools::TUserError("File '", filename, "' is expected to have 3 columns, but has ",
		                            file.numCols(), " !");
	}

	// One edge at a time, never a vector of them: the strings die at the end of each iteration.
	TPhylogenyBuilder builder;
	for (; !file.empty(); file.popFront()) {
		builder.add_edge(std::string(file.get(0)), std::string(file.get(1)), file.get<double>(2));
	}

	return builder.finish();
}
