"""Node ordering, replicated from the C++ tree loader.

The C++ stores nodes in a **canonical order** — all leaves, then all internal
non-root nodes in post-order, then all roots — rather than in the order the tree
file lists them. See `docs/adr/0004-nodes-are-stored-in-canonical-order.md`; the
permutation itself is `TPhylogenyBuilder::_canonical_order` in
`src/tree/TPhylogeny.cpp`. Within the leaf block and the root block the nodes
keep their relative *file* order, which is what this module has to reproduce
first: reading each `child<TAB>parent` row, the loader appends the parent if
unseen and then the child (`TPhylogenyBuilder::add_edge`).

This is not documentation of the C++ ordering, it is a second implementation of
it. Every positional file the harness writes or reads is addressed in one of the
three orderings below, so the two have to agree or every downstream file is
silently misaligned:

- **leaf order**: the leaf block, so leaf index *is* node index. Indexes the
  field's dimensions.
- **internal order**: the internal non-root block followed by the root block, so
  internal index is node index minus the leaf count. Indexes the internal-state
  dimension. Note that "internal" includes roots.
- **branch order**: the leaf block followed by the internal non-root block — every
  non-root node — so branch index *is* node index. Indexes branch lengths.

The consequence worth holding on to is that a child's index is always below its
parent's, the opposite of what file order gave.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class TreeIndex:
    """A tree's topology under the C++ node ordering."""

    names: list[str]
    """Node names, in C++ node order."""

    parent: np.ndarray
    """`parent[n]` is the node index of `n`'s parent, or -1 if `n` is a root."""

    leaves: np.ndarray
    """Node indices of the leaves, in C++ leaf order."""

    internals: np.ndarray
    """Node indices of the non-leaf nodes (roots included), in C++ internal order."""

    branch_nodes: np.ndarray
    """Node indices of the non-root nodes, in C++ branch order."""

    branch_of_node: np.ndarray
    """`branch_of_node[n]` is `n`'s index in `branch_nodes`, or -1 for a root."""

    @property
    def n_nodes(self) -> int:
        return len(self.names)

    @property
    def n_leaves(self) -> int:
        return len(self.leaves)

    @property
    def n_internals(self) -> int:
        return len(self.internals)

    @property
    def n_branches(self) -> int:
        return len(self.branch_nodes)

    def leaf_names(self) -> list[str]:
        return [self.names[i] for i in self.leaves]

    def internal_names(self) -> list[str]:
        return [self.names[i] for i in self.internals]

    def branch_names(self) -> list[str]:
        return [self.names[i] for i in self.branch_nodes]


def _file_order(edges: list[tuple[str, str]]) -> tuple[list[str], dict[int, int]]:
    """Node names as the C++ creates them, and `child -> parent` over their indices."""
    order: list[str] = []
    index_of: dict[str, int] = {}

    def intern(name: str) -> int:
        if name not in index_of:
            index_of[name] = len(order)
            order.append(name)
        return index_of[name]

    parent: dict[int, int] = {}
    for child, par in edges:
        # The C++ appends the parent before the child when both are new, and after
        # it when the child is a node some earlier row already created. Interning
        # the parent first reproduces both cases: in the second the child is
        # already in, so it keeps the lower index either way.
        par_ix = intern(par)
        child_ix = intern(child)
        if child_ix in parent:
            raise ValueError(f"Node {child!r} has more than one parent.")
        parent[child_ix] = par_ix

    return order, parent


def _canonical_order(parent: np.ndarray, children: list[list[int]]) -> np.ndarray:
    """`new_of_old[i]`: where node `i` of the file order lands in canonical order.

    Leaves keep their relative file order and come first; then the internal
    non-root nodes in post-order from each root, visiting children in file order;
    then the roots, again in their relative file order.
    """
    n = len(parent)
    is_leaf = [len(kids) == 0 for kids in children]
    is_root = parent < 0

    new_of_old = np.full(n, -1, dtype=np.int64)
    next_leaf = 0
    for node in range(n):
        if is_leaf[node]:
            new_of_old[node] = next_leaf
            next_leaf += 1

    next_internal = next_leaf
    next_root = n - int(is_root.sum())

    for root in range(n):
        if not is_root[root]:
            continue
        # Iterative post-order: `(node, how many of its children have been pushed)`.
        stack: list[list[int]] = [[root, 0]]
        while stack:
            node, seen = stack[-1]
            if seen < len(children[node]):
                stack[-1][1] += 1
                child = children[node][seen]
                # A leaf is already numbered and has no descendants.
                if not is_leaf[child]:
                    stack.append([child, 0])
                continue
            if is_root[node]:
                new_of_old[node] = next_root
                next_root += 1
            else:
                new_of_old[node] = next_internal
                next_internal += 1
            stack.pop()

    if (new_of_old < 0).any():
        missing = [i for i in range(n) if new_of_old[i] < 0]
        raise ValueError(
            f"Nodes {missing[:5]} are not reachable from a root: not a forest."
        )
    return new_of_old


def build_tree_index(edges: list[tuple[str, str]]) -> TreeIndex:
    """Build a `TreeIndex` from `(child, parent)` rows in tree-file order."""
    order, parent_of = _file_order(edges)
    n = len(order)

    parent_arr = np.full(n, -1, dtype=np.int64)
    for child_ix, par_ix in parent_of.items():
        parent_arr[child_ix] = par_ix

    # Children in the order the file introduced them, which is what the post-order
    # walk below visits them in (`children_of[parent].push_back(child)` per edge).
    index_of = {name: i for i, name in enumerate(order)}
    children: list[list[int]] = [[] for _ in range(n)]
    for child, par in edges:
        children[index_of[par]].append(index_of[child])

    new_of_old = _canonical_order(parent_arr, children)

    names = [""] * n
    parent = np.full(n, -1, dtype=np.int64)
    for old in range(n):
        at = new_of_old[old]
        names[at] = order[old]
        parent[at] = -1 if parent_arr[old] < 0 else new_of_old[parent_arr[old]]

    n_leaves = int(np.sum([len(kids) == 0 for kids in children]))
    n_roots = int(np.sum(parent_arr < 0))
    leaves = np.arange(n_leaves, dtype=np.int64)
    internals = np.arange(n_leaves, n, dtype=np.int64)
    branch_nodes = np.arange(n - n_roots, dtype=np.int64)

    branch_of_node = np.full(n, -1, dtype=np.int64)
    branch_of_node[branch_nodes] = branch_nodes

    # Post-order storage puts every child below its parent. The sampler relies on
    # this -- it walks the nodes backwards -- so assert it rather than assume it.
    non_roots = parent >= 0
    if not np.all(parent[non_roots] > np.arange(n)[non_roots]):
        raise ValueError("Node order is not post-order: a parent precedes its child.")

    return TreeIndex(
        names=names,
        parent=parent,
        leaves=leaves,
        internals=internals,
        branch_nodes=branch_nodes,
        branch_of_node=branch_of_node,
    )
