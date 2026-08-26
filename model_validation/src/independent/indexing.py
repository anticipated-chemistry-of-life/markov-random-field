"""Node ordering, replicated from the C++ tree loader.

The C++ assigns every node an index by *file appearance order*, with no sort:
reading each `child<TAB>parent` row, it appends the parent first if unseen, then
the child (`TPhylogenyBuilder::add_edge` in `src/tree/TPhylogeny.cpp`).
Three separate orderings are then derived from that node order, and the Python
must agree with all three or every downstream file is silently misaligned:

- **leaf order**: leaf nodes, in node order. Indexes the field's dimensions.
- **internal order**: non-leaf nodes *including roots*, in node order
  (the classification loop in `TPhylogenyBuilder::finish`). Indexes the
  internal-state dimension.
- **branch order**: non-root nodes, in node order. Indexes branch lengths.

Note the two easy mistakes: "internal" includes roots, and "branch" interleaves
leaves with internal non-roots rather than listing one group after the other.
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


def build_tree_index(edges: list[tuple[str, str]]) -> TreeIndex:
    """Build a `TreeIndex` from `(child, parent)` rows in tree-file order."""
    order: list[str] = []
    index_of: dict[str, int] = {}

    def intern(name: str) -> int:
        if name not in index_of:
            index_of[name] = len(order)
            order.append(name)
        return index_of[name]

    parent = {}
    for child, par in edges:
        # The C++ appends the parent before the child when both are new.
        par_ix = intern(par)
        child_ix = intern(child)
        if child_ix in parent:
            raise ValueError(f"Node {child!r} has more than one parent.")
        parent[child_ix] = par_ix

    n = len(order)
    parent_arr = np.full(n, -1, dtype=np.int64)
    for child_ix, par_ix in parent.items():
        parent_arr[child_ix] = par_ix

    # A node is a leaf iff it never appears as a parent.
    has_child = np.zeros(n, dtype=bool)
    has_child[parent_arr[parent_arr >= 0]] = True

    leaves = np.flatnonzero(~has_child)
    internals = np.flatnonzero(has_child)
    branch_nodes = np.flatnonzero(parent_arr >= 0)

    branch_of_node = np.full(n, -1, dtype=np.int64)
    branch_of_node[branch_nodes] = np.arange(len(branch_nodes))

    # Parents are appended before their children, so node order is already
    # topological. The sampler relies on this; assert it rather than assume it.
    if not np.all(parent_arr < np.arange(n)):
        raise ValueError("Node order is not topological: a child precedes its parent.")

    return TreeIndex(
        names=order,
        parent=parent_arr,
        leaves=leaves,
        internals=internals,
        branch_nodes=branch_nodes,
        branch_of_node=branch_of_node,
    )
