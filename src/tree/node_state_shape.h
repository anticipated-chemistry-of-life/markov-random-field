//
// The container space one tree's node state occupies.
//

#pragma once

#include "constants.h"
#include "tree/TPhylogeny.h"
#include <cstddef>

/// The field's leaf-space shape with `dimension`'s extent replaced by the node count of the tree
/// that owns it: every node, leaves included, against the other trees' leaves (ADR-0005).
///
/// This lives in one place because two callers need it and they must not drift. TTree sizes the
/// real node state from it, and the storage conformance tests take the shapes they assert over
/// from it. A test that restated the rule instead of calling it would keep passing after the rule
/// changed, which is exactly what it is there to catch.
///
/// Only the owning dimension moves. That is what keeps the node state in the field's orientation,
/// so one leaf pair is at the same (row, column) in both, and -- because a leaf's index in leaf
/// space is its node index (ADR-0004) -- the leaf block needs no conversion at all.
[[nodiscard]] inline IndexArray node_state_dimensions(IndexArray leaf_counts_per_tree,
                                                      size_t dimension,
                                                      const TPhylogeny &tree_of_that_dimension) {
	leaf_counts_per_tree[dimension] = tree_of_that_dimension.n_nodes();
	return leaf_counts_per_tree;
}
