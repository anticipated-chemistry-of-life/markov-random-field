# Nodes are stored in canonical order, so a node's index tells you what kind of node it is

`TPhylogeny` stores its nodes as three contiguous blocks: **all leaves, then all internal non-root nodes in post-order, then all roots.** Within the leaf block and the root block the nodes keep the relative order the tree file put them in; within the middle block the order is post-order from each root, visiting children in file order. The construction code is unaware of this — it builds in file order and one permutation is applied at the end.

The point is that the layout answers, by arithmetic alone, the three questions the sweep asks most often. A node is a leaf exactly when its index is below the leaf count. A node is a root exactly when its index is at or above the node count minus the root count. A leaf's index in leaf space *is* its node index; a branch's index *is* its child node's index; an internal node's index in internal-node space is its node index minus the leaf count. All of it derives from three integers — node count, leaf count, root count.

That is a structural claim standing in for a structural test, and it is the reason this record exists: a reader who finds `node < _n_leaves` where they expected "does this node have children" needs to know why the comparison is legal.

**The three blocks partition the nodes only because a root is required to have at least one child.** A childless node would be both a leaf and a root and would belong to two blocks at once. That rule is not decoration — it is enforced during construction, and this layout is what it is enforced *for*.

## What it buys

- Ten parallel vectors — five mapping node index forward into a category list, five mapping back — collapse into three counts. The sentinel that meant "not applicable" stops existing, so an unsigned index can no longer wrap around into something that reads like a valid one.
- The hot accessors stop touching memory to answer a question the index already answers. Asking whether a node is a leaf used to mean a bounds-checked index into the node array, a pointer chase into that node's child list to test whether it was empty, and then a second random lookup into a category index vector. It is now an integer comparison.
- Index-space conversions become the identity or a subtraction, so indexing the field, the internal state and the branch lengths needs no lookup table.
- Bottom-up traversal becomes linear. Children always precede parents, so initialising internal state from children is a single forward pass over the middle block and then the roots, rather than a fixed-point loop that re-sweeps a set of not-yet-ready nodes and is quadratic in the worst case.
- Any future bottom-up walk gets children-before-parents for free, as a property of the storage rather than something rediscovered at runtime.

## What it costs

- **Output row ordering changes relative to what the code produced before.** Anything written in node, leaf, internal-node or branch order comes out permuted — and, because branch indices change, the sampler forms different branch pairs, so the random stream changes too. Two runs either side of this decision are not comparable byte for byte, and a diff between them says nothing.
- **Reverting is expensive.** It means reintroducing the ten vectors and every call site that read them, including the ones in other classes that were deleted precisely because the ordering made them unnecessary.
- The conversions stay named methods on `TPhylogeny` rather than arithmetic inlined at call sites, so that the invariant keeps exactly one home. That is a small ongoing tax on anyone who would rather write the subtraction where they need it.

## Considered options

**Sorting node ids alphabetically** would be stabler across edits to an input tree and nicer for diffing two runs. It was rejected because it reorders output relative to file order for no functional gain: the arithmetic above needs the three blocks to be contiguous, and it does not care what the order *within* a block is. File order within each block is the decision, so a user who reorders the lines of their tree gets a correspondingly reordered output rather than an unrelated one.

**Keeping the index vectors and merely reordering** would have delivered the traversal win without touching the bookkeeping. It was rejected as a stopping point rather than as an option: the vectors are the thing worth removing, and the ordering is what makes removing them possible. The two are nevertheless split across two commits, because reordering changes the chain and can only be checked statistically, while deleting the vectors changes nothing and can be checked exactly.

## Consequences

**Previously written internal-state files keep working.** The reader that loads them resolves each row by node name rather than by raw linear index — a change made in phase 1 (#19) specifically so that this decision would need no migration. Node names do not change here. Without that change this record would have had to specify one.

**The independent validation harness replicates this ordering.** `model_validation/src/independent/indexing.py` derives leaf order, internal order and branch order the same way the C++ does, because it writes and reads the same positional files. It is not documentation of the ordering; it is a second implementation of it, and it has to move whenever this one does.

**The ordering is asserted, not assumed.** The property-based tests in `tests/TPhylogeny_Tests.cpp`, over the forests `tests/phylogeny_generators.h` generates, state each guarantee as something a caller can observe — every node in the middle block has a parent whose index is greater than its own, the three blocks partition the nodes, the conversions round trip — over multi-root forests, deep chains and wide stars, not only over balanced trees. A balanced fixture satisfies post-order by accident, which is exactly the way this would fail silently.
