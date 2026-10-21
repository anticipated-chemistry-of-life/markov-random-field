# Leaves are numbered depth-first, so siblings are adjacent

ADR-0004 fixed the three blocks of the canonical order — leaves, internal non-root nodes in post-order, roots — and left the order *inside* the leaf block to the file. That order does not matter to the arithmetic, and it turned out to matter to the cache.

The block update reads, for every leaf pair, the node state of the species leaf's parent. Since ADR-0013 the species node state keeps its cliques contiguous, so the parents of a band of 64 species leaves are read out of one column, and how many cache lines that costs is decided by how far apart those parents' indices are. In `inference_lotus/species.tsv` the leaves are listed in an order unrelated to the tree: along the 32,347 leaves the parent changes 32,312 times, and a band of 64 leaves has 62 distinct parents. Every parent read was a cache line of its own.

The leaves are now numbered by the same depth-first walk that numbers the internal nodes. When the walk first reaches a node it numbers that node's leaf children together, in file order, and then descends into its internal children in file order. So the leaves below any node are a run of consecutive indices, and a node's own leaf children are adjacent within it. On the LOTUS species tree a band of 64 leaves now has 16.8 distinct parents on average, and the parent changes 8,129 times along the leaves. Roots keep their file order, and the middle block is the post-order it always was.

## What it costs

ADR-0004 chose file order inside the leaf block so that "a user who reorders the lines of their tree gets a correspondingly reordered output rather than an unrelated one". That is given up for the leaves. The leaf order is now a function of the tree's shape and the order the file lists each node's children in. Reordering the lines of a tree file still reorders the output, but no longer line for line.

- **A different chain.** A cell's linear index is built from leaf indices, and its uniform is hashed from the linear index (ADR-0007). Runs either side of this change are equally valid draws from the same posterior; they are not comparable byte for byte.
- **Output row order.** Everything written in leaf order — the field posterior, the node-state files — comes out in the new order. The node-state and data readers resolve leaves by name and are unaffected. A downstream script that lined up rows with the tree file by position has to join on name.
- **Field files written before this change.** The field reader (`--set_Y`, `TMarkovField::_read_Y_from_file`) goes by the `position` column alone, so a field file written by an earlier binary, or by the validation harness before it moved, means something else under the new order and has to be written again.
- **The validation harness moves with it.** `model_validation/src/independent/indexing.py` reimplements the order (ADR-0004), and does it the same way. One small tree on which the two orders differ is pinned on both sides, `PhylogenyOrder.leaves_are_depth_first_with_siblings_adjacent` and `test_leaves_are_numbered_depth_first_with_siblings_adjacent`, so neither can move alone.

## Considered options

**Sorting the input instead.** A preprocessing step that rewrites the tree file in depth-first order gives the same speed with no change to the code, and keeps ADR-0004's promise intact. It was rejected because the speed would then depend on every user remembering to run it, on every tree, and the gain is not a property of anybody's data — it is what any tree file gets.

**Ordering only the species tree.** The block update's stride problem is the species tree's. But the walk is one function for both trees, the molecule tree gains too (below), and two orders would be one more thing to explain.

**Sorting siblings by subtree size, or by any key other than file order.** Any order that keeps each clade a run gives the same cache behaviour. File order among siblings keeps as much of the old promise as survives.

## Evidence

`bench/dense_storage_bench.cpp`, `own_clique`, on the full LOTUS trees, 32 threads, a shared machine. With `species.tsv` in its file order and the leaves in file order: block update 5.84 / 5.57 s, iteration 13.66 / 13.34 s. With a copy of `species.tsv` rewritten depth-first (the test before the code changed): block update 4.22 / 4.24 s, iteration 11.58 / 11.88 s — even though those runs saw the higher load. With this change and the unmodified file: block update 4.12 / 4.16 s, molecule walk 1.18 s → 0.74 s, iteration 11.04 / 11.02 s.

The production binary before and after, profiled with gperftools over a complete chain of 3 iterations on `inference_lotus/` (`--numBurnin 0 --thinning 1 --fixedSeed 42 --numThreads 32`), back to back on a shared machine: the whole run 1026 → 953 CPU s, the chain 3 min 4 s → 2 min 35 s, though the second run started at a load of 38 against 12. `TPhylogeny::parent_of`, where the parent misses were charged, went from 27 CPU s to below the profile's cut-off. The block update as a whole did not move (263.5 → 263.6 CPU s): in production it is computation, not memory — `draw_block`, `block_probabilities` and `states_of_index` are half of it — and the benchmark's reduced model, whose data factors are all 1, is what made the misses its largest cost. The clique traversals went 258 → 229 CPU s.

The chain start needed one fix to show this. `leaf_layer_start::hold_tree_fields_at_the_field` copied the field into both tree fields row by row and by linear index, which under ADR-0013's layout is a division and a cache line per cell of the species node state: 17 CPU s on the row-major binary, 48–115 s after ADR-0013. It now visits the field in bands of 64 species leaves, column by column, by subscript: 24 CPU s.
