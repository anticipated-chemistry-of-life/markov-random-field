# Each tree draws its own leaf layer, and the field draws its own

_The **decision** below is superseded by [ADR-0010](0010-the-block-draws-the-leaf-layer.md), which reinstates the block update. **The measurement is not superseded.** It is the evidence ADR-0010 stands on, and the `reference` column of the table in the consequences is the binary that runs today: reverting this record's commit leaves `src/` and `tests/` byte-identical to `a3d1ae6`. So read the decision and the considered options as the argument that was made, and the measurement as the price the decision was reversed on._

ADR-0005 gave each tree a leaf-level field of its own, and it argued for one **eight-state block update** over the triple `(Y, Z_s, Z_m)` at a leaf pair. This record retires the block. The leaf layer is drawn **one variable at a time**: each tree draws every node of its node state, its leaf block included, and the field has an update of its own.

ADR-0005 said the failure "would present as slow mixing rather than as a bug", which is a claim about a number nobody had. Issue #69 measured it. This record carries that measurement as its consequences, so the trade is priced here rather than argued here.

## What ADR-0005 keeps

The block was one paragraph of ADR-0005's consequences and not its argument. Everything the argument reached still stands, unchanged and still read from that record:

- **The factorisation.** `p(Z_s | theta_s) * p(Z_m | theta_m) * p(Y | Z_s, Z_m, omega) * p(L, D | Y)`. Every factor is a proper conditional density, `C` is identically 1, and ADR-0002's bias is still structurally absent. Nothing here touches the model.
- **The link table.** The AND over two independently corrupted tree fields, behind the same link-policy seam, with the same two parameter-free constraints traced as diagnostics.
- **The error probability's support.** Still the open interval `(0, 0.5)`, and still for both of ADR-0005's reasons. One of the two has to be re-read: ADR-0005 wrote the lower bound as "the link is the deterministic AND, the block update hits `log 0`, and a field cell at one pins both tree fields". There is no block to hit `log 0` now, but `log P_0 = 2 log omega` is the term a single-variable draw evaluates as well, so the bound holds for the same arithmetic. The upper bound is untouched: above `0.5` the tree fields are anti-correlated with the field, which is a genuine second mode whatever draws it.
- **The six-counter collapse.** The link's likelihood is still a function of `omega` and `n(bucket, field state)`, and the error probability's move still costs the same whatever the size of the field. What changed is who tallies them, which is below.
- **The tree-field-versus-field distinction.** A tree field is the leaf block of one tree's node state; the field is the noisy reconciliation of the two. That is the distinction this record acts *on*: the block was the one place the two were drawn by one owner.

ADR-0005's derivation 3 also survives with more force than before, because the measurement leans on it: the field's density is the product `a~_s * a~_m`, so a chain riding that ridge moves the two tree fields in opposite directions while the field's own posterior stands still. That is why the measurement below reads both tree field posteriors and not the field's.

## The decision

**The node-state walk covers every node.** A root is drawn from the stationary distribution and its children. An internal node is drawn from its parent and its children. A leaf is drawn from its parent and from the link. `node_state_walk::update_clique`, `src/tree/`.

**The field has a pass of its own.** It visits every leaf pair, draws that cell from the two tree field cells at that pair and from the data that observes the field, and retallies the six link counters as it goes. `field_update::run`, `src/field/`.

**An iteration runs species node state, molecule node state, then field.** The field update therefore reads two tree fields that have already moved this iteration, and the six counters it leaves describe the configuration the error probability then proposes against.

The reason is ownership, and it is ADR-0005's own reason turned one step further. ADR-0005 made a tree field the leaf block of that tree's node state. The leaf layer was then the one part of a node state that its own tree did not draw, and the block was a second owner for a variable that already had one. A tree that cannot draw its own leaves is a tree whose update stops for a reason no reader of ADR-0005 would predict.

**A leaf reaches the link through a seam the caller binds.** The walk asks for `{ P(Y | this leaf = 0), P(Y | this leaf = 1) }` at the leaf pair the leaf occupies, and the tree forwards what its caller bound and looks behind it nowhere. `TMarkovField` binds one `tree_field_link::TLeafLinks` per tree, and hands each tree the *other* tree's node state inside it. So a tree still names neither the field, nor the other tree, nor the error probability — which is what lets the walk run in a test with none of the three. `src/field/tree_field_link.h`.

That seam is the whole reason this is not a layering violation. A leaf is the one node that is *not* conditionally independent of everything outside its own tree, and the seam is where that one exception is stated.

## The leaf's uniform moved streams

ADR-0007 derives a cell's uniform from the seed, the stream, the tree, the iteration and the cell's linear index. The block update drew from `TCellStream::field`, at the field's linear index, and it drew once for the whole triple. A leaf now draws from `TCellStream::node_state`, at that leaf's linear index in its own tree's node-state container space, with the tree in the key. The field draws its own pass on `TCellStream::field`.

**The values change and the guarantee does not.** A cell's uniform is still a function of the cell's position and of nothing that counts the draws before it, so a cell still gets the same number whichever thread reaches it and whichever order the update visits the container in. What moved is which position names it. ADR-0007's property is the one thing the parity gate rests on, and this record leaves it exactly as it was: the gate is byte-identical across the change.

Two things follow, and both are stated so that an expected difference is not read as a bug.

The chain moves. Three single-variable draws replace one eight-state block, and a leaf takes a different number from a different stream, so no trace from before this change compares with one from after it. The target distribution is untouched. The chain says so in its log file.

And **the two revisions cannot be compared seed for seed.** One seed drives two different sequences of draws, so a seed-matched pair would differ for a reason that says nothing about mixing. The measurement below therefore fixes the data set and lets the seed go.

## What `--fix_Y` and `--fix_Z` now mean

Both are held by their own flag, `--Y.update false` and `--Z.update false`.

**`--fix_Y` holds the field alone.** Under the block it held the whole triple, because one draw owned all three. Both tree fields now move, each drawn by its own tree, so the bucket a leaf pair falls in has moved even where its field state has not. The tally the chain start built no longer describes the configuration, and the error probability is scored against six counters. So `field_update::tally` recounts all six every iteration, in place of the pass that used to build them as it went. A field that is empty *and* fixed is a user error, and the run stops before the log line reports a start.

**`--fix_Z` holds every node of both node states, leaves included.** A tree field is the leaf block of a node state, so fixing a node state fixes its tree field with it. There is no flag that holds the nodes above the leaves and lets the leaf block move; that would be a third meaning for a walk that has one.

**Passing both reproduces the old `--fix_Y`.** The field, both tree fields and every node above them stand still, which is the whole leaf layer and more. A run that wants the old flag's behaviour passes both, and the constructor's log line says so.

## Considered options

**Keep the block, and let each tree draw its internal nodes only.** This is the arrangement ADR-0005 left, and it is the one this record replaces. It is the best of these options for mixing and the worst for shape: the walk stops above the leaves for a reason internal to another update, `TCliqueView` has to refuse a leaf write, and a leaf's conditional lives in a third place that names the field, both trees and the error probability at once. Issue #57 recorded the trade as simplicity for mixing, deliberately, and the measurement below is what it costs.

**Fuse the two leaf blocks and the field into a single pass.** Not rejected — deferred. A quarter of the measured loss is not mixing at all, it is three passes over the leaf-pair space where there was one. This keeps the single-variable draws and gives that quarter back, and it is repair 2 of the four below.

**Put a block move back as a proposal rather than as the update.** Also deferred, and it is the option to reach for first if the cost is judged unacceptable. A Metropolis move that flips the triple at one leaf pair sits *beside* the single-variable walk. It needs no eight-state enumeration, no third traversal and no model seam, and it targets exactly the trap. Repair 3.

**Reinstate the block.** ADR-0005's own answer, and the most complexity for the most mixing. The measurement below prices what it buys, which is what makes this a decision rather than a preference. Repair 4.

## Consequences: the measured cost

`model_validation/mixing_cost/findings.md` is the write-up and `model_validation/mixing_cost/run.sh` reproduces it. **The cost is real, it is confined to the field, and it is a small-`omega` effect.**

Reference `a3d1ae6`, the commit before the block was deleted, against current `f13fef9`, the commit that deleted it. Both link the same coretools and stattools checkout. One data set, simulated by the reference binary on a 128 x 256 fixture — 32 768 leaf pairs — 20 000 iterations after 2 000 burn-in, `--numThreads 1`, four replicates per binary at `omega = 0.005` and two at `omega = 0.1`.

The joint density trace is read factor by factor, and that is what decides the result. Integrated autocorrelation time at `omega = 0.005`:

| factor                       | reference `tau` | current `tau` |    ratio |
| ---------------------------- | --------------: | ------------: | -------: |
| `data` — the field's own     |            63.2 |         295.1 | **4.67** |
| `link`                       |            1371 |          1604 |     1.17 |
| `species_node_state`         |            3630 |          3718 |     1.02 |
| `molecules_node_state`       |            1182 |          1377 |     1.17 |
| `joint_density` — the total  |            2295 |          2829 |     1.23 |

**Read the `data` row and not the total.** `data` is `log p(L, D | Y)`: a function of the field and of the two data-source parameters, and of nothing either phylogeny carries. It is the field's own instrument, and it is the only factor whose drift says the chain had reached the distribution its autocorrelation time describes. The two node-state factors are dominated by the phylogenetic parameters, which had not converged in either binary at 20 000 iterations — they drift about one trace standard deviation — and they dilute the total. The total's 1.23 is that dilution and not a finding.

The shape is the metastability ADR-0005 described. `acf(1)` on the field's factor goes 0.567 to 0.934, and at lag 200 the current binary is still 0.233 correlated where the reference has reached 0.052.

**The iteration also costs more, so the two costs compound.** 11.20 ms against 8.95 ms, a factor of 1.25, because three passes over the leaf-pair space cost more than one eight-state enumeration. Together that is **5.9 times fewer effective field draws per second of wall clock**.

At `omega = 0.1` the field's ratio is **1.23**, and the two autocorrelation functions are within noise of each other by lag 10. **This is the regime dependence ADR-0005 predicted, confirmed.** The trap is the deterministic AND at `omega = 0`, and it loosens as the link is allowed to be wrong more often. The block was worth its complexity exactly where `omega` is small.

**Both tree field posteriors agree.** Pooled over replicates, mean absolute difference 0.004 per cell and correlation 0.9998 for each tree. On ADR-0005's ridge, current against reference is −0.00087 **along** it and +0.00079 **across** it. The along shift is the coordinate the metastability would move, and it sits well inside either binary's own replicate scatter — 0.0025 for the reference, 0.0043 for the current — so neither rides the ridge. The across shift is about three times the reference's own spread, so it is not simply noise; in absolute terms it is a 0.08% difference in the field's density. Small, and recorded as a difference rather than as a zero.

**What the measurement does not settle**, stated so that nobody reads it for more than it says. The phylogenetic parameters had not converged in either binary, so the tree field posteriors are compared at a point the chain was still moving through. The `omega` dependence has two points and not a crossover. And at `omega = 0.005` the current binary's field factor leaves an effective sample of about 68 per chain, so the 4.7 is an order of magnitude and not three digits.

**If the cost is judged unacceptable**, four things could be done, in rising order of what they cost to build: run longer; recover the iteration time by fusing the three passes; add a block Metropolis proposal beside the single-variable update; reinstate the block. The write-up ranks them and says what each buys. The price this record accepts is bounded: it is confined to the field, it is a factor of about 5 at the pessimistic end of `omega`, it shrinks to about 1.2 as `omega` grows, and it does not move the answer.

## Other consequences

**A tree's update writes leaves now, so the clique view's write assertion inverted.** `TCliqueView` refused a leaf write, because a leaf was the block's. The view takes a write policy instead, and the chain start is now the one walk that stops above the leaves.

**The field update's model answers for the data alone.** The two tree factors left with the block, so what the field's model scores is `log p(L, D | Y)` at one cell. `TFieldModel`, `src/field/`.

**The eight-state math is deleted, and the link is untouched.** `TBlockUpdate` and its model go. The buckets, the six counters, the AND policy, the per-bucket probabilities and the diagnostic all stay exactly where ADR-0005 put them.

**The parity gate is byte-identical across the change**, and it stays the only thing that asserts two backends compute one chain.

**`--K` is gone.** It sized the sheet the field update cached under the block, and the field update reads a whole row at a time now. A run that still passes it stops with a message rather than reporting an unused argument at the end.

**CONTEXT.md keeps *Block update* as a retired term.** It names what ADR-0005 built and this record retired, so a reader who finds the phrase in an old commit or an old branch can still resolve it.
