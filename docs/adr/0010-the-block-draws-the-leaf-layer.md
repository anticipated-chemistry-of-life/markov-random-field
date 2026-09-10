# Each tree stops above its leaves, and the block draws the leaf layer

ADR-0008 retired the eight-state **block update** and drew the leaf layer one variable at a time. It also named what it did not know: "if the cost is judged unacceptable", four things could be done about it. Issue #69 measured the cost. This record is that judgement.

**The block comes back as the leaf-layer update.** The field, the species tree field and the molecule tree field at one leaf pair are drawn from all eight combinations at once, and a tree's node-state walk stops above its leaves.

This record supersedes the **decision** of ADR-0008 and none of its **measurement**. That table is the evidence this record stands on, and it is still read from ADR-0008 and from `model_validation/mixing_cost/findings.md`.

## What the measurement says

Integrated autocorrelation time by factor of the joint density, at `omega = 0.005`, four replicates per binary on a 128 x 256 fixture. `block` is ADR-0005's update, `single-variable` is ADR-0008's, and the ratio is the second against the first, as recorded.

| factor                      | `block` `tau` | `single-variable` `tau` |    ratio |
| --------------------------- | ------------: | ----------------------: | -------: |
| `data` — the field's own    |          63.2 |                   295.1 | **4.67** |
| `link`                      |          1371 |                    1604 |     1.17 |
| `species_node_state`        |          3630 |                    3718 |     1.02 |
| `molecules_node_state`      |          1182 |                    1377 |     1.17 |
| `joint_density` — the total |          2295 |                    2829 |     1.23 |

**Read the `data` row and not the total.** `data` is `log p(L, D | Y)`, a function of the field and of the two data-source parameters and of nothing either phylogeny carries. It is the field's own instrument, and it is the only factor whose drift says the chain had reached the distribution its autocorrelation time describes. The two node-state factors are dominated by phylogenetic parameters that had not converged in either binary, and they dilute the total.

The shape is the metastability ADR-0005 described. `acf(1)` on the field's factor goes 0.567 to 0.934: one pass of single-variable draws moves the field barely at all. At lag 200 the single-variable chain is still 0.233 correlated where the block has reached 0.052.

**The iteration also costs more, so the two costs compound.** 8.95 ms against 11.20 ms, a factor of 1.25, because three passes over the leaf-pair space cost more than one eight-state enumeration. Together that is **5.9 times fewer effective field draws per second of wall clock**.

**Both binaries agree about where the answer is.** Mean absolute difference 0.004 per cell on each tree field posterior, correlation 0.9998, and on ADR-0005's ridge a shift of −0.00087 along it against a replicate scatter of 0.0025 to 0.0043. So this is a decision about how fast the field is explored, not about what the chain converges to. That is what makes it a trade and not a bug.

## Why the small-`omega` regime decides it

At `omega = 0.1` the field's ratio is 1.23 against a 1.27 iteration cost, which is a wash, and the two autocorrelation functions are within noise of each other by lag 10. The trap is the deterministic AND at `omega = 0`, and it loosens as the link is allowed to be wrong more often.

**This model lives at the small end.** The default error probability is 0.005 (`cli.h`), under an exponential prior of rate 200 — a prior mean of 0.005. The chain is not merely started in that regime; the prior concentrates it there. So the operative number is 4.67 and not 1.23, and a leaf layer that mixes badly in the interval the prior concentrates on is a leaf layer that mixes badly.

**The condition for revisiting this record is therefore an `omega` regime an order of magnitude above the prior mean.** Real LOTUS data may put it there. If it does, the block costs about a quarter of the iteration time for a mixing gain of about the same size, and repairs 1 to 3 of `findings.md` are back on the table. The `omega` dependence is measured at two points, so where the crossover sits is not known.

## The decision

**The block draws the triple.** One thread takes a species leaf and walks the molecule leaves of its row. Per leaf pair it reads five cells — the three the draw moves and the two tree parents — draws all three from the eight combinations at once, writes back the ones the draw changed, and adds the pair to the six link counters. `block_update::run`, `src/field/TBlockUpdate.h`. The arithmetic stays in `field/TFieldMath.h`.

**The rows are conditionally independent.** A leaf pair's Markov blanket holds no cell of another leaf pair: it holds the two tree parents and the two data terms. A leaf is never a root and an internal node is never a leaf, so no thread writes a cell another row reads.

**A tree's walk covers roots and internal nodes, and `TCliqueView` refuses a leaf write.** The leaf layer is not a node state's own to draw.

**An iteration runs the block, then each tree's nodes above its leaves, then the parameters.** The six counters the block leaves describe the configuration the error probability then proposes against.

**A leaf pair takes one uniform, from the field stream at the field cell's linear index.** ADR-0007's property is untouched: a cell's uniform is a function of its position and of nothing that counts the draws before it. The values differ from a run under ADR-0008, where a leaf drew from its own tree's node-state stream, so no trace compares across either change.

**`--fix_Y` holds the field and both tree fields.** One draw owns all three, so holding one holds them all, and the tally the chain start built stands for the whole chain. **`--fix_Z` holds the nodes above the leaves**, the leaf layer staying the block's. This reverses #57's user stories 8 to 10, and the constructor's log line says what each flag holds.

**The block writes through `write_or_defer(storage, index, state, inserts)`, which locates the cell for it.** ADR-0009's handle is still the one write path and its lifetime convention still holds — the locate and the write are one call apart here. The block reads a cell with `is_one` and locates it again only if the draw moved it, which is a second lookup on the cells that changed. Under a hash map that is a hash, and the 8.95 ms iteration measured above already includes it.

## What ADR-0008 was right about

ADR-0008's argument was **ownership**, not mixing. This record does not answer that argument. It accepts it, and pays in three places:

1. **A tree's walk stops above its leaves**, for a reason internal to another update. ADR-0005 makes a tree field the leaf block of a node state, so a reader of that record would not predict the stop.
2. **`TCliqueView` refuses a leaf write**, which is a rule about a caller the clique view does not name.
3. **The leaf's conditional lives in a third place that names everything at once.** The block reads the field, both node states, both phylogenies and the error probability. ADR-0008's `tree_field_link` seam spread that dependency thinner, and it is deleted.

Nothing in that list is wrong, and none of it is cheap. What ADR-0008 lacked was the number, and the number is 4.67 in the variable the whole model exists to infer.

One of the seam's properties survives without the seam: a tree's update still names neither the field, nor the other tree, nor the error probability, and its walk still runs in a test with none of the three. Under the block it never reaches a leaf, so it never needs the link at all.

## The revert is exact, and that is the evidence

`git revert f13fef9` is the whole code change, and at that commit **`git diff a3d1ae6 -- src tests` is empty**. `TBlockUpdate.h` already addressed the hash-map storages by cell handle — those landed at `6d42857`, three commits before the deletion — so ADR-0009 needed nothing here and there was no port to do.

So the reinstated binary **is** the reference binary of #69's measurement, and that measurement's reference column is a prediction about this code rather than a historical note.

**No new chains were run**, and that is a decision rather than an omission. A confirming run would measure the same source with different seeds: it would reproduce the reference column, and it could only weaken the record by inviting a reader to ask why the third digit moved.

Three things after the revert do change the source, and none of them can change a chain: `field_update.h` and its tests are deleted, the `--K` message and the `--fix_Z` log line stop using a retired term, and that log line now says what the flag holds — every node above the leaves, with both tree fields still moving under the block.

## Considered options

**Run longer** — `findings.md` repair 1. Rejected. It is 5.9 times fewer effective field draws per second of wall clock, in the one variable the model exists to infer, and it is paid on top of a chain that is already long for another reason: the phylogenetic parameters had not converged at 20 000 iterations on a 128 x 256 fixture in either binary.

**Fuse the two leaf blocks and the field into a single pass** — repair 2. Rejected. It recovers the 1.25 in iteration time and leaves the 4.67 in mixing untouched. It is the right repair for the smaller quarter of the problem.

**A block Metropolis proposal beside the single-variable walk** — repair 3, and ADR-0008's own first choice. Rejected. Its benefit is argued and not measured, which is exactly the position ADR-0005 was in until issue #69 went and priced it. It also keeps three passes over the leaf-pair space and adds an accept/reject on top, so it cannot recover the iteration cost. The block, by contrast, exists, is tested, and its number is known. If the small-`omega` regime ever stops being the design point, this is still the option to reach for first.

**Keep the field's own pass beside the block, uncalled.** Rejected. The backend parity gate never runs it and no record argues for it, so a tested but unreachable update claims a support it does not have — and it would rot against the storage interface, whose handle lifetime convention is exactly the kind of rule a later change moves.

**`field_math::prob_field_cell_is_one` is kept, though the same deletion leaves it with no caller in `src/`.** The rot argument above is about a pass over storages, and this is arithmetic: it takes two probabilities and an error probability, it holds no storage and no handle, and five tests pin it. It is also the kernel any of repairs 1 to 3 would call if a larger `omega` ever reopens them. An unreachable traversal and an unreachable expression are not the same risk.

**Drop `f13fef9` from the branch instead of reverting it.** Rejected. The branch is pushed, so that is a force-push over a shared ref. Worse, `findings.md` and ADR-0008 both cite that sha, and dropping it would orphan the reference of the measurement this record stands on. The detour in the history is not noise; it is the evidence.

## Consequences

**ADR-0008 keeps its measurement and loses its decision.** Its banner names this record. A reader who arrives at that table for the price still gets the price; the verdict is here.

**ADR-0005's block paragraph is live again**, and its banner says so. Everything else that record reaches was never in question: the factorisation, the link table, the error probability's support, the six-counter collapse and the tree-field-versus-field distinction.

**ADR-0009 is untouched, and it is what made this a revert.** A storage is still addressed by cell handle, the sparse storages still hold their cells in a hash map, and both backends still traverse alike. Had the block been deleted before the window was, this record would have been a rewrite.

**The backend parity gate is unaffected.** Both backends still compute one chain, and the gate is still the only thing that asserts it.

**Issue #57 is closed with its ledger split.** Its storage half landed and stands (ADR-0009, user stories 12 and 24 to 32). Its leaf-layer half is reversed here (stories 1 to 6, 8 to 11, and 19 to 21). Issue #71 carries this record.

**The mixing harness's default pair is now source-identical.** `run.sh` builds its `current` binary from the working tree and defaults `ACOL_MIXING_REFERENCE` to `a3d1ae6`, which the working tree now matches, so a default run compares two identical binaries and measures nothing. A comparison from today's tree sets `ACOL_MIXING_REFERENCE=f13fef9`, and reproducing the recorded table means running from a tree at `f13fef9`. `model_validation/mixing_cost/README.md` records both. The script itself is left alone, so the command `findings.md` prints stays the command that was run.

**CONTEXT.md swaps two terms.** *Block update* is a live term again, and *Field update* is retired — it names what ADR-0008 built and this record retired, so the phrase can still be resolved where it survives in `findings.md` and in the commits of this branch. *Update*'s iteration order goes back to the block first.

**ADR-0006's and ADR-0007's mentions of the field update are left alone.** Both records carry their own context, ADR-0009 set the precedent of not sweeping such mentions, and a doc-wide rewrite is a large diff over a pointer that already works.
