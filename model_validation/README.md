This directory shows a couple of tests to validate the MRF model. See [here](https://anticipated-chemistry-of-life.github.io/acol-dws/notes/ufdb3a64xjinl5jopn58g8q/)

- In the [single root](./single_root) directory, the species tree has a single root and only leaves that are directly linked to the root.
- In the [balanced tree](./balanced_tree) directory, the species tree has a single root and always two children per node.
- In the [convergence_bigger_tree](./convergence_bigger_tree) directory, the species and the molecules tree are bit more complex. This directory was used to check if the model converges with a bigger tree and if it did, after how many iterations.

## Independent-field validation

`simulate_independent.py` draws a whole scenario from a Python reference
implementation that shares no code with the C++ binary, so a disagreement between
them is a bug in one of the two. **Both trees are active.** Each draws its own
node state in one pass down its nodes — root from the stationary distribution,
then every node given its parent, leaves included — and the field is a noisy AND
of the two leaf blocks. Nothing is neutralised: under
[ADR-0005](../docs/adr/0005-each-tree-owns-its-leaf-level-field.md) the reference
is exact with both trees running, which is precisely what
[ADR-0001](../docs/adr/0001-neutralise-molecules-dimension-for-validation.md)
could not do.

```bash
uv run python simulate_independent.py --seed 42     # writes a scenario directory
cd independent_y_s255_m255_seed42

bash rung1_pin_field_and_states.sh                  # each rung in order
uv run python validate_independent.py . rung1_pin_field_and_states
```

The neutralised rung ADR-0005 keeps as a cheap regression check is **not** here.
It needs a neutralised *scenario*, because a neutralised inference against a
scenario both trees drew is a misspecified fit rather than a regression check.
Issue #44 builds it, along with the rest of the ladder.

Run the rungs in order and stop at the first failure. Each pins strictly less
than the one before, so a failure localises the fault: rung 1 pins the field and
both trees' node states and is close to closed form; rung 2 adds the Z Gibbs
update; rung 3 infers everything from observations, against the simple error model
alone, then LOTUS alone, then both.

Rung 1 is also the empirical ceiling. Its scores are what the looser rungs should
be judged against, which is why `validate_independent.py` reports rather than
gates by default — pass `--gates rung1_.../validation_summary.json` to score a
later rung against it.

`visualize.py` plots a rung's inference output against the truth. Its filename
options default to what the C++ binary writes for a scenario simulated by that
same binary; an independent scenario names its simulation output differently, so
it passes them. Run it from `model_validation/`:

```bash
uv run python visualize.py independent_y_s255_m255_seed42/rung3_from_data_ls \
    --true-values simulated_parameters.txt \
    --true-branch-lengths simulated_parameters.txt \
    --true-y simulated_Y.txt
```

`simulated_parameters.txt` covers both trees, which is why it also serves as the
pooled `--true-branch-lengths` file — a value without a `{tree}` placeholder is
shared by every tree. It carries no `gamma` or `epsilon`, so those get no panel;
`--true-scalar name=value` supplies a truth that no file holds, where one applies.

Separately, `replicates.sh` runs the C++ simulator under the same parameters and
`compare_fields.py` compares both against the analytic prediction. That tests the
simulators; the rungs test inference. Five statistics: the field's density against
the product of the two adjusted rates, each tree field's per-clique density
against that clique's `alpha`, the rate at which the field reads 1 in each link
bucket against `P_k`, the sibling calibration of the species tree field, and the
six counters the C++ traced against a tally of the same cells.

`--gate N` turns the report into a check: each deviation is scored against the
scatter of the replicates themselves, and the command exits non-zero above `N`
standard errors. The scatter comes from the replicates and not from a binomial
formula, because the cells of one replicate are correlated along the trees.

### The enumerable case

A chain cannot reach the normalising constant of the whole model, so
`diagnose_normaliser.py` shrinks the field until every configuration can be
enumerated and computes it exactly:

```bash
uv run python diagnose_normaliser.py
```

It reports where the C++'s objective peaks against where the correctly
normalised one does, for a molecules dimension swept from neutral to strongly
non-neutral, and then reproduces the drift as an MCMC and removes it. See
[ADR-0002](../docs/adr/0002-the-two-tree-product-is-unnormalised.md). Repurposing
it to assert that the new joint sums to one is issue #43.

### Remarks

- **The branch-length budget is conserved.** Bins must sum to
  `n_branches * n_bins / 2`, enforced once at startup and preserved by every
  proposal thereafter, which only ever moves +1 on one branch and -1 on another.
  Branch lengths summing to anything else are _unreachable_, not merely unlikely.
- **Initial-value files are dispatched by filename.** A `name`/`value` file is
  only matched up by parameter name when the filename contains `trace`,
  `simulated`, `meanVar`, `statePosteriors` or `posteriorMode`. Otherwise it must
  be a bare one-column file of exactly the right length. Renaming
  `simulated_parameters.txt` breaks `replicates.sh`.
- **Research effort is driven by log paper counts, not raw ones.** The counts
  are read raw (`read_paper_counts` in `src/lotus/paper_counts.cpp`) and the
  `log(count + 1)` is applied by `lotus_math::TReportingModel`. Simulating LOTUS
  data from raw counts is indistinguishable from an inference bug.
- **`--numThreads all` is not reproducible for `infer`.** Two identical
  invocations under the same `--fixedSeed` give posterior means differing by ~0.2
  posterior standard deviations. The cell draws are hashed from their position
  (ADR-0007); what is left is the alpha and nu moves, which still draw from the
  thread-local generator inside the clique loop. Any test that compares `infer`
  runs exactly must pass `--numThreads 1`. A `simulate` run is reproducible at any
  thread count -- it is a forward draw and takes nothing from the thread-local
  generator -- and `just parity` gates it.
- **"internal nodes" means two different things.** The startup log line counts
  internal nodes _excluding_ roots, while `TPhylogeny::n_internal_nodes()`
  _includes_ them. Neither sizes the node state any more: since ADR-0005 that
  dimension spans every node of its tree, leaves included.
