This directory shows a couple of tests to validate the MRF model. See [here](https://anticipated-chemistry-of-life.github.io/acol-dws/notes/ufdb3a64xjinl5jopn58g8q/)

- In the [single root](./single_root) directory, the species tree has a single root and only leaves that are directly linked to the root.
- In the [balanced tree](./balanced_tree) directory, the species tree has a single root and always two children per node.
- In the [convergence_bigger_tree](./convergence_bigger_tree) directory, the species and the molecules tree are bit more complex. This directory was used to check if the model converges with a bigger tree and if it did, after how many iterations.

## Independent-field validation

`simulate_independent.py` samples the field from a Python reference implementation
that shares no code with the C++ binary, so a disagreement between them is a bug
in one of the two. The reference draws one pass down the species tree — root from
the stationary distribution, then every node given its parent — with no molecules
tree involved at all. The molecules dimension is then pinned neutral so the C++
model reduces to exactly that process; see
[ADR-0001](../docs/adr/0001-neutralise-molecules-dimension-for-validation.md).

```bash
uv run python simulate_independent.py --seed 42     # writes a scenario directory
cd independent_y_s255_m255_seed42

bash check_neutrality_invariant.sh                  # run this first

bash rung1_pin_field_and_states.sh                  # then each rung in order
uv run python validate_independent.py . rung1_pin_field_and_states
```

Run the rungs in order and stop at the first failure. Each pins strictly less
than the one before, so a failure localises the fault: rung 1 pins the field and
both trees' node states and is close to closed form; rung 2 adds the Z Gibbs
sweep; rung 3 infers everything from observations, against the simple error model
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
simulator; the rungs test inference.

### What the rungs cannot see

Pinning the molecules dimension neutral is what makes an independent reference
possible, and it also makes an entire class of faults invisible. `log_nu`
drifting downward on a doubly balanced tree is one of them.

`diagnose_normaliser.py` covers that case by shrinking the field until every
configuration can be enumerated, so the quantity the rungs cannot reach — the
normalising constant of the two-tree product — can be computed exactly rather
than estimated:

```bash
uv run python diagnose_normaliser.py
```

It reports where the C++'s objective peaks against where the correctly
normalised one does, for a molecules dimension swept from neutral to strongly
non-neutral, and then reproduces the drift as an MCMC and removes it. See
[ADR-0002](../docs/adr/0002-the-two-tree-product-is-unnormalised.md).

### Remarks

- **The branch-length budget is conserved.** Bins must sum to
  `n_branches * n_bins / 2`, enforced once at startup and preserved by every
  proposal thereafter, which only ever moves +1 on one branch and -1 on another.
  Branch lengths summing to anything else are _unreachable_, not merely unlikely.
- **Initial-value files are dispatched by filename.** A `name`/`value` file is
  only matched up by parameter name when the filename contains `trace`,
  `simulated`, `meanVar`, `statePosteriors` or `posteriorMode`. Otherwise it must
  be a bare one-column file of exactly the right length. Renaming
  `simulated_pinned_molecules.txt` breaks every run script.
- **Research effort is driven by log paper counts, not raw ones.** The counts
  are read raw (`read_paper_counts` in `src/lotus/paper_counts.cpp`) and the
  `log(count + 1)` is applied by `lotus_math::TReportingModel`. Simulating LOTUS
  data from raw counts is indistinguishable from an inference bug.
- **`--numThreads all` is not reproducible.** Two identical invocations under the
  same `--fixedSeed` give posterior means differing by ~0.2 posterior standard
  deviations. Any test that compares runs exactly must pass `--numThreads 1`;
  `check_neutrality_invariant.sh` does.
- **"internal nodes" means two different things.** The startup log line counts
  internal nodes _excluding_ roots, while `TPhylogeny::n_internal_nodes()`
  _includes_ them. Neither sizes the node state any more: since ADR-0005 that
  dimension spans every node of its tree, leaves included.
