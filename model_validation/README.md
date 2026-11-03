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
both trees' node states, leaves included, and is close to closed form; rung 2
frees the node states, so the block update draws both tree fields given the
pinned field and each tree walks the nodes above them; rung 3 infers everything
from observations, against the simple error model alone, then LOTUS alone, then
both. What `--Y.update false` and `--Z.update false` hold is ADR-0019's.

`validate_independent.py` reports by default; pass `--gates` to make it exit 1 if
a gated coverage lies outside its binomial band. The gate is absolute, not
relative to rung 1: rung 1 knows every node state, which no later rung can match.
Alpha is gated at every rung. Log ν is gated only at rung 1: once the node states
are inferred, one shared var log ν sets every clique's interval at once, and an
exact posterior falls outside the band on about half of all seeds (#112). Branch
lengths are never gated: once the node states are freed, a calibrated posterior
on seed 42 ranks them at a Spearman of about 0.1-0.2, against rung 1's 0.88 (#110).

`validate_replicates.py` is the log ν check for rungs 2 and 3. It scores one rung
in many scenarios, one per seed, and gates each tree's mean alpha and log ν
coverage in standard errors of the scenarios' own scatter, as `compare_fields.py`
does. A shortfall every scenario shares fails it; one unlucky seed does not. The
rung scripts read `ACOL_THREADS` (default `all`) and skip their build when
`ACOL_NO_BUILD` is set, so several can run side by side:

```bash
for s in $(seq 1 16); do uv run python simulate_independent.py --seed "$s"; done
pixi run build release-s
printf '%s\n' independent_y_s255_m255_seed{1..16} |
    xargs -P 16 -I{} env ACOL_NO_BUILD=1 ACOL_THREADS=4 bash {}/rung2_pin_field.sh
uv run python validate_replicates.py rung2_pin_field \
    independent_y_s255_m255_seed{1..16} --gate 3
```

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

Under [ADR-0005](../docs/adr/0005-each-tree-owns-its-leaf-level-field.md) the
joint sums to one by construction. A chain cannot check that. The sum runs over
every `(Z_s, Z_m, Y)`. So `src/independent/enumeration.py` shrinks both trees
until every configuration can be listed, and the suite asserts the total:

```bash
uv run pytest tests/test_independent.py
```

Both trees stay active, and a test pins that they are. The same enumeration used
to measure how far
[ADR-0002](../docs/adr/0002-the-two-tree-product-is-unnormalised.md)'s constant
moved with the parameters, and it ran by hand. It runs in the suite now.

Read the total for what it is. It reads the two node-state densities, and it
reads whether a variable is scored twice, which is ADR-0002's defect. It is blind
to the shape of the link, to the error probability and to the orientation of a
leaf block. The field rate and the reference draw pin those. The module header
says so at more length.

### The LOTUS+MS scenarios

```bash
uv run python simulate_independent.py --mass_spec dense     # independent_y_s255_m255_seed42_ms_dense
uv run python simulate_independent.py --mass_spec sparse    # independent_y_s255_m255_seed42_ms_sparse
```

Each is the scenario above plus simulated mass-spec data, drawn by
`src/independent/mass_spec.py` after everything else, so adding `--mass_spec`
leaves the rest of a seed's scenario unchanged. The story is the model's own
([ADR-0018](../docs/adr/0018-mass-spec-candidates-enter-by-likelihood-ratio.md)),
so a recovery failure points at the code and not at misspecification.

- **Layouts.** `dense`: two runs per species leaf. `sparse`: eight species, chosen
  at random, with one run each. Both use 255-node trees and one filter.
- **Detection.** A present molecule is detected in a run with its filter
  probability, an absent one with κ = 0.01. Filter probabilities are drawn uniformly
  on the byte grid, `i / 255` for `i` in 0..254. Each detected molecule produces
  exactly one feature.
- **Unknown features.** Each run also gets Poisson(λ) unknown features. λ is the
  expected number of molecule features per run under the true Y, filter
  probabilities and κ, so about half of all features are unknown. ρ = λ / the
  number of molecule leaves.
- **Candidates.** n_f ~ Uniform{1..20} for both kinds of feature. A molecule
  feature lists its true molecule plus other molecule leaves; an unknown feature
  lists molecule leaves only. Both draw uniformly without replacement, so every
  proposed structure is a molecule leaf.
- **Scores.** A wrong row scores x ~ Uniform(0,1), a true row x ~ Beta(5,1). The
  log correct density is ln 5 + 4 ln x and the log incorrect density is 0.

Truth files, all under the `simulated_` prefix that picks the stattools reader:

| File | Holds |
| --- | --- |
| `simulated_ms_data.tsv` | The data, in the loader's six columns. The only mass-spec file an inference run reads. |
| `simulated_ms_assignment.tsv` | The true assignment, in the layout of the C++ assignment posterior: `species extract feature molecule fraction`, one row per candidate and one for the unknown molecule (`<unknown>`), fraction 1 on the true explanation and 0 elsewhere. |
| `simulated_parameters.txt` | Also `filter_proba_<molecule leaf>` (the byte index the C++ stores; the probability is the index over 255), `contamination_proba` (κ), `ms_lambda` (λ) and `ms_unknown_rate` (ρ). |

The mass-spec truth is not added to `truth_<tree>.txt`, which `validate_independent.py`
reads its tree scores from, and `meta.json` gains a `mass_spec` entry with the counts.

The rungs are the same three, on the `lm` build, named
`rung1_pin_field_and_states`, `rung2_pin_field` and `rung3_from_data_lm`. Every one
passes `--ms_data` and `--ms_unknown_rate` at its true value, plus
`--write_assignment_posterior`. Assignments, filter probabilities and κ are pinned
in none of them, and ρ is the only mass-spec truth a run is told.
`validate_independent.py` still scores alpha, log nu and the branch lengths on them.

#### Mass-spec scores

On a LOTUS+MS scenario `validate_independent.py` adds a `-- mass spec --` block,
scored by `src/independent/recovery.py` from `simulated_Y.txt`, `simulated_ms_assignment.tsv`,
`simulated_parameters.txt` and the rung's `acol_Y_posterior.tsv`, `acol_assignment_posterior.txt` and `acol_trace.txt`:

- **Field.** AUC of P(Y=1) against the true Y, and a reliability table (ten equal bins:
  mean predicted, observed frequency, count), each for species with runs and species
  without. A cell the posterior file omits was never a one and counts as 0. Reported,
  no threshold.
- **Assignment.** From `acol_assignment_posterior.txt` against
  `simulated_ms_assignment.tsv`, over every (run, feature): the accuracy of the
  posterior mode against the true assignment (a candidate molecule, or the unknown
  molecule); the mean posterior mass on the true assignment; and the mean posterior
  on the unknown molecule, reported separately for truly unknown features and for
  features of a detected molecule (a good fit is high on the first, low on the
  second). Reported, no threshold. Before scoring, the file's format is checked
  against `simulated_ms_data.tsv`: every (run, feature) is present with all its
  candidates and `<unknown>`, nothing else is, every fraction lies in [0, 1] and
  each feature sums to 1. A malformed or incomplete file fails the validation
  with a message naming the feature.
- **Filter probabilities.** 90% credible-interval coverage over the molecule leaves,
  stratified by the number of runs the molecule is present in (the runs of the species
  whose true Y carries it): `0`, `1`, `2`, `3-9`, `10+`; empty strata are skipped.
- **κ.** Whether its 90% interval holds the truth.
- **Pass rule.** A filter-probability stratum passes when its coverage lies inside the
  95% binomial band around 0.9, the central 95% of Binomial(n, 0.9) over n. The band
  is two-sided, so over-wide intervals fail as well as over-confident ones. κ is a
  single interval, and the band at n = 1 is [0, 1], which a miss would also pass; κ
  therefore passes when its interval holds the truth.

Pass and fail are printed and written to `validation_summary.json` under `mass_spec`;
they do not set the exit code. Scenarios without mass spec are scored exactly as before.

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
  generator -- and `pixi run parity` gates it.
- **"internal nodes" means two different things.** The startup log line counts
  internal nodes _excluding_ roots, while `TPhylogeny::n_internal_nodes()`
  _includes_ them. Neither sizes the node state any more: since ADR-0005 that
  dimension spans every node of its tree, leaves included.

## The mixing cost of dropping the block update

ADR-0005 argued for the leaf layer's eight-state block update on the grounds
that the field and both tree fields are metastable under single-variable draws at
a small error probability, and that the failure would present as slow mixing
rather than as a bug. Issue #68 deleted that block update. `mixing_cost/`
measures what that costs, since slow mixing is not something a test suite trips
over.

```bash
bash mixing_cost/run.sh
```

Two binaries — the revision before the block update was deleted, and the working
tree — infer from one data set the reference binary simulated. The data set is
fixed and the seed is not: the pattern of random-number consumption changed, so
the two cannot be compared seed for seed.
`compare_mixing.py` reads the joint density trace's autocorrelation time and both
tree field posteriors, the second as a point on the ADR-0005 ridge.
[`mixing_cost/README.md`](./mixing_cost/README.md) says how, and
[`mixing_cost/findings.md`](./mixing_cost/findings.md) is the recorded result.
