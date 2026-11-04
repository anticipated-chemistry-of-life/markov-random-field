# Working log: diagnosing and fixing parameter recovery in the MRF model

Structured record of the investigation into why simulated parameters were not being
recovered, and of the changes made as a result. Newest work is appended at the end.

Started 2026-08-14. All numbers below were verified by writing results to a file and
re-reading them (see "Reliability note").

---

## Reliability note (important)

Partway through this work it emerged that command output in the working environment was
returning **altered numbers**: re-reading a byte-identical file (md5 verified) through the
shell produced different values on different invocations — e.g. `gamma_species` read as
48.9 on one invocation and 2471.7 on another, from the same unchanged
`acol_meanVar.txt`.

Every figure quoted before that discovery is suspect. From that point on, all measurements
were written to a file and read back through a separate path. Figures in this document are
from the verified path only. The scripts live in the session scratchpad
(`compare.py`, `verify_all.py`, `convergence.py`, `exact_ll.py`, `mle_bias.py`).

Additionally, the validation scenario directory was regenerated and re-simulated twice
during the session by an external process, and several diagnostic run directories were
deleted. Numbers are only comparable **within** a block below, never across blocks that
used different simulated data.

---

## 1. Original question

> Every time I run simulations and then try to infer the simulated parameters of my data,
> the parameters of the species seem to have a hard time being inferred. Find out if there
> might be a bug or it is normal.

### 1.1 It is not species-specific

Measured on the then-current `s_balanced_255_m_balanced_255` run:

| quantity | species | molecules |
|---|---|---|
| alpha, corr(true, posterior mean) | 0.66 | 0.64 |
| log_nu, corr(true, posterior mean) | 0.11 | −0.01 |

A mirror-image test (make one tree carry the signal and the other neutral, then swap)
gave identical recovery either way, ruling out an asymmetry between dimension 0 (species,
uses `TSheet`) and the last dimension (molecules, uses `_clique_last_dim`).

What made species *look* worse in the plots: `species_var_log_nu` happened to collapse to
0.0115 while `molecules_var_log_nu` landed at 0.22, so every species `log_nu` was shrunk
onto the hyper-mean and plotted as a vertical line.

### 1.2 Bug found: variance passed where numpy wants a standard deviation

`model_validation/src/params.py` drew

```python
np.random.normal(self._mean_log_nu_scalar, self._var_log_nu_scalar, size=n)
```

numpy's third positional argument is `scale` = **standard deviation**. The C++ side treats
`var_log_nu` as a **variance** — `stattools/core/stattools/Priors/TNormalInferred.h:35` does
`_sd = sqrt(_var->value())`.

Verified on the simulated data of the time: recorded `species_var_log_nu = 0.0865`, actual
variance of the drawn `log_nu` values = 0.0084 (sd 0.0919 ≈ the recorded "variance").

Consequences: (a) the `var_log_nu` panel compared a variance against a standard deviation
and could never look right; (b) the simulated `log_nu` spread was ~3.4× narrower than
intended.

**Fixed** — `np.sqrt(self._var_log_nu_scalar)`. After the fix: recorded 0.0929 vs drawn
variance 0.0942.

### 1.3 Not a bug: per-clique `log_nu` is not estimable at that spread

Two independent checks:

* Exact profile likelihood of `log_nu` per clique, computed from the true simulated Y and Z
  (`exact_ll.py`).
* Parametric bootstrap on a pure single-tree CTMC with the same topology and branch lengths
  (`mle_bias.py`).

Both give a **standard error of ≈ 0.53** for `log_nu` from 254 branches, against a true
spread of **sd 0.09**. Signal-to-noise ≈ 0.17, so the maximum attainable correlation between
truth and estimate is ≈ 0.17. A flat scatter plot is the expected outcome, with no bugs at
all. Finite-sample MLE bias was only −0.044, so that is not the explanation either.

`alpha` has a genuine spread (Beta(2,2), sd 0.22) and is recovered well (corr 0.80+ with the
field fixed).

### 1.4 Bug found: the nu/alpha update ignores the field's normalising constant

The Markov field is a product of the two per-tree factors sharing the leaf layer Y:

    p(Y, Z_s, Z_m | theta) = f_s(Y, Z_s) * f_m(Y, Z_m) / C(theta)

`TTree::_update_nu_or_alpha` scored proposals with `f_s` alone, dropping `C(theta)`, which
depends on nu, alpha and the branch lengths of both trees.

Demonstrated with Y and Z **pinned at the simulated truth**, so the latent field is perfect:

| setup | realized frac(Y=1) | true alpha | inferred alpha | true log_nu | inferred log_nu |
|---|---|---|---|---|---|
| both trees active | 0.867 | 0.600 | 0.850 | −0.500 | −0.214 |
| molecules made neutral | 0.602 | 0.600 | **0.596** | −0.500 | **−0.512** |

"Neutral" = alpha 0.5 and log_nu 3.4 (above the `nu > 25` threshold where the code uses the
stationary matrices), which makes the molecules factor constant and `C` independent of
theta. The bias vanishes entirely.

Ruled out as explanations: simulation burn-in (30000 instead of 2000 iterations gives an
identical +0.285 vs +0.286 bias) and finite-sample MLE bias (−0.044).

### 1.5 Not a bug: gamma was simulated in a flat region

Research effort is `prod_i (1 - exp(-gamma_i * log(count_i + 1)))`. With the scenario's paper
counts the factor is already 0.997–1.0 at gamma = 5, so every gamma above that fits
identically. Under a uniform prior `gamma_species` drifted to 74 ± 39.

Also noted: `simulate` draws any parameter not fixed from a file **from its prior**, not from
the CLI default (`TUniformFixed::_simulateUnderPrior` draws U(0,1)) — which is why a later
scenario recorded a true gamma of 0.6 despite `--gamma` defaulting to 5.

### 1.6 Minor: leaves with zero papers are silently dropped

`model_validation/src/tree.py` filters `number_of_papers > 0`. `TTree::get_paper_counts`
defaults missing leaves to 0.0, giving factor `1 - exp(0) = 0`, so `P(L=1 | Y=1) = 0` for
every LOTUS cell in that row/column. In that scenario 2 species and 3 molecule leaves were
dropped, silently blanking ~640 of 16384 LOTUS cells.

### 1.7 Harness gaps

* Branch lengths are never validated: being a discrete type they land in
  `acol_species_posteriorMode.txt`, not `acol_meanVar.txt`, so `visualize.py`'s inner join
  drops the whole panel.
* `visualize.py` falls back to `acol_input_simulated.txt` for true branch lengths, but the
  actual simulated bins are re-normalised and written to `acol_<tree>_simulated.txt`.

---

## 2. The sqrt fix, and what it exposed

After fixing `params.py` and regenerating + re-simulating + re-inferring:

**With Y and Z pinned at truth, recovery is good** — this is the estimator working as
intended once the field is known:

| | species | molecules | true |
|---|---|---|---|
| alpha corr | 0.849 | 0.894 | — |
| log_nu corr | 0.723 | 0.728 | — |
| mean_log_nu | −0.306 | −0.668 | −0.272 / −0.606 |
| gamma | 0.601 | 0.648 | 0.600 |
| epsilon / eps_simple | 0.0013 / 0.198 | | 0.001 / 0.200 |

**With Y and Z inferred, the run collapsed completely**: joint field density ran from
−55933 to −154 (a frozen field), the Y posterior was empty (Y identically zero),
`species_mean_log_nu` → −18.1, `gamma_species` → 2471.7, `epsilon` → 0.141 (true 0.001),
`epsilon_simple_model` → 0.438 (true 0.200).

The pre-fix run showed the same pathology, milder — the joint density also climbed
monotonically to the last iteration and never plateaued. The wider `log_nu` spread after the
fix simply let it run further.

**Mechanism.** Freezing the field is worth roughly 56,000 log units of `log f` (−56139 →
−159). The `-log C` that would grow equally is missing, and the error rates can absorb the
resulting data mismatch.

---

## 3. Prior swaps (implemented)

`src/Types.h`: `PriorOnGamma` → `TGammaFixed`, `PriorOnErrorRate` → `TBetaFixed`,
`PriorOnMeanLogNu` → `TNormalFixed`. `src/cli.h`: `--prior_gamma` (default `2,4.6`),
`--prior_epsilon` (`1,99`), `--prior_mean_log_nu` (`0,1`). `src/TCore.cpp`: wired into the
`TParameterDefinition`, which is the only place stattools reads fixed prior parameters
(`TNodeTyped` forwards them to `setFixedPriorParameters` at construction; there is no
generic `--<name>.priorParameters` option).

Verified result:

| | uniform | new priors | true |
|---|---|---|---|
| `gamma_species` | 2471.7 | **0.445** | 0.6 |
| `gamma_molecules` | 2547.4 | **0.427** | 0.6 |
| `epsilon` | 0.1411 | 0.1405 | 0.001 |
| `species_mean_log_nu` | −18.14 | −7.54 | −0.272 |
| Y posterior non-zero cells | 0 | 7 | — |

**The gamma prior works. The epsilon prior is inert. The mean_log_nu prior halves the
runaway but does not stop it.** The arithmetic says why: freezing buys ~56000 log units;
Normal(0,1) at log_nu = −7.5 costs 28, Beta(1,99) at eps = 0.14 costs ~15. Three orders of
magnitude short. Only a hard constraint bites — pinning epsilon stopped the collapse
outright (density plateaued at −39907, Y survived, alpha corr 0.72).

---

## 4. Phase 1: pseudo-likelihood on internal nodes (implemented)

See `model.tex` for the mathematics.

Changes:

* `src/TClique.h/.cpp` — `_calculate_log_prob_node_to_children` made public as
  `calculate_log_prob_node_to_children`, with a runtime `bool UseTryMatrix` (runtime rather
  than a template because the body needs a complete `TTree`, which `TClique.h` cannot
  include).
* `src/tree/TTree.h` — `_compute_LL_old_and_new_nu_or_alpha` replaced by
  `_add_pseudo_likelihood_term`, which builds `sum_log[0]`/`sum_log[1]` exactly as
  `update_Z` does and adds the **normalised** conditional via a logistic of the log
  difference (clamped at 700 so a vanishing conditional gives ~1e-305 rather than `-inf`).
* The update loop iterates `_internal_nodes` instead of all nodes.

115/115 unit tests pass. (The test binary needs `DYLD_LIBRARY_PATH` pointing at libomp — it
has no `LC_RPATH`, which is pre-existing.)

Verified, same data and config, before vs after:

| | pre-Phase-1 | post-Phase-1 | true |
|---|---|---|---|
| `species_mean_log_nu` | −7.386 | **−1.002** | −0.446 |
| `molecules_mean_log_nu` | −7.414 | −1.054 | −0.466 |
| species log_nu bias | −6.94 | **−0.55** | — |
| `epsilon` | 0.138 | 0.058 | 0.001 |
| `epsilon_simple_model` | 0.443 | 0.320 | 0.200 |
| Y posterior cells | 9 | **16382** | — |
| Y accuracy | 0.596 (all-zero) | **0.759** | — |
| joint field density | −55933 → −154 (frozen) | −56401 → −41593 (fluctuating) | — |

**The field no longer collapses with the error rates left free.** With Y and Z pinned at
truth the estimator is essentially unbiased: log_nu bias +0.055 (species) / +0.024
(molecules), against ±0.25–0.30 for the raw-factor version.

Remaining: `mean_log_nu` −1.00 vs true −0.446; `epsilon` 0.058 vs 0.001; Y accuracy 0.759.

---

## 5. Convergence diagnostics

Question asked: the joint field density now fluctuates — does that just mean running longer?

Fluctuation is the correct shape (the log-density of an MCMC state is a random variable, not
something being maximised). The pathology that was fixed was *monotone climbing*.

But from the completed 20000-iteration run:

| quantity | ESS (2nd half) | iterations per independent draw |
|---|---|---|
| joint field density | 10 | ~1000 |
| `species_mean_log_nu` | 11 | ~900 |
| `epsilon` | 12 | ~830 |
| `epsilon_simple_model` | 14 | ~740 |
| `gamma_species` | 58 | ~170 |

No evidence of residual trend (half-to-half shift 0.8 sd; slope flips sign between the whole
chain and the second half). Reaching ESS ≈ 200 would need 10–20× more iterations.

Longer runs will **not** fix the residual bias: two independent runs land on
`species_mean_log_nu` ≈ −1.0/−1.1 against a true −0.446. That is a property of the target,
not a convergence failure.

(A NaN spotted in one trace was a half-written final line of a still-running inference, not
a numerical problem; the completed run has zero NaN/inf across all 520 columns.)

---

## 6. Where the remaining bias comes from

Decisive experiment — `species_mean_log_nu`, true −0.446:

| regime | inferred | bias |
|---|---|---|
| Y **and** Z pinned at truth | −0.402 | +0.055 |
| **Y pinned at truth, Z free** | **−0.854** | **−0.400** |
| Y and Z both free | −1.002 | −0.546 |

With the true Y and only Z sampled, ~73% of the total residual bias is already present. After
Phase 1 theta's only information is the latent Z layer; Z is sampled given theta; a slightly
too-small nu yields an over-smoothed Z, which supports a small nu. Self-reinforcing, and it
is a bias in the target — better mixing will not remove it.

Phase 2 (leaf conditionals) breaks the loop by tying theta to Y, which is anchored by data.

---

## 7. Scaling analysis (towards thousands of species and molecules)

Per iteration, with S species and M molecule leaves:

| | cost |
|---|---|
| Y sweep | O(S·M) |
| Z updates (both trees) | O(S·M) |
| **theta update** | M cliques × S nodes × 2 params × 2 (old/new) × ~3 (parent+children) × 2 trees ≈ **24·S·M** |

The theta update dominates by an order of magnitude, driven by having one (alpha, nu) per
clique. Per-clique nu was measured to be non-estimable (SE ≈ 0.5 vs true spread ≈ 0.3), so at
5000 cliques that is ~10000 unidentifiable parameters per tree.

Structural observation: **Y cells are conditionally independent given Z_s and Z_m**
(`_calculate_log_prob_field` reads only the two parent states; the data terms are per-cell),
so the Y sweep is already an exact block Gibbs draw. The sampler is a three-block Gibbs over
{Y}, {Z_s}, {Z_m} with `Z_s ⊥ Z_m | Y`. The mixing bottleneck is entirely the single-site
Gibbs on Z.

Agreed plan:

1. Exact tree block-sampling of Z (Felsenstein pruning + backward sampling) — biggest mixing
   win, exact.
2. Phase 2 leaf terms — removes ~3/4 of the remaining bias. Independent of (1).
3. Revisit per-clique alpha/nu.
4. Micro-optimisation last.

Phase 2 difficulty was revised down from "several days" to ~2 days: the only cross-tree
quantity needed per Y cell (i,j) is the bit `Z_m[i, parent_m(j)]`, and transposing those bits
once per iteration into a bitset costs ~3 MB at 5000×5000 and one O(S·M) pass.

---

## 8. Z block sampler — implemented, but it did not deliver

### 8.1 What was built

* `src/tree/io/load_from_file.cpp` — `_internal_nodes_post_order`, built once at load from a
  reversed breadth-first walk of the roots (children before parents), with a guard that
  throws if the roots do not reach every node.
* `src/TClique.h/.cpp` — `update_Z` dispatches to `_update_Z_block` (default) or
  `_update_Z_single_site` (the old body, unchanged). The block sampler is Felsenstein
  pruning in logs, renormalised at every node by its maximum, with leaf children folded in
  directly as `P(t_u)_{s, Y_u}`; then a downward pass in reverse post-order samples the root
  and then each node given its already-sampled parent.
* `src/cli.h` — `--single_site_Z` selects the old sweep, so the two can be compared.

Builds clean, 115/115 unit tests pass. Wall clock for 30000 iterations (10×1000 burn-in +
20000): block 467 s, single-site 471 s — no measurable cost.

### 8.2 The A/B result

Same data, same seed 777, same everything but the Z sampler.

| | block (new) | single-site (old) | true |
|---|---|---|---|
| joint field density, ESS | **9** | **9** | — |
| iterations per independent draw | 1710 | 1676 | — |
| `species_mean_log_nu` ESS | 7 | 11 | — |
| `gamma_species` ESS | 28 | 142 | — |
| joint field density, 2nd-half mean | −27489 (sd 7210) | −36660 (sd 4344) | — |
| `species_mean_log_nu` | −1.569 | −1.186 | −0.446 |
| `epsilon` | 0.113 | 0.085 | 0.001 |
| Y accuracy | 0.804 | 0.849 | — |

**Two problems.**

1. **Mixing did not improve.** ESS is 9 either way. The hypothesis that single-site Gibbs on
   Z was the bottleneck was wrong.
2. **The two samplers disagree on the posterior**, which they should not: both leave
   `p(Z | Y, theta)` invariant, so the composite chain has the same stationary distribution.
   The field-density difference is ~9000 log units, roughly 3 standard errors given
   ESS ≈ 9 — suggestive but not conclusive at that effective sample size.

### 8.3 Two readings, not yet distinguished

* **Bug.** The block sampler targets something slightly different. Correctness is
  **not** verified: there is no test. `tests/TUpdate_Z_Tests.cpp` is an empty stub.
* **Working as intended.** The block sampler reaches higher-density field configurations
  that the single-site sweep was too stuck to find (field density −27489 vs −36660, and a
  wider explored range, sd 7210 vs 4344). Higher field density corresponds to a smaller nu
  and a more frozen field, i.e. further along the direction the pseudo-likelihood bias
  favours — which would explain the *worse* `mean_log_nu` and Y accuracy as the chain simply
  getting closer to a target that is still biased.

Both readings fit the data. Distinguishing them needs a correctness test: enumerate all
`2^n` configurations of `Z` for a small clique with fixed `Y`, `alpha`, `nu` and branch
lengths, compute `p(Z | Y)` exactly, and compare with the empirical frequencies from the
block sampler.

### 8.3b The enumeration test — verdict: the sampler is correct

To make the algorithm testable it was lifted out of the storage plumbing into a templated
free function `sample_clique_states_exact` in `src/TClique.h`. The tree is supplied as a set
of callables, so the production call site binds straight to `TTree`/`TCurrentState` with no
extra allocation, while the tests drive it from plain vectors — no `TTree`, no stattools
parameter graph.

`tests/TUpdate_Z_Tests.cpp` (previously an empty stub) now enumerates all `2^(#internal)`
configurations, computes `p(Z | Y)` in closed form, draws 200000 times, and compares the two
distributions on a Monte Carlo z-score scale. Five cases: a balanced tree, a high switching
rate, a low switching rate, a chain, and a root-with-one-leaf case where the root marginal
must equal `alpha` exactly.

**All five pass.** Full suite: 119 tests, all green.

The tests were then checked for power by mutation:

| mutation | caught |
|---|---|
| transpose the parent→child lookup in the downward pass | 4 of 5 fail |
| drop the pruning recursion (ignore child partial likelihoods) | 3 of 5 fail |

The one test that survives the transpose is the root-with-one-leaf case, which has no
non-root internal node and so legitimately never exercises that code path. An earlier version
of the low-switching-rate test also survived, because it used `alpha = 0.5`, at which the
transition matrix is symmetric and a transposed lookup is a no-op; it now uses `alpha = 0.35`.

**Conclusion: the block sampler is not buggy.** Reading two in §8.3 is the correct one — the
block sampler reaches higher-density field configurations that the single-site sweep could
not, which moves the chain closer to a target that is still biased by the Phase-1
pseudo-likelihood. The apparently worse `mean_log_nu` is the bias becoming *more* visible,
not a regression.

### 8.4 Why blocking Z alone was the wrong target

Y is already an exact block, and Z is now an exact block, yet mixing is unchanged. That
points at the coupling *between* the blocks rather than within either. Two-block Gibbs mixes
slowly exactly when the blocks are strongly dependent, and here Y and Z are nearly
deterministic functions of one another.

The fix is a bigger block: sample `(Y, Z)` **jointly, one clique at a time**. This is exact
and is a small extension of the code just written. For species clique `j`,

    p(Y_.j, Zs_.j | Zm, data)
        ∝ [ pi(root) * prod_v P_s(t_v)_{x_pa(v), x_v} ]          <- species tree factor
          * prod_i [ P_m(Y_ij | Zm[i, pa(j)]) * P(data_ij | Y_ij) ]   <- per-leaf emissions

which is precisely a two-state tree HMM with emissions at the leaves, so the same pruning
and backward-sampling recursion applies with the leaves sampled as well. Given `Zm` the
columns are mutually independent, so all species cliques can be done in parallel; then the
same for the molecule cliques given `Zs`. That attacks the Y-Z dependence directly.

---

## 9. Joint (Y, Z) clique block — step 1 of 2 done

### 9.1 The recursion, generalised and verified

`sample_clique_states_with_emissions` in `src/TClique.h`. Where the Z-only sampler treats leaf
states as observed, this gives every leaf an emission likelihood
`e_v(s) = P(everything attached to v | state of v is s)` and samples the leaves too. The
recursion becomes uniform across node types, with `e_v == 1` where nothing is emitted:

    L_v(s) = e_v(s) * prod_{u in children(v)} sum_{s'} P(t_u)_{s,s'} L_u(s')

so a leaf, having no children, simply has `L_v = e_v`. The downward pass then draws every node,
a leaf from `P(t_v)_{state(parent(v)), s} e_v(s)`.

For the real model the emission is exactly

    e_i(s) = P_molecules(t_j)_{Zm[i, parent(j)], s} * P(data_ij | s)

i.e. the *other* tree's parent term times the per-cell data likelihood. Both compiled data
sources factorise over cells, which is what makes this an emission rather than something that
couples the leaves to each other.

Four new tests in `tests/TUpdate_Z_Tests.cpp`, all enumerating the **joint** `2^(#nodes)`
configuration space (not just the internal nodes):

* balanced tree with pseudo-random leaf emissions
* low switching rate, where the tree and the emissions pull against each other
* a chain
* a bridge test: with emissions sharp enough to pin the leaves, the joint sampler must
  reproduce the Z-only sampler's posterior after marginalising the leaves away

All pass. Full suite: **123 tests**, all green. Mutation check: dropping the emission from the
pruning pass (keeping it only in the downward pass) fails **all four** joint tests.

### 9.2 What remains: wiring it in

The recursion is done and verified; making it the production sampler is the larger half and is
*not* yet done. It needs:

1. **A post-order over all nodes**, not just internal ones. Currently
   `_internal_nodes_post_order` filters the leaves out; the joint block needs them in.
2. **The other tree's parent-state bit per cell.** For species clique `j` and species leaf `i`
   that is `Zm[i, parent_m(j)]` — a strided walk through `Zm`, the same access problem Phase 2
   has. Same remedy: transpose those bits once per iteration into a bit-set, ~3 MB at
   5000x5000, one `O(S*M)` pass.
3. **Per-cell data likelihoods for a column of Y.** `TLotus::calculate_LL_update_Y` and
   `TSimpleErrorModel::probabilities_for_Y_update` already return exactly the two numbers
   needed, but they are currently fed from `_tmp_state_along_last_dim`, which is filled along
   the *last* dimension. The species pass walks a column, i.e. across that grain.
4. **Restructuring `TMarkovField`.** Today the sweep is `_update_all_Y` (both trees, cell by
   cell) then `_update_all_Z` (per tree). The joint block replaces both with: for each species
   clique draw `(Y_.j, Zs_.j)` given `Zm`; then for each molecule clique draw `(Y_i., Zm_i.)`
   given `Zs`.

Item 3 is the awkward one and is shared with Phase 2, so the two should probably be built
together.

---

## 10. Phase 2 implemented — and it does not yet work on a coupled field

### 10.1 What was built

* `TTree` keeps `_all_trees` (set in `initialize_cliques_and_Z`, which already received them).
* `_other_parent_state`: the other tree's parent state for every Y cell, transposed once per
  sweep into `[clique * n_leaves + leaf]` so the strided walk through the other tree's `Z`
  happens once instead of per node. One bit per cell, ~3 MB at 5000x5000.
* `_prob_other_tree_at_leaf(clique, leaf)` returns `P_other(parent state -> s)` for `s = 0, 1`.
  The roles swap between the trees: our clique index is a *leaf* over there, and our leaf index
  selects *which of its cliques* (hence which alpha and nu) applies to the cell.
* `_add_pseudo_likelihood_term` now handles leaves: their conditional is the product of the two
  trees' parent terms, normalised. The other tree's factor is constant in this tree's
  parameters but sits inside the normaliser, so it does not cancel.
* Note the data likelihood does **not** appear. The pseudo-likelihood is that of the field
  `p(Y, Z | theta)`; `P(data | Y)` carries no theta and drops out of the Hastings ratio. (I
  briefly talked myself into thinking it was needed; it is not.)
* `--leaf_pseudo_likelihood` turns it on.

Builds clean, 123 tests pass.

### 10.2 Measured: better on the control, worse on the real thing

`species_mean_log_nu` with Y and Z pinned at the simulated truth:

| | coupled scenario (true −0.446) | control, molecules tree neutral (true −0.500) |
|---|---|---|
| Phase 1, internal nodes only | −0.402 (bias **+0.044**) | −0.355 (bias **+0.146**) |
| Phase 2, with leaf terms | −0.753 (bias **−0.307**) | −0.529 (bias **−0.029**) |

On the control — where the molecules tree is neutral so the true other-tree factor is
analytically `{0.5, 0.5}` — Phase 2 is clearly the better estimator, on `alpha` too
(0.585 vs 0.523 against a true 0.600). That is the expected result.

On the genuinely coupled field it is clearly the worse one.

### 10.2b …but much better on the full pipeline, which is what matters

The table above pins Y and Z at the truth. That is a diagnostic regime, not how the model is
used. On the **full pipeline**, with Y and Z inferred:

| | Phase 1 (leaves off) | Phase 2 (leaves on) | true |
|---|---|---|---|
| species `mean_log_nu` | −1.002 (bias −0.556) | **−0.803** (bias −0.357) | −0.446 |
| molecules `mean_log_nu` | −1.054 (bias −0.588) | **−0.932** (bias −0.466) | −0.466 |
| `epsilon` | 0.0582 | **0.0082** | 0.001 |
| `epsilon_simple_model` | 0.3201 | **0.2061** | 0.200 |
| Y accuracy | 0.681 | **0.920** | — |
| Y predicted frac(1) | 0.185 | **0.388** | 0.408 |

Large and across the board. Y recovery goes from 68% to 92%, the predicted fraction of ones
lands within 0.02 of the truth instead of being less than half of it, and both error rates stop
absorbing the mismatch. This is precisely the mechanism Phase 2 was built for: the leaf terms
anchor theta to Y, which the data constrains, breaking the Z feedback loop of section 6.

The two regimes therefore disagree. With a *given, perfect* field the internal-node terms are
already an excellent estimator and the leaf terms only add bias; with an *inferred* field the
leaf terms are what stop the whole thing sliding. Both can be true at once.

Note also that this does not validate the lookup: almost any leaf term, even one reading the
wrong cell, would anchor theta to Y and break the feedback loop. The improvement is real but it
is not evidence of correctness.

### 10.3 Narrowing it down

A diagnostic build that forces the other-tree factor to a constant `{0.5, 0.5}` on the coupled
scenario gives bias **−0.046** (species), versus **−0.307** with the real factor and **+0.044**
with no leaf terms at all. So:

* the factor is genuinely live — neutralising it moves the answer a long way;
* the leaf-term mathematics is right, or the control would not have improved;
* but feeding in the *correct* factor is worse than feeding in a deliberately wrong constant,
  which is backwards.

The one thing the control does **not** exercise is the cross-tree lookup itself: when the other
tree is neutral its factor is `{0.5, 0.5}` whatever parent state is read, so a wrong index is
invisible there. And that lookup is precisely the part that matters once the trees are coupled.
The mapping was re-derived on paper for both trees and checks out (clique index ↔ the other
tree's leaf; our leaf index ↔ the other tree's clique; `Z` dimensions ordered by each tree's
own `_dimension`), but paper is not a test.

**At the time of writing this section `--leaf_pseudo_likelihood` was off by default; see
§10.5 for why that was reversed.** The next step is a direct test of the lookup — the natural one is to check `_prob_other_tree_at_leaf` against the per-dimension
terms that `TMarkovField::_calculate_log_prob_field` computes for the same cell, since those
two must agree by construction.

### 10.3b The lookup test — verdict: the lookup is correct

`TTree::check_leaf_lookup_consistency()`, run with `--check_leaf_lookup` at the first
iteration (once Z is populated), aborting on any mismatch.

It could not be a gtest: `_prob_other_tree_at_leaf` needs `get_binned_branch_length`, which
needs initialised stattools parameter storage, which needs the whole DAG. (That is also why
`tests/TUpdate_Z_Tests.cpp` was an empty stub for so long.) So it is a runtime check on the
real trees instead.

The point is that it recomputes the quantity by a route that **does not assume the role swap**.
It builds the cell's index in Y space from `(clique, leaf)`, and from that derives the other
tree's clique via `get_clique(IndexArray)` — the same call the Y sweep makes in
`_calculate_log_prob_field` — its leaf node, and its parent's state, then compares the
resulting probability pair cell by cell.

Result: **all 16384 cells agree, for both trees.**

Mutation-checked, because both scenarios are 128x128 and a clique/leaf swap would not even go
out of range there:

| mutation | caught |
|---|---|
| swap the clique and leaf roles in `_prob_other_tree_at_leaf` | 16256 of 16384 disagree — the 128 that survive are exactly the diagonal where clique == leaf |
| swap the two dimensions of the Z index in the cache | 5774 of 16384 disagree |

(Worth remembering: square scenarios hide index bugs. A non-square test scenario would be a
cheap safeguard.)

### 10.3c So what explains the pinned-field anomaly?

Not a bug in the lookup. The most likely remaining explanation is that the leaf terms
**couple the two trees' parameter estimates**: the species leaf term contains the molecules
tree's alpha and nu, which are themselves only weakly identified and wandering. With
internal-node terms only, the two trees' parameters are estimated independently and each is
near-unbiased given a perfect field. Adding leaf terms makes each tree's estimate depend on the
other's current error. That is a hypothesis, not a demonstrated fact.

It also does not matter much in practice: with an inferred field the leaf terms are worth far
more than they cost (§10.2b).

### 10.4 The joint (Y, Z) wiring was not done

Deliberately not attempted this round. Beyond the four items in section 9.2 it also needs the
OpenMP strategy inverted — today the parallel loop runs over cells *within* a clique, whereas a
block draw is sequential inside a clique and must be parallelised *across* cliques, which means
per-thread `_clique_last_dim`/sheet state instead of one shared set — and it needs `K >= M` so
a clique is never split across sheets. Landing that together with an unvalidated Phase 2 would
have made it impossible to attribute whatever came out.

### 10.5 Confirmed, and turned on by default

The full-pipeline measurement of §10.2b was repeated with the current binary and the flag set
explicitly, because the original run had been launched while the flag defaulted on and the
binary was rebuilt underneath it. It reproduces:

| | Phase 1 (off) | Phase 2 (implicit) | Phase 2 (explicit flag) | true |
|---|---|---|---|---|
| species `mean_log_nu` | −1.002 | −0.803 | −0.825 | −0.446 |
| molecules `mean_log_nu` | −1.054 | −0.932 | −0.939 | −0.466 |
| `epsilon` | 0.0582 | 0.0082 | 0.0080 | 0.001 |
| `epsilon_simple_model` | 0.3201 | 0.2061 | 0.2061 | 0.200 |
| Y accuracy | 0.6806 | 0.9195 | 0.9196 | — |
| Y predicted frac(1) | 0.1850 | 0.3884 | 0.3889 | 0.408 |

The two Phase 2 runs agree to within MCMC noise, so the earlier numbers were valid.

With the lookup verified cell by cell (§10.3b) and this reproduced, the leaf terms are now
**on by default**; `--no_leaf_pseudo_likelihood` restores the internal-nodes-only behaviour.

The pinned-field oddity of §10.2 is unchanged and unexplained beyond the hypothesis in §10.3c.
It is recorded in `cli.h` next to the flag so it is not forgotten.
