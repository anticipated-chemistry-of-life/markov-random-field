# Mass-spec candidates enter by likelihood ratio

_Supersedes [ADR-0012](0012-mass-spec-scores-enter-only-relative-to-their-features-best-row.md). The 1/16-nat grid and the argument that every assignment move reads only ratios within one feature carry over from it; the reference and the unknown molecule's value do not._

A candidate's weight is its **likelihood ratio**, its correct density divided by its incorrect density. The **unknown weight** of a feature is $\rho\,n_f$, where $\rho$ is the unknown rate and $n_f$ the number of structures proposed for the feature, candidates or not. ADR-0012 divided the correct density by the feature's best row and gave the unknown molecule one fixed value. That is exact only if the incorrect density is flat in CSI:FingerID score, and it is not. feature-probabilities now fits and exports both densities (its ADR 0004).

## Why the ratio, and why $n_f$

Take a feature with proposed structures $j$ and scores $x_j$. Suppose the true structure sits among them like any other row and rows are independent given the truth. Then

$$P(A=m\mid x)\propto \pi_m\,\frac{g_1(x_m)}{g_0(x_m)}\prod_j g_0(x_j),\qquad P(A=U\mid x)\propto \pi_U\prod_j g_0(x_j).$$

The product cancels. Each candidate's weight is its prior times its likelihood ratio, and the unknown molecule's weight is a prior odds. Dividing by the best row instead makes the unknown molecule's effective weight $c$ times that row, so it changes from feature to feature for no reason in the data.

The odds scale with $n_f$. Proposed structures drawn uniformly from $M$ give $P(\text{list}\mid U)/P(\text{list}\mid m)=n_f/M$. A uniform prior over the listed structures gives $n_f\,\pi_U/(1-\pi_U)$. So the unknown weight is $\rho\,n_f$ with one $\rho$ for every feature, not one constant.

The held-out comparison in feature-probabilities (`evaluate-calibration`, issue #38) measured this. Over 73,768 ground truth features, the likelihood ratio raised top-1 from 0.358 to 0.496 and lowered log-loss from 2.04 to 1.55. The best-row weighting put every one of 6.45 M candidate rows below a posterior of 0.1.

## Bins stay per feature

Each weight is binned relative to the largest weight in its feature, the unknown weight included: $b=\mathrm{round}(16(\ln w_{\max,f}-\ln w))$. The grid, its 8 bits and its floor are ADR-0012's. The unknown molecule's bin is now per feature instead of global. A likelihood ratio above 1 needs no special case, and a weight more than $255/16$ nats below its feature's largest is clamped to the floor.

## Rows outside the molecule tree

A proposed structure that is not a molecule leaf counts toward $n_f$ and nothing else. Exactly, its likelihood ratio belongs to the unknown weight, since the unknown molecule stands for every explanation outside the tree. Adding it needs a prior for such rows, which is not settled. Simulated data has no such rows, so the limitation costs nothing there. On real data it understates the unknown weight whenever a strong match lies outside the tree. It must be resolved before the model is run on real data.

## Considered options

**Keep the best-row division.** Rejected on the evidence above.

**One offset for the whole file.** Bins would sit below the largest weight in the file. Rejected: the moves never compare weights across features, so a global scale only wastes the grid's range.

**One constant unknown weight.** Rejected: it ignores that a feature with more proposed structures gives each a smaller share of the prior. In the held-out comparison, one fitted constant came out at 32, while $\rho\,n_f$ at the observed odds and mean $n_f$ would be about 122.

**Estimate $\rho$.** Deferred. $\rho$ is an input (`--ms_unknown_rate`, required with `--ms_data`, no default). Estimating it would add a parameter and a move.

## Consequences

- The loader reads `log_correct_density` and `log_incorrect_density` instead of `calibration_score`, and `inference_lotus/ms_prep.py` passes them through.
- `--ms_unknown_fraction` is removed, not aliased.
- The feature term is still a weight, not a probability, and it still cancels from every field draw and every filter-probability move, as in ADR-0012.
