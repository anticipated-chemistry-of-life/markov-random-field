# Outside structures enter the unknown weight by likelihood ratio

_Supersedes the section "Rows outside the molecule tree" of [ADR-0018](0018-mass-spec-candidates-enter-by-likelihood-ratio.md). The rest of ADR-0018 stays valid._

An outside structure is a proposed structure that is not a molecule leaf. ADR-0018 counted it toward $n_f$ and discarded its likelihood ratio. Now the unknown weight of feature $f$ is

$$W_{U,f}=\rho\,n_f+\tau\sum_{j\in\text{outside}(f)}\frac{g_1(x_j)}{g_0(x_j)},$$

where $\tau$ is the outside rate. The input is `--ms_outside_rate`. It is required with `--ms_data`, has no default, and must be $\ge 0$. With $\tau=0$ the model is ADR-0018 exactly.

## Why

On the real export, 82% of proposed structures are outside structures, and in 73% of features the best-scored one is an outside structure. With ADR-0018, the model did not see these matches. The unknown weight was therefore too low, and weak candidates took features that a strong outside match explains. That puts false presence into the field.

ADR-0018's product argument extends with no change. If the true structure is outside structure $j$, its weight is its prior times $g_1(x_j)/g_0(x_j)$, as for a candidate. A candidate's prior in `acol` is its detection odds: $f/(1-f)$ if $Y=1$ and $\kappa/(1-\kappa)$ if $Y=0$. An outside structure has no field cell, so its detection odds cannot come from $Y$. The outside rate $\tau$ is one fixed value for all outside structures. It is the odds that the molecule of an outside structure is present in the species and passes the filter.

$\rho\,n_f$ now covers only the explanations that are not in the list: an unlisted structure, or no molecule at all.

## Consequences

- $W_{U,f}$ is still one number per feature, and it is known at load time. The bins, the moves and the field draws do not change.
- An outside structure stays inside the unknown molecule, so it is exempt from the rule "one molecule explains at most one feature per run". Two features of one run can share an outside structure. This is a small misspecification, and we accept it.
- The assignment posterior divides the unknown molecule's fraction. Each outside structure gets a row with $P(U)\,\tau\,\mathrm{LR}_j/W_{U,f}$. The remainder, $P(U)\,\rho\,n_f/W_{U,f}$, goes in a row under the reserved name `<unlisted>`. The rows of a feature still add up to 1.
- Source of $\rho$: a fit of $c=\rho\,n_f$ on held-out ground-truth features in feature-probabilities (`evaluate-calibration`), with outside structures counted as listed structures. All ground-truth features have a molecule, so the fit is a lower bound on $\rho$. The fit also gives every listed structure a prior of 1, so it agrees with `acol` only where the detection odds are about 1.
- Source of $\tau$: no data set measures it. It is a plug-in value: the mean detection odds of the in-tree candidates in a pilot run, checked against the values $\{0.01, 0.1, 1\}$.

## Considered options

- **$\tau = 1$, no input.** Rejected. It treats every outside molecule as present and detected, so outside structures would get more weight than in-tree candidates with a low $P(Y=1)$.
- **Estimate $\tau$.** Deferred, for the same reason as $\rho$ in ADR-0018.
- **Attach outside structures to the molecule tree.** Rejected. About 80% of rows would need a placement, and the tree would change.
