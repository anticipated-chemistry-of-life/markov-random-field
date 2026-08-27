# The two-tree product is unnormalised, which biases every parameter toward small nu

The field's joint density is the product of the two trees' likelihoods — `TMarkovField::_calculate_complete_joint_density` sums the two trees' log-densities and nothing else, and `TTree::_update_nu_or_alpha` forms its acceptance ratio from one tree's likelihood alone. Each factor is a properly normalised tree likelihood over its own variables. Their product is not normalised, because both are functions of the same field:

```
sum over Y, Z_s, Z_m of  p_s(Y, Z_s | theta_s) p_m(Y, Z_m | theta_m)
    = sum over Y of  p_s(Y | theta_s) p_m(Y | theta_m)
    = C(theta_s, theta_m)
```

`C` is not 1 and, more importantly, it moves with the parameters. A sampler that treats the product as a density is therefore targeting `p_correct(theta) * C(theta)`, and is pulled toward whatever maximises `C`. Since `C` is an overlap between the two trees' field distributions, it grows as either tree becomes more strongly correlated down its branches — that is, as `nu` falls. Nothing in the acceptance ratio pushes back, so `log_nu` and with it `mean_log_nu` drift downward without bound.

This is the doubly intractable problem: `C` is a sum over `2^(n_species_leaves * n_molecule_leaves)` fields.

## Why validation did not catch it

Under a neutral molecules dimension (ADR-0001) `p_m(Y)` is the same for every field, so `C` factors into a constant and the bias is exactly zero. The independent-field harness pins molecules neutral by construction, so it is blind to this by design, not by oversight. This is the concrete instance of the blind spot ADR-0001 records.

## Evidence

`model_validation/diagnose_normaliser.py` enumerates every field for two balanced trees, so `C` is computed exactly rather than estimated. With the molecules dimension neutral, the C++'s objective and the correctly normalised one peak at the same place. As molecules leaves neutrality the two separate, and the gap grows with both the strength of the molecules dimension and the size of the field: at `log_nu_molecules = -3` the peak sits 0.58 below the truth for a 2x2 field and 1.38 below for a 4x4 one. A Metropolis chain started at the truth reproduces the drift and stops drifting when `-log C` is added to the ratio, with nothing else changed.

## Consequences

The drift is a property of the target, not of the sampler, so it will not yield to more iterations, better proposals, or a tighter prior on `mean_log_nu`. Any run with two active trees is affected; the more phylogenetic signal the molecules tree carries, the worse it is.

Options, none yet chosen: estimate the gradient of `log C` by sampling from the model (a stochastic-approximation / persistent-contrastive-divergence scheme); replace the product with a genuinely normalised model, for instance a single process on the product tree; or accept the product as a composite likelihood and stop interpreting `nu` as a rate, reporting only quantities that are invariant to the bias. The first keeps the current model, the second changes it, the third changes what is claimed about the output.
