# Neutralise the molecules dimension to validate against a single-tree reference

_Superseded by [ADR-0005](0005-each-tree-owns-its-leaf-level-field.md), which gives each tree its own leaf-level field. Under that model the independent reference is exact with **both** trees active, so neutralisation is no longer the price of having a reference at all. It survives as a cheap regression rung, and everything below still describes what that rung does and what it must pin. What it can no longer be is the rung that exercises the whole model: by derivation 3 of ADR-0005 a neutralised tree contributes nothing to the new error probability, so such a rung reaches that parameter through one tree where the model has two._

To check the inference code against an independent implementation, we need a reference simulator that is simple enough to be obviously correct — a one-pass-down sample of a single tree — but the model requires exactly two trees, and the field conditional is a product over both. Rather than adding a single-tree mode to the C++, we make the molecules dimension **neutral**: `alpha = 0.5` and `log_nu = 5.0`, which puts nu past the threshold above which `TMatrices::set_lambda` switches to the stationary distribution, giving transition rows of exactly (0.5, 0.5). The molecules factor is then a constant that cancels out of every field conditional, and the model reduces exactly — not approximately — to the species tree process the reference implements.

## Considered options

Adding a genuine single-tree code path to the C++ would test a code path that production never runs, so a passing test would prove little about the two-tree model that actually ships. Neutralisation exercises the real code path with real two-tree machinery.

## Consequences

The reduction is exact only while the pinning holds, so validation runs must fix `molecules_alpha`, `molecules_log_nu`, `molecules_mean_log_nu`, `molecules_var_log_nu` *and* `molecules_branch_lengths` via `--<name>.update false`. Leaving branch lengths free is harmless to correctness but produces a meaningless random walk in the traces, since under neutrality every bin yields the same matrix.

The neutrality assumption is load-bearing for the entire validation result, so it is asserted rather than trusted: the generator checks `nu_molecules > 25` explicitly, and the harness runs one rung twice with different `Z_molecules` draws and requires the species posteriors to come out identical.
