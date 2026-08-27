# The transition grid is built by recursion, not evaluated bin by bin

`TTransitionGrid` evaluates the closed form of `expm(Lambda * t)` exactly twice — at `Delta / 2` and at `Delta` — and then reaches every other bin by repeated multiplication, `P_k = P_{k-1} * P_Delta`. Evaluating the same closed form directly at each bin's centre `Delta * (k + 0.5)` would be shorter, faster, and produce the same numbers. It must not be done, and the reason is not visible from inside `src/`.

`model_validation/src/independent/field.py` evaluates the closed form at bin centres **because** the C++ does not. That is the entire content of the independence claim the validation harness rests on: two implementations that reach the same quantity by different routes, so that an error in one is exposed rather than mirrored. Collapse the C++ onto the closed form and the Python stops being a check and becomes a copy — a mistake in the discretisation would then be reproduced identically on both sides and pass every rung. ADR-0001 argues that the reference must be independent; this is where that independence actually lives, now that both sides share the closed-form expression itself.

## Considered options

Keeping `arma::expmat` in the implementation would preserve the divergence more obviously, since the two sides would then share no formula at all. It was dropped because armadillo's only use in `src/` was these two calls, and a 2x2 generator with eigenvalues `0` and `-nu` has a closed-form exponential that needs no numerical routine. The cross-check is preserved instead by keeping `arma::expmat` in `TTransitionGrid_Tests` as the numerical oracle: the implementation walks the grid by multiplication, the test evaluates a numerical matrix exponential at each bin's own branch length, and the two are required to agree.

Evaluating the closed form per bin and deleting the recursion was rejected for the reason above.

## Evidence

The recursion is not, at present, accumulating error worth guarding against: over `alpha` in `[0.01, 0.99]` and `log_nu` in `[-8, 3.2]`, the grid built by recursion agrees with direct evaluation at the bin centres to `1.6e-15` at `n_bins = 10`. The guard is insurance against a larger grid — `--n_bins` is a command-line option and the accumulated error grows with the number of multiplications — rather than protection against a live problem.

Replacing `arma::expmat` with the closed form for the two base matrices moved the grid by at most `2.1e-15` across that same sweep. At the parameters of the standard validation scenario the base matrices differ bit-for-bit in 250 of 256 cases, by at most `2.2e-16`, and a 200-iteration rung-1 chain at `--numThreads 1 --fixedSeed 42` nonetheless produced byte-identical traces, posteriors and `Y`: the numbers changed, no accept/reject decision did.

## Consequences

`arma::expmat` survives in `TTransitionGrid_Tests` although nothing in `src/` uses it. That is deliberate and should not be tidied away.

A future reader who wants to simplify `TTransitionGrid` by looping `transition(alpha, nu, bins.grid_branch_length(k))` over the bins is making a change to the *validation architecture*, not to an implementation detail, and it cannot be made without deciding what `field.py` should do instead.
