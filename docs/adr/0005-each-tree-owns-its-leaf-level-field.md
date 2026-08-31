# Each tree owns its leaf-level field, and the field is a noisy AND of the two

ADR-0002 records that the field's joint density is formed as the product of the two trees' likelihoods, that the product is not normalised because both factors are functions of the same field, and that the missing constant `C(theta_s, theta_m)` drags `log_nu` downward without bound. The cause is one variable serving two processes: the field is simultaneously the leaf state of the species tree process and the leaf state of the molecule tree process.

So stop sharing it. Each tree's node state extends down to its own leaves, giving that tree a complete leaf-level view of the field — its **tree field**. The field becomes a noisy reconciliation of the two, with a single **error probability** `omega`. The joint then factors as

```
p(Z_s | theta_s) * p(Z_m | theta_m) * p(Y | Z_s, Z_m, omega) * p(L, D | Y)
```

Every factor is a proper conditional density and nothing is conditioned on its own descendants, so this is a directed factorisation and it sums to 1 by construction. `C` is identically 1 — for every parameter value, not merely at a neutral point. **ADR-0002's bias is structurally gone, not mitigated.** And because `theta_s`'s Markov blanket is now `Z_s` alone, the alpha and nu acceptance ratios read one tree, which is the independence ADR-0002 says the current model lacks.

ADR-0002 lists three remedies: estimate the gradient of `log C` by sampling, replace the product with a genuinely normalised model, or accept the product as a composite likelihood and stop interpreting nu as a rate. **This is the second.** ADR-0002 sketches that option as "a single process on the product tree"; the construction here reaches a normalised model by a different route, keeping two independent tree processes and paying for normalisation with one extra latent layer and one extra parameter instead.

Three derivations follow. Each is load-bearing for a decision made in the spec, and each is stated here rather than left to be rediscovered from the code.

## Derivation 1 — the link is an AND over independently corrupted tree fields

The link in the factorisation above is defined cell by cell: `p(Y | Z_s, Z_m, omega)` is a product, over leaf pairs, of the table below. Throughout these derivations `Z_s` and `Z_m` name the *cell* of each tree field at one leaf pair — one entry of the leaf block of the node state written with the same letter — and `Y` names the field's cell there. Corrupt each tree field independently with probability `omega`, then take the AND of the results. Since a corrupted cell reads 1 with probability `1 - omega` when it is truly 1 and `omega` when it is truly 0, and the two corruptions are independent, the AND is a product of two independent factors:

```
P(Y = 1 | Z_s = 1, Z_m = 1) = (1 - omega)^2
P(Y = 1 | Z_s = 1, Z_m = 0) = (1 - omega) * omega
P(Y = 1 | Z_s = 0, Z_m = 1) = omega * (1 - omega)
P(Y = 1 | Z_s = 0, Z_m = 0) = omega^2
```

The AND is a real modelling commitment and not an algebraic consequence of anything above it. It says a molecule occurs in a species when *both* phylogenies allow it — the lineage can make the compound and the compound is producible — which makes the field **sparser than either tree field**: a one needs both trees to agree, a zero is cheap. It is stated rather than inherited, and it sits behind a link-policy seam, because the normalisation win comes from the *structure* — each tree owning its own leaf-level field — and not from this particular table. The link can be changed later without revisiting anything argued above.

Two boundaries follow from the table and both are why `omega` is constrained to the open interval `(0, 0.5)`. At `omega = 0` the link is the deterministic AND, the block update hits `log 0`, and a field cell at one pins both tree fields. Above `0.5` the tree fields are anti-correlated with the field, which is meaningless as an error model and is a genuine second mode for a sampler to find.

## Derivation 2 — the link likelihood collapses to six counters

The table depends on `(Z_s, Z_m)` only through the sum `k = Z_s + Z_m`, so there are three distinct probabilities rather than four:

```
P_k = (1 - omega)^k * omega^(2 - k)        P_2 = (1-omega)^2, P_1 = omega(1-omega), P_0 = omega^2
```

The link's contribution to the log-likelihood is therefore a function of `omega` and six integers, `n(k, y)` counting cells with `Z_s + Z_m = k` and `Y = y`:

```
log p(Y | Z_s, Z_m, omega) = sum over k in {0,1,2} of  [ n(k,1) * log P_k  +  n(k,0) * log(1 - P_k) ]
```

This is what makes the `omega` Metropolis move O(1) in the number of cells rather than O(cells) — the same trick the simple error model's disagreement count already plays for its error rate. The block update knows its own old and new `(k, y)`, so the counters are maintained incrementally.

The collapse also leaves something testable behind. Three Bernoulli probabilities are pinned by one parameter, so two parameter-free constraints hold for **every** `omega`:

```
P_1^2 = P_0 * P_2                 (1)
sqrt(P_0) + sqrt(P_2) = 1         (2)
```

Both are checkable directly against the counters, using `p_k = n(k,1) / (n(k,0) + n(k,1))`, with no parameter estimated first. Equivalently, (1) says `log P_k` is affine in `k` and its second difference vanishes.

The two test different assumptions, and knowing which is which is the point of writing them down. Allow the two trees distinct error probabilities `omega_s` and `omega_m` and work with the four-cell table instead, writing `P_{z_s z_m}` for `P(Y = 1 | Z_s = z_s, Z_m = z_m)`. Then

```
P_11 * P_00 = P_10 * P_01
```

still holds identically — it is the signature of *independent corruption followed by an AND*, and it survives any pair of error probabilities. What fails is (2): by Cauchy–Schwarz, `sqrt(P_00) + sqrt(P_11) <= 1` with equality exactly when `omega_s = omega_m`. So (2) is the test of the **one shared `omega`** assumption, and the four-cell product identity is the test of the **link structure**.

Constraint (1) is the six-counter shadow of the four-cell identity, and bucketing by `k` pools the two middle cells, which costs more than it first looks. A violation of (1) means *either* that the link is not an independent-corruption AND *or* that the two trees do not share one `omega`, and six counters cannot say which; separating those needs the middle cells counted apart, which is eight counters rather than six. Worse, (1) is necessary but not sufficient. The pooled middle rate is a weighted average of `P_10` and `P_01`, and their geometric mean always lies between them, so for **any** unequal pair of error probabilities there is a split of the cells across the two middle buckets at which (1) holds exactly — at `omega_s = 0.05` against `omega_m = 0.2`, a `0.31 / 0.69` split does it. Constraint (2) has no such blind spot and is the one to watch: on that same pair its left side misses 1 by `2.8e-2`, and at `0.02` against `0.45` by `1.7e-1`.

Both ship as **traced diagnostics, not as failing tests**. A violation means the link is wrong, which is a finding rather than a regression.

## Derivation 3 — the marginal field rate is a product, and predicts a ridge

Take one leaf pair with the tree fields at their stationary rates, `P(Z_s = 1) = alpha_s` and `P(Z_m = 1) = alpha_m`. A corrupted tree field cell reads 1 with probability — call this the tree's **corrupted rate** —

```
a~ = omega + (1 - 2*omega) * alpha
```

and the two corrupted cells are independent, so the expected density of the field is a **product**:

```
P(Y = 1) = a~_s * a~_m
```

The overall density of the field therefore says nothing about how the rate splits between the two trees. Only the phylogenetic correlation structure — how cells covary along each tree — separates `alpha_s` from `alpha_m`. When that signal is weak the posterior develops a ridge along constant `a~_s * a~_m`, a hyperbola in `(alpha_s, alpha_m)` at fixed `omega`.

This is the same shape as the ADR-0002 defect in one respect and unlike it in another, and the difference matters. Both are the two trees trading against each other through a single quantity the field pins down. But ADR-0002's `C` biases the target — the posterior peaks in the wrong place, and more iterations make the answer more confidently wrong. This one does not: the posterior is correct and merely flat in one direction. It is an identifiability limit to be reported, not a defect to be fixed, and the remedy is a rung that deliberately weakens both trees and asserts that the *product* is recovered while the individual alphas are not.

`omega` sits inside both factors, which extends the flat direction to three parameters rather than two. The degenerate case is worth naming because it is the old validation setup: at `alpha = 0.5` the expression above gives `a~ = 0.5` for **every** `omega`. A neutralised tree's own factor therefore carries no information about the error probability. `omega` is not thereby lost — it stays identified through the *other* tree, since with molecules neutral `P(Y = 1 | Z_s = 1) = 0.5 * (1 - omega)` against `P(Y = 1 | Z_s = 0) = 0.5 * omega`, a contrast of `0.5 * (1 - 2*omega)` — but one of the two trees stops speaking to it. That is one more reason neutralisation (ADR-0001) is retired as the primary rung.

## Evidence

Each derivation was checked against a brute-force enumeration sharing no algebra with the closed forms above: the link table rebuilt by summing over both corruption events, the six-counter log-likelihood against a cell-by-cell sum over a random configuration, and the marginal rate against an explicit sum over the four `(Z_s, Z_m)` states. The three agree to floating point — worst absolute differences `0`, `0` and `1.1e-16` over `omega` in `[1e-6, 0.499]` and alphas in `[0.05, 0.9]` — and the joint sums to `1.000000000000000` on enumerable `(Z_s, Z_m, Y)` cases at several parameter settings.

The constraints behave as claimed. Under one shared `omega`, `P_1^2 - P_0 * P_2` is at worst `1.4e-17` and `sqrt(P_0) + sqrt(P_2) - 1` is exactly `0`. Under unequal error probabilities the four-cell identity survives to `8.7e-19` while (2) misses 1 by `2.8e-2` at `(0.05, 0.2)` and `1.7e-1` at `(0.02, 0.45)`; the middle-bucket split at which (1) holds regardless is `0.3145 / 0.6855` for the first of those pairs.

These check the algebra in this record, not an implementation. The implementation-side versions — the block update against brute force, the six counters against naive recomputation, and the identity as a traced diagnostic — belong to the field-math header and are specified there.

## Considered options

**Estimating the gradient of `log C` by sampling** — ADR-0002's first option — keeps the current model and the current interpretation of every parameter, which is its attraction. It was rejected because it buys a biased-but-familiar model at the cost of a permanently harder sampler: an inner sampling loop inside every parameter move, a stochastic approximation to tune, and a correctness argument that is asymptotic rather than structural. The defect would become a controlled error rather than an absent one, and every future run would carry that qualification.

**Accepting the product as a composite likelihood** — ADR-0002's third option — is the cheapest, and it is honest as far as it goes. It was rejected because it forfeits the quantity the project exists to report: `mean_log_nu` as a switching rate. A model that cannot be asked its own headline question is not the model to keep.

**A single process on the product tree**, the construction ADR-0002 names for its second option, normalises correctly but merges the two phylogenies into one object. It was rejected because it destroys the per-tree parameters — there is no longer a `nu` belonging to the species tree — and because the product tree is quadratically larger than either factor. Tree fields keep both trees intact and pay instead with one latent layer and one scalar.

**A shared field with an explicit normalising term** was not seriously available: `C` is a sum over `2^(n_species_leaves * n_molecule_leaves)` fields, which is the doubly intractable problem ADR-0002 names.

## Consequences

**The field is no longer where the two trees meet; the link is.** Anything that read the field as "the species tree's leaf state" is now reading the wrong variable, and there are three variables at each leaf pair where there was one. This is the change that makes the sweep a single 8-state block rather than three single-site updates: with a small `omega`, a field cell at one pins both tree fields, and given both at one the field stays at one with probability near one. That triple is metastable under single-site Gibbs and would present as slow mixing rather than as a bug.

**One new parameter, and it is confounded with two things already in the model.** `omega` trades against the alphas through derivation 3, and against the simple error model's own misreport probability, since both describe noise standing between a latent one and an observed one. Both limits are known in advance and are watched rather than discovered — which is the whole reason the six counters, the posterior of each tree field and the full joint density are traced.

**ADR-0001 and ADR-0002 are superseded, not deleted.** ADR-0002 is the reason this work happened and its evidence is what justifies the cost; ADR-0001 records a validation architecture that was correct for the model it was built against. Neutralisation is retired as the primary rung because this model admits an exact independent reference with **both** trees active — two independent top-down draws and a four-cell table — but a neutralised rung is kept as a cheap regression check.

**ADR-0003 is unaffected.** The recursion in `TTransitionGrid` and the independence it protects are properties of one tree's process, and each tree still runs exactly that process.

**The old model is deleted when this branch merges, and the merge-base is tagged.** Git is the flag. Keeping a known-biased code path compiled and tested to preserve the option of comparing against it pays a permanent cost for what a checkout gives free.

**Whether the new model predicts better is not settled here.** The merge gate is correctness: the joint sums to 1 on enumerable cases, a chain started at the truth shows no drift, and the rung ladder recovers what it should. Held-out predictive performance on real LOTUS records is a later question, and it decides whether the old model stays deleted — not whether this branch merges.
