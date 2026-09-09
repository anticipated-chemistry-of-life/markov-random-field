# The mixing cost of dropping the block update, measured

ADR-0005 built the leaf layer's eight-state block update on an argument. At a small error
probability a field cell at one pins both tree fields; given both tree fields at one the field
stays at one with probability near one; that triple is metastable under single-variable draws, and
the failure "would present as slow mixing rather than as a bug". Issue #68 deleted the block update
and draws a leaf with its tree, then the field on its own. This is the measurement that replaces
the argument.

It is written to be carried: ADR-0008 is issue #70's, and this is the measurement it reports as its
consequences.

**The cost is real, it is confined to the field, and it is a small-`omega` effect.** At
`omega = 0.005` the field's own factor of the joint density takes **4.7 times longer** to
decorrelate. At `omega = 0.1` the same number is **1.2**. Nothing else in the model mixes
measurably differently, and both tree field posteriors land in the same place under both binaries.

## What was run

`run.sh` in this directory, at two error probabilities. Two binaries: `reference` is `a3d1ae6`, the
commit before the block update was deleted, and `current` is `f13fef9`, the commit that deleted it.
Both link the same coretools and stattools checkout.

| | |
| --- | --- |
| fixture | `../s_balanced_255_m_balanced_511`, 128 x 256 = 32 768 leaf pairs |
| data | one set, simulated by the reference binary at seed 20260909 |
| chains | 20 000 iterations after 2 000 burn-in, `--numThreads 1` |
| replicates | 4 per binary at `omega = 0.005`, 2 per binary at `omega = 0.1` |

**The data is fixed and the seed is not.** The two binaries cannot be compared seed for seed: a
leaf's cell uniform used to come from the field stream, drawn by the block update, and now comes
from the node-state stream at the leaf's own linear index, with the field drawing a pass of its
own. One seed therefore drives two different sequences of draws, and a seed-matched pair would
differ for a reason that says nothing about mixing.

Be precise about which of issue #69's two reasons applies here. It names both "the storage layer
and the pattern of random-number consumption", and against the *branch point* both did change. But
the reference chosen here is `a3d1ae6`, which already carries the hash-map storages — those landed
at `6d42857`, three commits earlier, and `git diff a3d1ae6 f13fef9 -- src/storages/` is empty. So
against **this** reference it is the draw pattern alone, and that on its own is enough: it is why
the seed is deliberately let go, with every number below a median or a spread over replicates
rather than a chain-against-chain difference.

## The cost, at `omega = 0.005`

Integrated autocorrelation time, by factor of the joint density trace. `ratio` is current against
reference; `drift` is how far the second half of the trace sits from the first, in trace standard
deviations.

| factor | reference `tau` | current `tau` | ratio | reference drift | current drift |
| --- | ---: | ---: | ---: | ---: | ---: |
| `data` — the field's own | 63.2 | 295.1 | **4.67** | −0.03 | 0.05 |
| `link` | 1371 | 1604 | 1.17 | −0.25 | 0.49 |
| `species_node_state` | 3630 | 3718 | 1.02 | −1.12 | −1.00 |
| `molecules_node_state` | 1182 | 1377 | 1.17 | 0.08 | −0.42 |
| `joint_density` — the total | 2295 | 2829 | 1.23 | −0.95 | −0.85 |

**Read the `data` row and not the total.** `data` is `log p(L, D | Y)`: the LOTUS and simple error
model likelihoods, which are functions of the field and of their own two parameters and of nothing
either phylogeny carries. It is the field's own instrument, and it is the only factor whose drift
says the chain had reached the distribution its autocorrelation time describes. The two node-state
factors are dominated by the phylogenetic parameters, which had not converged in either binary at
20 000 iterations — their drift of about one standard deviation says so — and they drag the total
along with them. The total's 1.23 is that dilution, not a finding.

The autocorrelation function of the `data` factor, averaged over the replicates:

| lag | 1 | 2 | 5 | 10 | 25 | 50 | 100 | 200 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| reference | 0.567 | 0.527 | 0.451 | 0.369 | 0.251 | 0.170 | 0.108 | 0.052 |
| current | 0.934 | 0.888 | 0.799 | 0.714 | 0.586 | 0.476 | 0.362 | 0.233 |

The shape is the metastability ADR-0005 described. One pass of single-variable draws moves the
field barely at all — `acf(1)` is 0.93 where the block update left it at 0.57 — and the trace is
still 0.23 correlated 200 iterations later, where the reference has reached 0.05.

**And it does not buy iteration time.** The current binary is **1.25 times slower per iteration**
(11.20 ms against 8.95 ms), because three passes over the leaf-pair space — the species tree's
leaves, the molecule tree's leaves, then the field — cost more than one eight-state enumeration.
So the two costs compound rather than trade:

| | reference | current |
| --- | ---: | ---: |
| seconds per chain | 196.9 | 246.3 |
| effective field draws per second | 1.615 | 0.276 |

**5.9 times fewer effective draws of the field per second of wall clock**, which is the 4.7 in
mixing multiplied by the 1.25 in iteration cost.

## The same measurement at `omega = 0.1`

Two replicates per binary, so only the stationary factors carry a claim.

| factor | reference `tau` | current `tau` | ratio | reference drift | current drift |
| --- | ---: | ---: | ---: | ---: | ---: |
| `data` — the field's own | 71.4 | 88.0 | **1.23** | −0.08 | 0.06 |
| `link` | 864 | 817 | 0.95 | −0.26 | −0.48 |
| `species_node_state` | 1389 | 1525 | 1.10 | 0.64 | −0.87 |
| `molecules_node_state` | 2862 | 909 | 0.32 | 0.11 | 0.28 |
| `joint_density` — the total | 2893 | 1609 | 0.56 | 0.07 | −0.71 |

The field's cost has all but gone: `acf(1)` is 0.70 against the reference's 0.64, and the two
autocorrelation functions are within noise of each other by lag 10. Per second of wall clock the
current binary gives 1.57 times fewer effective field draws, which is the 1.27x iteration cost here
multiplied by the 1.23x in mixing — the mixing has shrunk to a quarter of what it was at
`omega = 0.005`, but it has not vanished.

**Do not read the other rows.** `molecules_node_state` at 0.32 is what drags the total to 0.56, and
it is two chains per binary on a factor neither had converged. Those rows are noise, and no claim
rests on them.

**This is the regime dependence ADR-0005 predicted, confirmed.** The metastability is a
small-`omega` effect: it is the deterministic AND at `omega = 0` that pins the triple, and the trap
loosens as the link is allowed to be wrong more often. The block update was worth its complexity
exactly where `omega` is small.

## Both tree field posteriors

The field's own posterior cannot see the failure that shows up first. ADR-0005, derivation 3: the
field's density is the product `a~_s * a~_m`, so `omega` and the two alphas trade at constant
product, and a chain riding that ridge moves the two tree fields in opposite directions while the
field stays put. So each run is a point on the ridge, in a coordinate **along** it
(`log a~_s − log a~_m`, which moves when the trees trade) and one **across** it
(`log a~_s + log a~_m`, which moves only when the field's density does).

At `omega = 0.005`, replicates pooled:

| tree | reference density | current density | mean abs diff | max abs diff | correlation |
| --- | ---: | ---: | ---: | ---: | ---: |
| species | 0.50414 | 0.50412 | 0.00404 | 0.06955 | 0.9998 |
| molecules | 0.49651 | 0.49692 | 0.00414 | 0.08615 | 0.9998 |

Current against reference, averaged over seeds: **−0.00087 along the ridge and +0.00079 across
it.**

- **Along** — this is the coordinate the metastability would move, and the shift is well inside
  the scatter of either binary's own replicates (0.0025 for the reference, 0.0043 for the current).
  Neither binary rides the ridge, and they do not sit at different places on it.
- **Across** — the shift is about three times the reference's own replicate spread (0.00026) and
  about equal to the current binary's (0.00082), so it is not simply noise. In absolute terms it is
  a 0.08% difference in the field's density, which is far below anything a posterior would be read
  at, but it is a difference and not a zero.

**The two agree about where the tree fields are.** The measured cost is in how fast the field is
explored, and — to within a tenth of a percent on the field's density — not in where the chain ends
up.

## What this measurement does not settle

- **The phylogenetic parameters had not converged**, in either binary, at 20 000 iterations on this
  fixture. The two node-state factors drift by about a standard deviation across the trace. The
  leaf-layer question is answered anyway, because the factor that answers it is stationary — but
  the tree field posteriors above are compared at a point the chain was still moving through. They
  agree there; that they would still agree after the tree parameters settle is an inference, not a
  measurement. Settling it costs a chain long enough to converge 380 branch lengths.
- **One fixture, one data set, one link policy.** The `omega` dependence is measured at two points,
  which is enough to say the effect has the shape ADR-0005 predicted and not enough to say where
  the crossover is.
- **`tau` is estimated, not known.** At `omega = 0.005` the current binary's field factor leaves an
  effective sample of about 68 per chain, which makes its own `tau` a noisy number. The four
  replicates do not overlap the reference's four, so the 4.7 is a ratio to trust as an order of
  magnitude and not to three digits. At `omega = 0.1` there are only two replicates per binary,
  which is why only the two stationary factors are read there.

## If the cost is judged unacceptable

It may well be, for a production run at a small `omega`: 5.9 times fewer effective field draws per
second is a real bill, and it is paid in the one variable the whole model exists to infer. Four
things could be done about it, in rising order of what they cost to build.

1. **Do nothing, and run longer.** The posteriors agree, so the chain is not wrong, it is slow.
   Where `omega` is not small this is already close to the right answer — at 0.1 what is left is
   mostly the iteration cost.
2. **Recover the iteration time.** A quarter of the loss is not mixing at all, it is three passes
   over the leaf-pair space where there was one. Fusing the two trees' leaf blocks and the field
   into a single pass would keep the single-variable draws and give the 1.25x back.
3. **Put a block move back as a proposal, not as the update.** The single-variable walk stays the
   update, and a Metropolis move that flips the triple `(Y, Z_s, Z_m)` at one leaf pair together is
   added beside it. That is far less machinery than the deleted block update — no eight-state
   enumeration, no third traversal, no model seam — and it targets exactly the trap.
4. **Reinstate the block update.** ADR-0005's own answer. It is the most complexity for the most
   mixing, and this measurement says what it buys: a factor of about 5 in the field's
   autocorrelation time at `omega = 0.005`, and about 1.2 at `omega = 0.1`.

The trade recorded in issue #57 was simplicity for mixing, deliberately. This measurement prices
it, and the price is bounded: it is confined to the field, it is a factor of about 5 at the
pessimistic end of `omega`, it shrinks to about 1.2 as `omega` grows, and it does not move the
answer.

## Reproducing

```bash
bash model_validation/mixing_cost/run.sh                          # omega = 0.005, 4 replicates
ACOL_MIXING_OMEGA=0.1 ACOL_MIXING_REPLICATES=2 \
    ACOL_MIXING_DIR="$PWD/build/mixing_cost/run_omega_0.1" \
    bash model_validation/mixing_cost/run.sh                      # omega = 0.1, 2 replicates
```

Each writes `mixing_cost.json` beside its runs, which is where every number above comes from.
