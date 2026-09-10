# The mixing cost of dropping the block update

ADR-0005 built the leaf layer's eight-state block update on an argument rather than a measurement:
at a small error probability a field cell at one pins both tree fields, the triple is metastable
under single-variable draws, and the failure would present as **slow mixing rather than as a bug**.
Issue #68 deleted that block update. A cost that only shows up as slow mixing does not announce
itself, so this looks for it on purpose.

```bash
bash model_validation/mixing_cost/run.sh
```

`run.sh` builds two binaries, simulates one data set, runs several chains under each and hands the
result to `../compare_mixing.py`. `findings.md` is the write-up of the run that was recorded, in
the form the superseding architecture decision record carries.

## What is compared

| binary      | revision                                                                |
| ----------- | ----------------------------------------------------------------------- |
| `reference` | the commit before the block update was deleted, built from a git worktree |
| `current`   | the working tree                                                        |

Both link the same `coretools` and `stattools` checkout — the one at the root of this repository,
pinned with `FETCHCONTENT_SOURCE_DIR_*`. A reference binary built against a different dependency
revision would put the difference somewhere the measurement cannot see.

**Since ADR-0010 the default pair is source-identical.** The block update was reinstated by reverting
the commit that deleted it, so the working tree now matches `a3d1ae6` in `src/` and `tests/`, and
`current` is always the working tree. A default run therefore compares two identical binaries and
measures nothing.

Two comparisons are worth naming. To measure the block against the single-variable walk **from
today's tree**, set `ACOL_MIXING_REFERENCE=f13fef9`: the reference is then the revision that drew
the leaf layer one variable at a time, and the ratios read the other way round from `findings.md`.
To **reproduce the recorded table**, run from a tree at `f13fef9` with the default reference, which
is the pair that produced it.

## Why the data is fixed and the seed is not

**The two binaries cannot be compared seed for seed.** A leaf's cell uniform used to come from the
field stream, drawn by the block update; it now comes from the node-state stream, at the leaf's own
linear index in node-state container space, and the field then draws its own pass. So one seed
drives two different sequences of draws, and a seed-matched pair of runs would differ for a reason
that says nothing at all about mixing.

Issue #69 names two reasons, the storage layer and the draw pattern. Against the default reference
only the second applies: `a3d1ae6` already carries the hash-map storages, which landed three
commits earlier, and `git diff a3d1ae6 f13fef9 -- src/storages/` is empty. The draw pattern alone
is enough. `findings.md` records this.

The **data set** is therefore the control variable. It is simulated once, by the reference binary,
and both binaries infer from those same bytes. The **seed** is deliberately let go: several chains
per binary, and every comparison is between the two binaries' medians and spreads rather than
between one chain and one other chain.

## The two instruments

**The joint density trace** says how fast a chain forgets where it was. Every factor of the
ADR-0005 factorisation is a proper density, so the total moves only with the configuration and the
parameters — it sees the leaf layer along with everything else. Its integrated autocorrelation time
`tau` is what one independent draw costs, and `ESS/s` is what that costs in wall-clock time. The
second number matters as much as the first: the block update buys mixing with arithmetic, so a
`tau` that doubles while the iteration halves in cost is not a loss.

**Both tree field posteriors** say where the chain went, and the field's own posterior cannot tell
you. ADR-0005, derivation 3: the field's density is the product `a~_s * a~_m` of the two corrupted
tree field rates, so `omega` and the two alphas trade against each other at constant product. A
chain riding that ridge moves the two tree fields in opposite directions and leaves the field
exactly where it was. So each run is reported as a point on the ridge, in two coordinates:

- **along** — `log a~_s - log a~_m`, which moves when the two trees trade against each other;
- **across** — `log a~_s + log a~_m`, which moves only when the field's own density does.

A binary whose replicates scatter **along** the ridge while the product stays put is riding it.
That is the failure the block update existed to prevent, and it is what shows up first.

## Knobs

`run.sh` documents them in its header. The ones worth reaching for:

| variable                | default                                    | why you would change it                         |
| ----------------------- | ------------------------------------------ | ----------------------------------------------- |
| `ACOL_MIXING_OMEGA`     | `0.005`                                    | the metastability is a small-`omega` effect, so this is the knob that sets how hard the measurement is |
| `ACOL_MIXING_ITERATIONS`| `20000`                                    | a `tau` of `t` needs a chain many multiples of `t` long before it is estimated rather than guessed |
| `ACOL_MIXING_REPLICATES`| `4`                                        | the spread along the ridge is a statistic over these |
| `ACOL_MIXING_FIXTURE`   | `../s_balanced_255_m_balanced_511`         | 128 x 256 leaf pairs; a smaller one runs faster and says less |

Everything runs at `--numThreads 1`. An `infer` run is only reproducible at one thread: every cell
draw is hashed from the cell's position (ADR-0007), but the `alpha` and `nu` moves still draw from
coretools' thread-local generator inside a dynamic loop over cliques.
