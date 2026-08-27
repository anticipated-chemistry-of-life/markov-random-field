# Acceptance baseline for the rung-1 run

Some changes to the model are supposed to leave its output alone, bit for bit.
This directory holds the baseline that claim is measured against: a SHA-256 per
output file of the rung-1 pinned run.

```bash
./verify_rung1.sh          # rebuild, rerun, compare against the manifest
./verify_rung1.sh --record # establish a new baseline
```

Run it after every commit in a behaviour-preserving stretch of work, not once at
the end. Nine changes verified together tell you that something broke; nine
verified one at a time tell you which.

`--record` is for the commit that deliberately moves the output, and for nothing
else -- during a behaviour-preserving stretch the whole point is that the
manifest does not move.

## What the manifests are

- **`rung1_canonical_order.sha256`** is the live one, and what `verify_rung1.sh`
  checks. It was recorded at the commit that put nodes in canonical order
  (ADR-0004), which is the last change that was *meant* to move the output.
- **`rung1_pre_split.sha256`** is the phase-1 record: the `TTree` split (#10) was
  defined as behaviour-preserving and every one of its nine commits was checked
  against this. It is kept as the evidence for that claim and is no longer
  checked by anything. Reproducing it needs the pre-split node ordering on
  *both* sides -- the binary and `src/independent/indexing.py`, which generates
  the scenario -- so it cannot be re-run against today's tree.

That last point generalises: **a manifest is a statement about a binary and a
scenario together.** The scenario is generated, so a change to the harness that
generates it invalidates the manifest exactly as surely as a change to the model.

## What is committed, and what is not

The scenario directory `independent_y_s255_m255_seed42/` is gitignored on purpose
(`model_validation/.gitignore`), because it is generated. So the ~8.6 MB of run
output cannot be committed, and `*_trace.txt` is ignored repo-wide besides.

What is committed is the manifest. For a byte-identical gate a hash is exactly as
strong as the file it stands for, and it is about a kilobyte.

If the scenario directory is missing, regenerate it first:

```bash
cd model_validation && uv run python simulate_independent.py --seed 42
```

## Why the run is single-threaded

`--numThreads 1` is not a preference. coretools' `TRandomGenerator` is
`static thread_local`, so each thread owns its own stream, and the sweep over
cliques is `schedule(dynamic)` -- so which clique is drawn by which thread varies
between runs. Under `--numThreads all` the same `--fixedSeed 42` therefore
produces a different chain every time; verified empirically, the entire species
trace diverges from line 3 onwards, while the pinned molecules side stays put.

Single-threaded, two runs agree on every output file. `acol.log` is excluded from
the manifest regardless: it carries a per-run ntfy topic UUID and wall-clock
timings, neither of which say anything about the model.

A consequence worth knowing outside this directory: **a multi-threaded acol run
is not reproducible from its seed.** That is a property of the harness, not of
any one refactor.

## Scope of the gate

The build is `release` with flags `s` -- the simple error model only. LOTUS is
compiled out, so anything guarded by `USE_LOTUS` is invisible to this gate. Only
`infer` is run, so the simulation path is invisible to it too.
