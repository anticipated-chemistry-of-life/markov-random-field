# Phase-1 acceptance baseline

Phase 1 of the TTree split (issue #10) is defined as behaviour-preserving: every
ticket in it must leave the model's output byte-identical. This directory holds
the baseline that claim is measured against.

## Usage

```bash
./verify_rung1.sh          # rebuild, rerun, compare against the manifest
```

Run it after every commit in phase 1, not once at the end. Nine changes verified
together tell you that something broke; nine verified one at a time tell you
which.

`--record` rewrites the manifest. It is for establishing a new baseline, and
must not be used during phase 1 -- the point is that the manifest does not move.

## What is committed, and what is not

The scenario directory `independent_y_s255_m255_seed42/` is gitignored on purpose
(`model_validation/.gitignore`), because it is generated. So the ~8.6 MB of run
output cannot be committed, and `*_trace.txt` is ignored repo-wide besides.

What is committed is `rung1_pre_split.sha256`, a SHA-256 per output file. For a
byte-identical gate a hash is exactly as strong as the file it stands for, and it
is about a kilobyte.

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
this refactor.

## Scope of the gate

The build is `release` with flags `s` -- the simple error model only. LOTUS is
compiled out, so anything guarded by `USE_LOTUS` is invisible to this gate.
