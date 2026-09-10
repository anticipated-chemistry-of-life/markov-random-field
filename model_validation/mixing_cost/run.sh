#!/usr/bin/env bash
#
# The mixing cost of dropping the block update, measured rather than assumed.
#
# ADR-0005 built the leaf layer's eight-state block update on an argument: at a small error
# probability a field cell at one pins both tree fields, the triple is metastable under single-site
# draws, and the failure would present as slow mixing rather than as a bug. Issue #68 deleted that
# block update and draws the leaf layer with its tree instead. Slow mixing does not announce itself,
# so this looks for it on purpose.
#
# Two binaries, one data set:
#
#   reference   the revision before the block update was deleted, built from a git worktree
#   current     the working tree
#
# WHY THE DATA IS FIXED AND THE SEED IS NOT
#
# The two binaries cannot be compared seed for seed. A leaf's cell uniform used to come from the
# field stream and now comes from the node-state stream, at the cell's own linear index in
# node-state container space; the field then draws its own pass. So the same seed drives two
# different sequences of draws, and a seed-matched pair of runs would differ for a reason that says
# nothing about mixing. The default reference already carries the hash-map storages, so against it
# the draw pattern is the whole of the difference -- findings.md says so at more length.
#
# The data set is therefore the control variable: it is simulated once, by the reference binary,
# and both binaries infer from those same bytes. The seed is deliberately let go -- several chains
# per binary, and every comparison is between the two binaries' spreads and medians rather than
# between one chain and one other chain.
#
# Both binaries link the same coretools and stattools checkout, the one at the root of this
# repository. FETCHCONTENT_SOURCE_DIR_* is what pins that, and it matters: a reference binary built
# against a different dependency revision would put the difference somewhere the measurement cannot
# see.
#
# Usage:  bash model_validation/mixing_cost/run.sh
#
# Environment:
#   ACOL_MODE                    debug | release          (default release)
#   ACOL_ENV                     micromamba environment   (default acol_env)
#   MAMBA_EXE                    path to micromamba       (default: the one on PATH)
#   ACOL_MIXING_DIR              where to run             (default build/mixing_cost/run)
#   ACOL_MIXING_REFERENCE        the pre-refactor revision (default a3d1ae6, the parent of
#                                "Draw the leaf layer with its tree, and the field on its own")
#   ACOL_MIXING_FIXTURE          directory holding the four tree and paper-count files
#                                (default model_validation/s_balanced_255_m_balanced_511)
#   ACOL_MIXING_ITERATIONS       chain length             (default 20000)
#   ACOL_MIXING_BURNIN           burn-in length           (default 2000)
#   ACOL_MIXING_REPLICATES       chains per binary        (default 4)
#   ACOL_MIXING_SIMULATE_SEED    the one seed that is fixed (default 20260909)
#   ACOL_MIXING_OMEGA            the simulated error probability (default 0.005, the model's own
#                                default and the regime ADR-0005 warns about)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

MODE="${ACOL_MODE:-release}"
CONDA_ENV="${ACOL_ENV:-acol_env}"
MAMBA="${MAMBA_EXE:-micromamba}"
REFERENCE_REV="${ACOL_MIXING_REFERENCE:-a3d1ae6}"
FIXTURE_DIR="${ACOL_MIXING_FIXTURE:-$ROOT/model_validation/s_balanced_255_m_balanced_511}"
ITERATIONS="${ACOL_MIXING_ITERATIONS:-20000}"
BURNIN="${ACOL_MIXING_BURNIN:-2000}"
NUM_BURNIN=1
REPLICATES="${ACOL_MIXING_REPLICATES:-4}"
SIMULATE_SEED="${ACOL_MIXING_SIMULATE_SEED:-20260909}"
OMEGA="${ACOL_MIXING_OMEGA:-0.005}"
WORKDIR="${ACOL_MIXING_DIR:-$ROOT/build/mixing_cost/run}"
WORKTREE="$ROOT/build/mixing_cost/reference-src"

# Both data sources, as everywhere else in this repository: LOTUS and the simple error model are
# the two things that read the field, and a measurement of the field's mixing wants both of them
# pulling on it.
FLAGS="ls"

# The directory is wiped below, so a stray or empty one would take something else with it.
case "$WORKDIR" in
    /*/?*) ;;
    *) echo "error: ACOL_MIXING_DIR must be an absolute path, but is '$WORKDIR'" >&2; exit 1 ;;
esac
for name in ITERATIONS REPLICATES; do
    value="${!name}"
    if ! [[ "$value" =~ ^[0-9]+$ ]] || ((10#$value == 0)); then
        echo "error: ACOL_MIXING_$name must be a positive integer, but is '$value'" >&2
        exit 1
    fi
done
# Zero is a legitimate burn-in: the comparison drops the rows it names, and naming none is allowed.
for name in BURNIN SIMULATE_SEED; do
    value="${!name}"
    if ! [[ "$value" =~ ^[0-9]+$ ]]; then
        echo "error: ACOL_MIXING_$name must be a non-negative integer, but is '$value'" >&2
        exit 1
    fi
done
# Checked here rather than left to the binary, because it is also written into meta.json as JSON.
if ! [[ "$OMEGA" =~ ^0*\.[0-9]+$ ]]; then
    echo "error: ACOL_MIXING_OMEGA must be a decimal in (0, 0.5), but is '$OMEGA'" >&2
    exit 1
fi

FIXTURE=(species.txt molecules.txt species_papers.txt molecules_papers.txt)
for file in "${FIXTURE[@]}"; do
    [[ -f "$FIXTURE_DIR/$file" ]] || {
        echo "error: the fixture directory '$FIXTURE_DIR' has no $file" >&2; exit 1
    }
done

command -v "$MAMBA" >/dev/null 2>&1 || {
    echo "error: micromamba not found (set MAMBA_EXE to its path)" >&2; exit 1
}
"$MAMBA" run -n "$CONDA_ENV" true >/dev/null 2>&1 || {
    echo "error: micromamba environment '$CONDA_ENV' is missing or broken; run 'just setup'" >&2
    exit 1
}

# ---------------------------------------------------------------------------
# The two binaries
# ---------------------------------------------------------------------------

cd "$ROOT"

echo "==> building the current binary"
just build "$MODE" "$FLAGS" >/dev/null
CURRENT_BIN="$ROOT/$(just bin "$MODE" "$FLAGS")"
CURRENT_REV="$(git rev-parse --short HEAD)"

REFERENCE_REV="$(git rev-parse --short "$REFERENCE_REV")"
echo "==> building the reference binary at $REFERENCE_REV"
# A worktree rather than a checkout, so the working tree the current binary was built from is left
# alone. Detached, because the revision is a point in history and not a branch to move.
if [[ -d "$WORKTREE" ]]; then
    git -C "$WORKTREE" checkout --detach --force "$REFERENCE_REV" >/dev/null || {
        echo "error: '$WORKTREE' exists but is not a worktree this can move to $REFERENCE_REV." >&2
        echo "       Delete it and run again." >&2
        exit 1
    }
else
    git worktree prune
    git worktree add --detach "$WORKTREE" "$REFERENCE_REV" >/dev/null
fi

ACOL_FLAG_SUFFIX="-$FLAGS" "$MAMBA" run -n "$CONDA_ENV" bash -eu -c '
    if [[ -z "${CXX:-}" ]]; then
        case "$(uname -s)" in
            Darwin) export CC="$CONDA_PREFIX/bin/clang" CXX="$CONDA_PREFIX/bin/clang++" ;;
            *)      export CC="$CONDA_PREFIX/bin/gcc"   CXX="$CONDA_PREFIX/bin/g++" ;;
        esac
    fi
    cd "$2"
    # BUILD_TESTING is off because only the binary is wanted here. The test suite of the
    # reference revision is already recorded in its own history.
    cmake --preset "$1" -DLOTUS=ON -DSIMPLE_DATA=ON -DUSE_MS_DATA=OFF -DBUILD_TESTING=OFF \
          -DFETCHCONTENT_SOURCE_DIR_CORETOOLS="$3/coretools" \
          -DFETCHCONTENT_SOURCE_DIR_STATTOOLS="$3/stattools" >/dev/null
    exec cmake --build "$4" --target acol
' _ "$MODE" "$WORKTREE" "$ROOT" "build/$MODE-$FLAGS" | tail -1
REFERENCE_BIN="$WORKTREE/build/$MODE-$FLAGS/acol"

for binary in "$CURRENT_BIN" "$REFERENCE_BIN"; do
    [[ -x "$binary" ]] || { echo "error: $binary was not built" >&2; exit 1; }
done

# ---------------------------------------------------------------------------
# One data set, simulated by the reference binary
# ---------------------------------------------------------------------------

rm -rf "$WORKDIR"
mkdir -p "$WORKDIR/simulate"
cp "${FIXTURE[@]/#/$FIXTURE_DIR/}" "$WORKDIR/"

echo "==> simulating one data set with the reference binary, seed $SIMULATE_SEED"
(cd "$WORKDIR" && "$REFERENCE_BIN" simulate \
    --out simulate/acol \
    --tree_species species.txt --tree_molecules molecules.txt \
    --species_paper_counts species_papers.txt \
    --molecules_paper_counts molecules_papers.txt \
    --iterations 1 \
    --error_probability "$OMEGA" --epsilon_simple_model 0.1 --gamma 1.1 \
    --numThreads 1 --fixedSeed "$SIMULATE_SEED" \
    --write_Y --write_Z --write_joint_log_prob_density >/dev/null)

# ---------------------------------------------------------------------------
# Both binaries, on those same bytes
#
# One thread, because an infer run is only reproducible at one. Every cell draw is hashed from the
# cell's position (ADR-0007), but the alpha and nu moves still draw from coretools' thread-local
# generator inside a dynamic loop over cliques, so at more than one thread the same seed gives a
# different chain each time. A measurement of mixing wants the chain the seed names.
# ---------------------------------------------------------------------------

timings=""
run_chain() {
    local label="$1" binary="$2" seed="$3"
    local out="$WORKDIR/$label/seed_$seed"
    mkdir -p "$out"
    local started
    started="$SECONDS"
    (cd "$WORKDIR" && "$binary" infer \
        --out "$label/seed_$seed/acol" \
        --tree_species species.txt --tree_molecules molecules.txt \
        --species_paper_counts species_papers.txt \
        --molecules_paper_counts molecules_papers.txt \
        --lotus simulate/acol_simulated_lotus.tsv \
        --simple_data simulate/acol_simulated_simple_data.tsv \
        --iterations "$ITERATIONS" --burnin "$BURNIN" --numBurnin "$NUM_BURNIN" \
        --epsilon_simple_model 0.1 --gamma 1.1 \
        --numThreads 1 --fixedSeed "$seed" \
        --write_joint_log_prob_density --write_tree_field_posteriors >/dev/null)
    # Whole seconds: a chain worth measuring runs for minutes, and this keeps the run free of an
    # interpreter outside the environment everything else uses.
    timings="$timings$label $seed $((SECONDS - started))
"
}

for seed in $(seq 1 "$REPLICATES"); do
    echo "==> reference, seed $seed"
    run_chain reference "$REFERENCE_BIN" "$seed"
    echo "==> current, seed $seed"
    run_chain current "$CURRENT_BIN" "$seed"
done

# ---------------------------------------------------------------------------
# What the comparison needs to know about the runs
# ---------------------------------------------------------------------------

# A leaf is a node no row names as a parent. The field's container space is the product of the two
# leaf counts, which is what every posterior mean below divides by.
count_leaves() {
    awk -F'\t' 'NR > 1 { child[$1]; parent[$2] }
                END { n = 0; for (c in child) if (!(c in parent)) n++; print n }' "$1"
}
N_SPECIES_LEAVES="$(count_leaves "$WORKDIR/species.txt")"
N_MOLECULES_LEAVES="$(count_leaves "$WORKDIR/molecules.txt")"

# The joint density trace is written from the first iteration on: unlike the parameter traces,
# which stattools withholds until --writeBurnin says otherwise, nothing in that writer knows the
# chain has not started yet. So the comparison is told how many rows to drop.
BURN_IN_ROWS=$((BURNIN * NUM_BURNIN))

# The quoted heredoc delimiter is the point: every value below arrives through the environment, so
# a fixture path with a quote in it is a string and never a line of source.
ACOL_META_TIMINGS="$timings" \
ACOL_META_REFERENCE="$REFERENCE_REV" ACOL_META_CURRENT="$CURRENT_REV" \
ACOL_META_FIXTURE="$FIXTURE_DIR" \
ACOL_META_SPECIES_LEAVES="$N_SPECIES_LEAVES" ACOL_META_MOLECULES_LEAVES="$N_MOLECULES_LEAVES" \
ACOL_META_ITERATIONS="$ITERATIONS" ACOL_META_BURNIN="$BURNIN" ACOL_META_NUM_BURNIN="$NUM_BURNIN" \
ACOL_META_BURN_IN_ROWS="$BURN_IN_ROWS" ACOL_META_REPLICATES="$REPLICATES" \
ACOL_META_SIMULATE_SEED="$SIMULATE_SEED" ACOL_META_OMEGA="$OMEGA" \
python3 - "$WORKDIR/meta.json" <<'PYTHON'
import json, os, sys

timings = {}
for line in os.environ["ACOL_META_TIMINGS"].strip().splitlines():
    label, seed, seconds = line.split()
    timings.setdefault(label, {})[seed] = float(seconds)

n_species = int(os.environ["ACOL_META_SPECIES_LEAVES"])
n_molecules = int(os.environ["ACOL_META_MOLECULES_LEAVES"])
n_replicates = int(os.environ["ACOL_META_REPLICATES"])

meta = {
    "reference_revision": os.environ["ACOL_META_REFERENCE"],
    "current_revision": os.environ["ACOL_META_CURRENT"],
    "fixture": os.environ["ACOL_META_FIXTURE"],
    "prefix": "acol",
    "tree_names": ["species", "molecules"],
    "n_species_leaves": n_species,
    "n_molecules_leaves": n_molecules,
    "n_cells": n_species * n_molecules,
    "iterations": int(os.environ["ACOL_META_ITERATIONS"]),
    "burnin": int(os.environ["ACOL_META_BURNIN"]),
    "num_burnin": int(os.environ["ACOL_META_NUM_BURNIN"]),
    "burn_in_rows": int(os.environ["ACOL_META_BURN_IN_ROWS"]),
    "n_replicates": n_replicates,
    "seeds": list(range(1, n_replicates + 1)),
    "simulate_seed": int(os.environ["ACOL_META_SIMULATE_SEED"]),
    "error_probability": float(os.environ["ACOL_META_OMEGA"]),
    "seconds": timings,
}
with open(sys.argv[1], "w") as handle:
    json.dump(meta, handle, indent=2, sort_keys=True)
    handle.write("\n")
PYTHON

echo
echo "==> comparing"
cd "$ROOT/model_validation"
uv run python compare_mixing.py "$WORKDIR" --json "$WORKDIR/mixing_cost.json"
