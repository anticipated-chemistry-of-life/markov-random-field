#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$SCRIPT_DIR/.."

MODE="${ACOL_MODE:-release}"
FLAGS="${ACOL_DATA_FLAGS:-lm}"
# The unknown rate rho: the unknown molecule's weight per structure proposed for a feature
# (ADR-0018). 1.4 is 122, the unknown weight ADR-0018 reads off the held-out odds, over the mean
# of about 87 rows per feature there. It is an input, not an estimate.
MS_UNKNOWN_RATE="${ACOL_MS_UNKNOWN_RATE:-1.4}"

cd "$ROOT"
pixi run build "$MODE-$FLAGS"
ACOL="$ROOT/build/$MODE-$FLAGS/acol"

cd "$SCRIPT_DIR"
# The binary is linked against the environment's libstdc++ and carries no rpath, so it has to be
# launched inside the environment or the loader takes the system one.
#
# `samply record --save-only` writes ./profile.json.gz and opens nothing. The profile is only as
# useful as the symbols it is read with, so do not rebuild the binary before reading it: samply
# stores addresses, and `pixi run build` relinks them somewhere else.
pixi run "$ACOL" infer \
    --out ./test_out/acol \
    --tree_species species.tsv \
    --tree_molecules molecules.tsv \
    --lotus lotus.tsv \
    --species_paper_counts species_paper_counts.tsv \
    --molecules_paper_counts molecules_paper_counts.tsv \
    --iterations 5000 \
    --numThreads 32 \
    --ms_data ms_data.tsv \
    --ms_unknown_rate "$MS_UNKNOWN_RATE"
