#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$SCRIPT_DIR/.."

MODE="${ACOL_MODE:-release}"
FLAGS="${ACOL_DATA_FLAGS:-l}"
MAMBA="${MAMBA_EXE:-micromamba}"
CONDA_ENV="${ACOL_ENV:-acol_env}"

cd "$ROOT"
just build "$MODE" "$FLAGS"
ACOL="$ROOT/$(just bin "$MODE" "$FLAGS")"

cd "$SCRIPT_DIR"
# The binary is linked against the environment's libstdc++ and carries no rpath, so it has to be
# launched inside the environment or the loader takes the system one.
#
# `samply record --save-only` writes ./profile.json.gz and opens nothing. The profile is only as
# useful as the symbols it is read with, so do not rebuild the binary before reading it: samply
# stores addresses, and `just build` relinks them somewhere else.
"$MAMBA" run -n "$CONDA_ENV" "$ACOL" infer \
    --out ./test_out/acol \
    --tree_species species.tsv \
    --tree_molecules molecules.tsv \
    --lotus lotus.tsv \
    --species_paper_counts species_paper_counts.tsv \
    --molecules_paper_counts molecules_paper_counts.tsv \
    --iterations 5000 \
    --numThreads 32 \
    --write_joint_log_prob_density
