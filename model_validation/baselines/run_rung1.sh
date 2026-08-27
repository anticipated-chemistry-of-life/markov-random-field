#!/usr/bin/env bash
# Runs the rung-1 pinned scenario into $SCENARIO/rung1_gate, single-threaded.
#
# Single-threaded is not a preference, it is what makes the run reproducible.
# coretools' TRandomGenerator is `static thread_local`, so every thread owns its
# own stream; the sweep over cliques is `schedule(dynamic)`, so which clique is
# drawn by which thread varies between runs. With --numThreads all the same seed
# therefore yields a different chain each time -- verified empirically, the whole
# species trace diverges. With --numThreads 1 the run is byte-stable.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$SCRIPT_DIR/../.."
SCENARIO="$ROOT/model_validation/independent_y_s255_m255_seed42"

MODE="${ACOL_MODE:-release}"
FLAGS="s"

[ -d "$SCENARIO" ] || {
    echo "error: scenario '$SCENARIO' is missing (it is gitignored)." >&2
    echo "       regenerate it: cd model_validation && uv run python simulate_independent.py --seed 42" >&2
    exit 1
}

cd "$ROOT"
just build "$MODE" "$FLAGS"
ACOL="$ROOT/$(just bin "$MODE" "$FLAGS")"

cd "$SCENARIO"
rm -rf rung1_gate && mkdir -p rung1_gate
"$ACOL" infer \
    --out ./rung1_gate/acol \
    --tree_species species.txt \
    --tree_molecules molecules.txt \
    --species_paper_counts species_papers.txt \
    --molecules_paper_counts molecules_papers.txt \
    --lotus simulated_lotus.tsv \
    --simple_data simulated_simple_data.tsv \
    --epsilon_simple_model 0.05 \
    --gamma 1.1 \
    --epsilon 0.001 \
    --iterations 10000 \
    --burnin 1000 \
    --numBurnin 10 \
    --fixedSeed 42 \
    --numThreads 1 \
    --writeBurnin \
    --write_joint_log_prob_density \
    --molecules_branch_lengths simulated_pinned_molecules.txt --molecules_branch_lengths.update false \
    --molecules_mean_log_nu simulated_pinned_molecules.txt --molecules_mean_log_nu.update false \
    --molecules_var_log_nu simulated_pinned_molecules.txt --molecules_var_log_nu.update false \
    --molecules_log_nu simulated_pinned_molecules.txt --molecules_log_nu.update false \
    --molecules_alpha simulated_pinned_molecules.txt --molecules_alpha.update false \
    --set_Y simulated_Y.txt --Y.update false \
    --set_species_Z simulated_Z_species.txt --set_molecules_Z simulated_Z_molecules.txt --Z.update false
