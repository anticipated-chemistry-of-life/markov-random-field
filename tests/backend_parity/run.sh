#!/usr/bin/env bash
#
# The dense-versus-sparse gate: two binaries, byte-identical output.
#
# Which storage backs the field, and which backs the node state, are two aliases in
# src/storages/storage_backend.h. The comparison here is therefore between two *builds* rather than
# between two runtime modes, which means it exercises exactly what ships, with no dispatch layer
# that exists only for testing. Both binaries come from the same sources and differ only in the
# defines that override those two aliases.
#
# The gate pairs sparse with sparse and dense with dense. Why those two and not the mixed pairs is
# in ADR-0006, and src/storages/storage_backend.h records which pairs it gates.
#
# Each backend runs the same two chains from the same seed, one after the other, in the *same*
# working directory with the fixture copied in and every argument spelled identically -- the
# outputs are moved aside between the two. Reusing the directory is what makes even
# `acol.parameters` comparable: it echoes the command line with every path resolved against the
# working directory, so two runs in two directories would differ there and nowhere else. Then every
# file one run wrote is compared with the other's byte for byte:
#
#   simulate -> the field and both node states, in full, plus the LOTUS and simple-error data
#               drawn from them and the per-iteration traces.
#   infer    -> the parameter traces, the field and node-state traces, the joint density and
#               the posterior field.
#
# Only `*.log` is left out: it carries a fresh ntfy topic UUID and wall-clock timings, so it
# differs between two runs of the *same* binary.
#
# Usage:  bash tests/backend_parity/run.sh          (or: just parity)
#
# Environment:
#   ACOL_MODE        debug | release          (default release)
#   ACOL_ENV         micromamba environment   (default acol_env)
#   MAMBA_EXE        path to micromamba       (default: the one on PATH)
#   ACOL_PARITY_DIR  where to run             (default build/parity)
#   ACOL_PARITY_SEED fixed seed               (default 42)
#   ACOL_PARITY_ITERATIONS  chain length      (default 400)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

MODE="${ACOL_MODE:-release}"
CONDA_ENV="${ACOL_ENV:-acol_env}"
MAMBA="${MAMBA_EXE:-micromamba}"
SEED="${ACOL_PARITY_SEED:-42}"
ITERATIONS="${ACOL_PARITY_ITERATIONS:-400}"
WORKDIR="${ACOL_PARITY_DIR:-$ROOT/build/parity}"

# Both data sources, and not a knob: LOTUS and the simple error model are the two things that read
# the field through the storage, so compiling in either alone would leave one of those read paths
# out of the comparison entirely.
FLAGS="ls"

if ! [[ "$ITERATIONS" =~ ^[0-9]+$ ]] || ((10#$ITERATIONS == 0)); then
    echo "error: ACOL_PARITY_ITERATIONS must be a positive integer, but is '$ITERATIONS'" >&2
    exit 1
fi
if ! [[ "$SEED" =~ ^[0-9]+$ ]]; then
    echo "error: ACOL_PARITY_SEED must be a non-negative integer, but is '$SEED'" >&2
    exit 1
fi
# The directory is wiped before the runs, so a stray or empty ACOL_PARITY_DIR would take something
# else with it.
case "$WORKDIR" in
    /*/?*) ;;
    *) echo "error: ACOL_PARITY_DIR must be an absolute path below the root, but is '$WORKDIR'" >&2
       exit 1 ;;
esac

# Above 32767 iterations the two fields thin their posterior counters differently -- the sparse
# counter is 15 bits and the dense one 16 -- and the thinning factor is also what decides which
# iterations get a trace line. The two would then write traces of different lengths, which is a
# property of the counter widths and not a regression. Keep the gate below that.
readonly MAX_IDENTICAL_THINNING_ITERATIONS=32767
if ((ITERATIONS > MAX_IDENTICAL_THINNING_ITERATIONS)); then
    echo "error: ACOL_PARITY_ITERATIONS=$ITERATIONS exceeds $MAX_IDENTICAL_THINNING_ITERATIONS," >&2
    echo "       above which the two backends thin their traces differently by design." >&2
    exit 1
fi

# Indexed rather than associative arrays, and indices rather than names throughout: macOS ships
# bash 3.2, which has no `declare -A`.
#
# One gated pair per index: its label, the field's storage and the node state's. `sparse` names the
# pair whose storages are both sparse, which is what the label meant before the two could differ.
BACKENDS=(sparse dense)
FIELD_STORAGES=(TStorageYMatrix TStorageYDense)
NODE_STATE_STORAGES=(TStorageZMatrix TStorageZDense)
BINARIES=()

# Copied into the working directory rather than referred to, so that both runs spell every argument
# identically -- see the note on `acol.parameters` above.
FIXTURE=(species.txt molecules.txt species_papers.txt molecules_papers.txt)

# ---------------------------------------------------------------------------
# Build one binary per gated pair
#
# The gate drives cmake itself rather than going through `just`. `just` knows nothing about the
# storages, because nothing in the build system chooses them any more, and this is the one build
# that overrides the aliases.
# ---------------------------------------------------------------------------

command -v "$MAMBA" >/dev/null 2>&1 || {
    echo "error: micromamba not found (set MAMBA_EXE to its path)" >&2
    exit 1
}
"$MAMBA" run -n "$CONDA_ENV" true >/dev/null 2>&1 || {
    echo "error: micromamba environment '$CONDA_ENV' is missing or broken; run 'just setup'" >&2
    exit 1
}

# One directory per pair, beside the ordinary ones, so a rerun rebuilds only what changed.
#
# Everything below runs inside the environment, because the flags are built from CXXFLAGS. The
# conda compiler packages export CC, CXX and the matching sysroot flags from their activation
# scripts. The defines are appended to CXXFLAGS rather than replacing it, so the sysroot flags
# survive. The plain compiler names are the fallback when nothing exported them.
#
# Configure only when the cache does not already hold exactly these flags. Re-running cmake
# regenerates armadillo's headers, which invalidates every object that includes them. Ninja
# re-runs cmake by itself when CMakeLists.txt or the presets change, so skipping it is safe.
build_binary() {
    local suffix="$1" field="$2" node_state="$3"
    local defines="-DACOL_FIELD_STORAGE=${field} -DACOL_NODE_STATE_STORAGE=${node_state}"

    ACOL_FLAG_SUFFIX="$suffix" "$MAMBA" run -n "$CONDA_ENV" bash -eu -c '
        if [[ -z "${CXX:-}" ]]; then
            case "$(uname -s)" in
                Darwin) export CC="$CONDA_PREFIX/bin/clang" CXX="$CONDA_PREFIX/bin/clang++" ;;
                *)      export CC="$CONDA_PREFIX/bin/gcc"   CXX="$CONDA_PREFIX/bin/g++" ;;
            esac
        fi
        flags="${CXXFLAGS:-} $3"
        if ! grep -qxF "CMAKE_CXX_FLAGS:STRING=$flags" "$2/CMakeCache.txt" 2>/dev/null; then
            cmake --preset "$1" -DLOTUS=ON -DSIMPLE_DATA=ON -DUSE_MS_DATA=OFF \
                  -DCMAKE_CXX_FLAGS="$flags"
        fi
        exec cmake --build "$2" --target acol
    ' _ "$MODE" "build/${MODE}${suffix}" "$defines"
}

cd "$ROOT"
for index in "${!BACKENDS[@]}"; do
    backend="${BACKENDS[$index]}"
    suffix="-${FLAGS}-parity-${backend}"
    echo "==> building the $backend-backed binary"
    build_binary "$suffix" "${FIELD_STORAGES[$index]}" "${NODE_STATE_STORAGES[$index]}"
    BINARIES+=("$ROOT/build/${MODE}${suffix}/acol")
done

# ---------------------------------------------------------------------------
# Run the same two chains under each backend
# ---------------------------------------------------------------------------

RUNDIR="$WORKDIR/run"
rm -rf "$WORKDIR"

# Both chains below pass `--numThreads 1`: the field update is parallel and the random stream is
# drawn from a shared generator, so a run is only reproducible at a fixed thread count (issue #38).
# Pinning it keeps the gate about the storage backend and nothing else. What that costs is the
# multi-batch commit of the deferred inserts, which one thread never produces -- that path is
# covered at the storage seam instead, by StorageEquivalence in tests/TStorageConformance_Tests.cpp.
run_acol() {
    local index="$1"; shift
    (cd "$RUNDIR" && "${BINARIES[$index]}" "$@" >/dev/null)
}

for index in "${!BACKENDS[@]}"; do
    backend="${BACKENDS[$index]}"

    rm -rf "$RUNDIR"
    mkdir -p "$RUNDIR/simulate" "$RUNDIR/infer"
    cp "${FIXTURE[@]/#/$SCRIPT_DIR/}" "$RUNDIR/"

    echo "==> $backend: simulate"
    run_acol "$index" simulate \
        --out simulate/acol \
        --tree_species species.txt --tree_molecules molecules.txt \
        --species_paper_counts species_papers.txt \
        --molecules_paper_counts molecules_papers.txt \
        --iterations "$ITERATIONS" --n_bins 6 \
        --epsilon_simple_model 0.1 --gamma 1.1 \
        --numThreads 1 --fixedSeed "$SEED" \
        --write_Y --write_Z --write_Y_trace --write_Z_trace

    # Inference reads the field and data this same run just simulated, so both backends infer from
    # bytes the simulate comparison below has already proven identical.
    echo "==> $backend: infer"
    run_acol "$index" infer \
        --out infer/acol \
        --tree_species species.txt --tree_molecules molecules.txt \
        --species_paper_counts species_papers.txt \
        --molecules_paper_counts molecules_papers.txt \
        --lotus simulate/acol_simulated_lotus.tsv \
        --simple_data simulate/acol_simulated_simple_data.tsv \
        --iterations "$ITERATIONS" --burnin 50 --numBurnin 2 --n_bins 6 \
        --epsilon_simple_model 0.1 --gamma 1.1 \
        --numThreads 1 --fixedSeed "$SEED" \
        --write_Y_trace --write_Z_trace --write_joint_log_prob_density

    # Everything compared below is written under the two --out prefixes, so anything that appears
    # beside them is an output the gate would not see. Today nothing does; `--write_branch_lengths`
    # is the one flag that would, because write_branch_length_grid names its file from a literal
    # rather than from the prefix -- which is why it is not passed above, and why this says so
    # rather than leaving the next person to find out.
    stray="$(comm -13 <(printf '%s\n' simulate infer "${FIXTURE[@]}" | sort) \
                       <(cd "$RUNDIR" && ls -1 | sort))"
    if [[ -n "$stray" ]]; then
        echo "FAIL: $backend wrote output outside the --out prefix, where the comparison below" >&2
        echo "      will not look at it:" >&2
        echo "$stray" | sed 's/^/        /' >&2
        exit 1
    fi

    mkdir -p "$WORKDIR/$backend"
    mv "$RUNDIR/simulate" "$RUNDIR/infer" "$WORKDIR/$backend/"
done
rm -rf "$RUNDIR"

compare_backends() {
    local subdir="$1"
    local left="$WORKDIR/${BACKENDS[0]}/$subdir"
    local right="$WORKDIR/${BACKENDS[1]}/$subdir"

    # A file one backend wrote and the other did not is a divergence in its own right, so the file
    # lists are compared before the contents.
    # `|| true` because `set -o pipefail` would otherwise make a directory of nothing but logs --
    # an empty grep -- abort the script before the "nothing to compare" check below can say so.
    local left_list right_list
    left_list="$(cd "$left" && ls -1 | grep -v '\.log$' | sort || true)"
    right_list="$(cd "$right" && ls -1 | grep -v '\.log$' | sort || true)"
    if [[ "$left_list" != "$right_list" ]]; then
        echo "FAIL: $subdir: the two backends wrote different files" >&2
        diff <(echo "$left_list") <(echo "$right_list") >&2 || true
        return 1
    fi

    local n_files
    n_files="$(printf '%s\n' "$left_list" | grep -c . || true)"
    if ((n_files == 0)); then
        echo "FAIL: $subdir: neither backend wrote anything to compare" >&2
        return 1
    fi

    local failed=0
    while IFS= read -r file; do
        [[ -n "$file" ]] || continue
        if ! cmp -s "$left/$file" "$right/$file"; then
            echo "FAIL: $subdir/$file differs between ${BACKENDS[0]} and ${BACKENDS[1]}" >&2
            diff "$left/$file" "$right/$file" | head -20 >&2 || true
            failed=1
        fi
    done <<<"$left_list"

    ((failed == 0)) || return 1
    echo "  $subdir: $n_files files identical"
}

# Both comparisons run even when the first fails, so a divergence is reported in full rather than
# one phase at a time.
divergences=0
echo "==> comparing the simulated field, node states and data"
compare_backends simulate || divergences=1
echo "==> comparing the parameter, field and node-state traces"
compare_backends infer || divergences=1

((divergences == 0)) || exit 1
echo "the dense and sparse backends agree byte for byte"
