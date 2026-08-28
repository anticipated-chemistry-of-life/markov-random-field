# acol -- build/run helpers.
#
#   just setup                      create or update the micromamba environment
#   just build [MODE] [FLAGS]       configure + compile ./acol
#   just run   [MODE] [FLAGS] ARGS  build, then run ./acol ARGS
#   just test  [MODE] [FLAGS] ARGS  build, then run ./acol_unitTests ARGS
#
# MODE  is `debug` (the default) or `release`.
# FLAGS is any combination of the letters below; it defaults to `ls`.
#
#     l -> -DLOTUS=ON        (LOTUS data)
#     s -> -DSIMPLE_DATA=ON  (simple error model data)
#     m -> -DUSE_MS_DATA=ON  (mass spec data)
#
# Anything left after MODE and FLAGS is forwarded verbatim to the executable:
#
#     just run release lsm --input file.txt --scalar 0.6
#
# is `./acol --input file.txt --scalar 0.6` compiled at -O3 with all three data
# sources.
#
# Which storage backs the field and the internal state is a third axis, set by
# ACOL_BACKEND (`sparse`, the default, or `dense`) rather than by an argument,
# because it changes no behaviour -- that is what `just parity` checks.
#
# Every combination gets its own build directory
# (build/<MODE>-<FLAGS>-<BACKEND>), so switching back and forth does not force a
# rebuild.
#
# All cmake/compiler/executable invocations run inside the micromamba
# environment, so there is nothing to activate by hand.

set positional-arguments

# Name of the micromamba environment; override with ACOL_ENV=... just build
conda_env := env('ACOL_ENV', 'acol_env')

# Data sources used when the recipe is called without a FLAGS argument
default_flags := env('ACOL_DEFAULT_FLAGS', 'ls')

# Storage backend for the field and the internal state: `sparse` or `dense`
backend := env('ACOL_BACKEND', 'sparse')

[private]
default:
    @just --list --unsorted

# Create (or update) the micromamba environment used to build acol.
setup:
    #!/usr/bin/env bash
    set -euo pipefail
    mamba="${MAMBA_EXE:-micromamba}"
    command -v "$mamba" >/dev/null 2>&1 || {
        echo "error: micromamba not found (set MAMBA_EXE to its path)" >&2
        exit 1
    }

    # git is deliberately not in here: FetchContent uses the system git, which
    # you already need to have cloned this repository.
    packages=(
        "cmake>=3.30" ninja
        cxx-compiler
        zlib fmt armadillo openssl nlohmann_json
        "libopenblas=*=*openmp*"
    )
    # On Linux the conda gcc brings its own libgomp; forcing llvm-openmp there
    # would flip _openmp_mutex to the llvm variant and break it.
    if [[ "$(uname -s)" == "Darwin" ]]; then
        packages+=(llvm-openmp)
    fi

    if "$mamba" run -n "{{ conda_env }}" true >/dev/null 2>&1; then
        echo "==> updating environment '{{ conda_env }}'"
        "$mamba" install -y -n "{{ conda_env }}" -c conda-forge "${packages[@]}"
    else
        echo "==> creating environment '{{ conda_env }}'"
        "$mamba" create -y -n "{{ conda_env }}" -c conda-forge "${packages[@]}"
    fi

# Configure the build directory without compiling anything.
configure *args:
    @just _drive configure "$@"

# Compile ./acol.
build *args:
    @just _drive build "$@"

# Compile ./acol and run it with the remaining arguments.
run *args:
    @just _drive run "$@"

# Compile and run the unit tests. Extra arguments go to gtest (e.g. --gtest_filter=...).
test *args:
    @just _drive test "$@"

# Print the path of the acol binary for a given MODE/FLAGS (used by scripts).
bin *args:
    @just _drive bin "$@"

# Print the build directory for a given MODE/FLAGS.
dir *args:
    @just _drive dir "$@"

# Build both storage backends and check they produce byte-identical output.
parity *args:
    @bash tests/backend_parity/run.sh "$@"

# Open a shell inside the micromamba environment.
shell:
    @"${MAMBA_EXE:-micromamba}" run -n "{{ conda_env }}" "${SHELL:-bash}"

# Delete every build directory.
clean:
    rm -rf build

# Delete build directories *and* the fetched coretools/stattools checkouts.
distclean: clean
    rm -rf coretools stattools

# Shared entry point: parses [MODE] [FLAGS] off the front of the argument list,
# configures + builds the matching directory, then performs `action`.
[private]
_drive action *args:
    #!/usr/bin/env bash
    set -euo pipefail

    action="$1"; shift

    mode="debug"
    case "${1-}" in
        debug|release) mode="$1"; shift ;;
    esac

    flags="{{ default_flags }}"
    if [[ "${1-}" =~ ^[lsm]+$ ]]; then flags="$1"; shift; fi
    # Allows `just run release -- --gtest_filter=...`-style separation when an
    # argument would otherwise be mistaken for MODE or FLAGS.
    if [[ "${1-}" == "--" ]]; then shift; fi

    lotus=OFF; simple=OFF; ms=OFF; key=""
    case "$flags" in *l*) lotus=ON; key="${key}l" ;; esac
    case "$flags" in *s*) simple=ON; key="${key}s" ;; esac
    case "$flags" in *m*) ms=ON; key="${key}m" ;; esac

    if [[ "$lotus" == OFF && "$simple" == OFF ]]; then
        echo "error: flags '$flags' compile in no source for Y; add 'l' and/or 's'" >&2
        exit 1
    fi

    backend="{{ backend }}"
    case "$backend" in
        dense|sparse) ;;
        *) echo "error: ACOL_BACKEND='$backend' is neither 'dense' nor 'sparse'" >&2; exit 1 ;;
    esac

    # Must agree with binaryDir in CMakePresets.json.
    export ACOL_FLAG_SUFFIX="-${key}-${backend}"
    build_dir="build/${mode}${ACOL_FLAG_SUFFIX}"

    case "$action" in
        dir) echo "$build_dir"; exit 0 ;;
        bin) echo "$build_dir/acol"; exit 0 ;;
    esac

    mamba="${MAMBA_EXE:-micromamba}"
    command -v "$mamba" >/dev/null 2>&1 || {
        echo "error: micromamba not found (set MAMBA_EXE to its path)" >&2
        exit 1
    }
    "$mamba" run -n "{{ conda_env }}" true >/dev/null 2>&1 || {
        echo "error: micromamba environment '{{ conda_env }}' is missing or broken; run 'just setup'" >&2
        exit 1
    }
    in_env() { "$mamba" run -n "{{ conda_env }}" "$@"; }

    # Only configure when there is a reason to. Re-running cmake regenerates
    # armadillo's headers, which invalidates every object that includes them --
    # a ~20 file rebuild on every `just run`. Ninja re-runs cmake by itself when
    # CMakeLists.txt or the presets change, so skipping it here is safe.
    cache="$build_dir/CMakeCache.txt"
    if [[ "$action" == "configure" ]] || ! [[ -f "$cache" ]] \
       || ! grep -qx "LOTUS:BOOL=$lotus" "$cache" \
       || ! grep -qx "SIMPLE_DATA:BOOL=$simple" "$cache" \
       || ! grep -qx "USE_MS_DATA:BOOL=$ms" "$cache" \
       || ! grep -qx "ACOL_STORAGE_BACKEND:STRING=$backend" "$cache"; then
        # The conda compiler packages export CC/CXX (and the matching sysroot
        # flags) from their activation scripts; fall back to the plain names if
        # they did not.
        in_env bash -eu -c '
            if [[ -z "${CXX:-}" ]]; then
                case "$(uname -s)" in
                    Darwin) export CC="$CONDA_PREFIX/bin/clang" CXX="$CONDA_PREFIX/bin/clang++" ;;
                    *)      export CC="$CONDA_PREFIX/bin/gcc"   CXX="$CONDA_PREFIX/bin/g++" ;;
                esac
            fi
            exec cmake --preset "$1" -DLOTUS="$2" -DSIMPLE_DATA="$3" -DUSE_MS_DATA="$4" \
                 -DACOL_STORAGE_BACKEND="$5"
        ' _ "$mode" "$lotus" "$simple" "$ms" "$backend"
    fi

    if [[ "$action" != "configure" ]]; then
        if [[ "$action" == "test" ]]; then target="acol_unitTests"; else target="acol"; fi
        in_env cmake --build "$build_dir" --target "$target"
    fi

    # .clangd points at build/compile_commands.json; keep it on the last build.
    ln -sfn "${mode}${ACOL_FLAG_SUFFIX}/compile_commands.json" build/compile_commands.json

    case "$action" in
        run)  exec "$mamba" run -n "{{ conda_env }}" "./$build_dir/acol" "$@" ;;
        test) exec "$mamba" run -n "{{ conda_env }}" "./$build_dir/acol_unitTests" "$@" ;;
    esac
