#!/usr/bin/env bash
#
# The one entry point behind the configure/build/run/test/bin/dir tasks in
# pixi.toml. It parses [MODE] [FLAGS] off the front of the argument list,
# configures + builds the matching directory, then performs `action`:
#
#   bash scripts/acol.sh ACTION [MODE] [FLAGS] [ARGS...]
#
# ACTION is configure, build, run, test, bin or dir; MODE is debug (the default)
# or release; FLAGS is any combination of l, s and m (default `ls`, or
# ACOL_DEFAULT_FLAGS). Everything left over goes to the executable.
#
# Call it through pixi (`pixi run build release ls`), not directly: cmake, ninja
# and the compiler come from the pixi environment, and the check below is what
# says so when they do not.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Data sources used when a task is called without a FLAGS argument
default_flags="${ACOL_DEFAULT_FLAGS:-ls}"

if (($# == 0)); then
    echo "error: no action; expected one of configure, build, run, test, bin, dir" >&2
    exit 1
fi
action="$1"; shift

mode="debug"
case "${1-}" in
    debug|release) mode="$1"; shift ;;
esac

flags="$default_flags"
if [[ "${1-}" =~ ^[lsm]+$ ]]; then flags="$1"; shift; fi
# Allows `pixi run test release -- --gtest_filter=...`-style separation when an
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

# The tasks run from the workspace root, but a direct call need not, and every
# path below is relative to it.
cd "$ROOT"

# Must agree with binaryDir in CMakePresets.json.
export ACOL_FLAG_SUFFIX="-${key}"
build_dir="build/${mode}${ACOL_FLAG_SUFFIX}"

case "$action" in
    configure|build|run|test) ;;
    dir) echo "$build_dir"; exit 0 ;;
    bin) echo "$build_dir/acol"; exit 0 ;;
    *) echo "error: unknown action '$action'" >&2; exit 1 ;;
esac

# pixi exports CONDA_PREFIX when it activates the environment, so an unset one
# means the toolchain this build needs is not on PATH.
if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "error: no pixi environment; run this through pixi, e.g. 'pixi run $action'" >&2
    exit 1
fi

# Only configure when there is a reason to. Re-running cmake regenerates
# armadillo's headers, which invalidates every object that includes them --
# a ~20 file rebuild on every `pixi run run`. Ninja re-runs cmake by itself when
# CMakeLists.txt or the presets change, so skipping it here is safe.
cache="$build_dir/CMakeCache.txt"
if [[ "$action" == "configure" ]] || ! [[ -f "$cache" ]] \
   || ! grep -qx "LOTUS:BOOL=$lotus" "$cache" \
   || ! grep -qx "SIMPLE_DATA:BOOL=$simple" "$cache" \
   || ! grep -qx "USE_MS_DATA:BOOL=$ms" "$cache"; then
    # The conda compiler packages export CC/CXX (and the matching sysroot
    # flags) from their activation scripts; fall back to the plain names if
    # they did not.
    if [[ -z "${CXX:-}" ]]; then
        case "$(uname -s)" in
            Darwin) export CC="$CONDA_PREFIX/bin/clang" CXX="$CONDA_PREFIX/bin/clang++" ;;
            *)      export CC="$CONDA_PREFIX/bin/gcc"   CXX="$CONDA_PREFIX/bin/g++" ;;
        esac
    fi
    # The conda toolchain is pinned to the SDKs it knew about when it was built. A macOS SDK
    # newer than that -- the weeks after an OS or Xcode upgrade, before conda-forge catches up --
    # breaks it two ways. Its linker fails to link even an empty program. Its libc++ headers fail
    # to compile against the new `math.h`: `<random>` at C++20 stops at an undeclared `INFINITY`.
    # Probe both, at the standard this project compiles with, and fall back to the Xcode SDK,
    # which always matches its own compiler, rather than fail every build until conda-forge ships
    # a fix.
    if [[ "$(uname -s)" == "Darwin" ]] && command -v xcrun >/dev/null 2>&1; then
        probe="$(mktemp -d)"
        printf "int main(){return 0;}" > "$probe/t.c"
        printf "#include <random>\nint main(){return 0;}" > "$probe/t.cpp"
        if { ! "$CC" "$probe/t.c" -o "$probe/t_c" >/dev/null 2>&1 \
             || ! "$CXX" -std=c++20 "$probe/t.cpp" -o "$probe/t_cxx" >/dev/null 2>&1; } \
           && xcrun -f clang++ >/dev/null 2>&1; then
            export CC="$(xcrun -f clang)" CXX="$(xcrun -f clang++)"
            export SDKROOT="$(xcrun --show-sdk-path)"
        fi
        rm -rf "$probe"
    fi
    cmake --preset "$mode" -DLOTUS="$lotus" -DSIMPLE_DATA="$simple" -DUSE_MS_DATA="$ms"
fi

if [[ "$action" != "configure" ]]; then
    if [[ "$action" == "test" ]]; then target="acol_unitTests"; else target="acol"; fi
    cmake --build "$build_dir" --target "$target"
fi

# .clangd points at build/compile_commands.json; keep it on the last build.
ln -sfn "${mode}${ACOL_FLAG_SUFFIX}/compile_commands.json" build/compile_commands.json

case "$action" in
    run)  exec "./$build_dir/acol" "$@" ;;
    test) exec "./$build_dir/acol_unitTests" "$@" ;;
esac
