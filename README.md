# Metabolite inference

## Compiling the repo

The build is driven by [`just`](https://github.com/casey/just), which wraps cmake and takes care of
the micromamba environment for you. Install the two prerequisites once:

```bash
brew install just          # or: cargo install just
"${SHELL}" <(curl -L micro.mamba.pm)
```

Then clone the repo and create the build environment:

```bash
git clone https://github.com/anticipated-chemistry-of-life/markov-random-field
cd markov-random-field
just setup
```

`just setup` creates (or updates) a micromamba environment called `acol_env` with cmake, ninja, a
compiler toolchain and the third-party libraries. Every other recipe runs inside that environment,
so there is nothing to `micromamba activate` by hand.

`coretools` and `stattools` are no longer git submodules: cmake checks them out via `FetchContent`
on the first configure.

### Building and running

```bash
just build                 # debug build of ./acol
just build release         # release build
just run                   # build and run ./acol
just test                  # build and run the unit tests
just test release          # ... in release mode
```

Everything after the build mode is forwarded verbatim to the executable, so

```bash
just run release --out results/acol --numThreads all
```

is the same as running `./acol --out results/acol --numThreads all` from a release build.

### Choosing the data sources

Which sources of information get compiled in is a compile-time decision. Pass any combination of
these letters right after the build mode:

| letter | cmake option       | data source             |
| ------ | ------------------ | ----------------------- |
| `l`    | `-DLOTUS=ON`       | LOTUS data              |
| `s`    | `-DSIMPLE_DATA=ON` | simple error model data |
| `m`    | `-DUSE_MS_DATA=ON` | mass spec data          |

The default is `ls`. At least one of `l` and `s` is required — with neither, nothing informs `Y` and
`src/Types.h` fails a `static_assert`.

```bash
just build l               # debug, LOTUS only
just run release lsm --out results/acol --numThreads all
```

Each combination gets its own build directory (`build/<mode>-<letters>`, e.g. `build/release-ls`),
so switching back and forth does not trigger a rebuild.

### Choosing the storage

Which storage backs the field, and which backs the node state, are two aliases in
`src/storages/storage_backend.h`. The two are chosen independently, and changing one is an edit to
one line. The build system takes no part.

| storage | field             | node state        |
| ------- | ----------------- | ----------------- |
| sparse  | `TStorageYMatrix` | `TStorageZMatrix` |
| dense   | `TStorageYDense`  | `TStorageZDense`  |

The interface an alias has to satisfy is the pair of concepts in
`src/storages/storage_concepts.h`, checked with `static_assert` rather than through virtual calls,
so nothing on a storage access path pays for the choice.

The default pairs a **sparse field with a dense node state**. Fill decides that pairing, not size.
`docs/adr/0006-each-storage-brings-its-own-traversal.md` gives the argument, and the header records
which pairings CI gates.

An external define overrides either alias:

```bash
cmake --preset debug -DCMAKE_CXX_FLAGS="${CXXFLAGS:-} -DACOL_FIELD_STORAGE=TStorageYDense"
```

That is how `just parity` builds two binaries from one source tree. It gates two of the four
pairings, sparse against sparse and dense against dense. It runs the same simulation and the same
chain under each from a fixed seed, then compares every file they write byte for byte. It runs in
CI on every push. See `tests/backend_parity/`.

Other recipes: `just configure` (configure only), `just bin` / `just dir` (print the binary or build
directory path), `just shell` (a shell inside the environment), `just clean`, `just distclean`.
Run `just` with no arguments for the full list.

### Using cmake directly

`just` is a convenience wrapper; the presets in `CMakePresets.json` work on their own. Activate the
environment first, since the presets do not do it for you:

```bash
micromamba activate acol_env
cmake --preset debug -DLOTUS=ON -DSIMPLE_DATA=ON
cmake --build build/debug
```

The presets put their output in `build/<preset>$ACOL_FLAG_SUFFIX`; `just` sets `ACOL_FLAG_SUFFIX` to
the data-source letters, and it is empty when you invoke cmake yourself.
