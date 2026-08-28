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

Each combination gets its own build directory (`build/<mode>-<letters>-<backend>`, e.g.
`build/release-ls-sparse`), so switching back and forth does not trigger a rebuild.

### Choosing the storage backend

Which implementation backs the field and the internal state is likewise a compile-time decision,
made by `-DACOL_STORAGE_BACKEND=<name>`. It selects the `using` aliases in
`src/storages/storage_backend.h`; the interface those aliases have to satisfy is the pair of
concepts in `src/storages/storage_concepts.h`, checked with `static_assert` rather than through
virtual calls, so nothing on a storage access path pays for the choice.

| value    | field             | internal state    |
| -------- | ----------------- | ----------------- |
| `sparse` | `TStorageYMatrix` | `TStorageZMatrix` |
| `dense`  | `TStorageYDense`  | `TStorageZDense`  |

`sparse` is the default and is the path to larger runs; `dense` stores the whole container space
and is the obviously-correct implementation to check the other against. Anything else fails at
configure time. Unlike the data sources this is not a `just` letter but an environment variable,
because it is a different kind of choice: the letters change what the program computes, this does
not.

```bash
ACOL_BACKEND=dense just build release ls    # build/release-ls-dense/acol
cmake --preset debug -DACOL_STORAGE_BACKEND=dense
```

That the two agree is checked, not assumed. `just parity` builds both binaries, runs the same
simulation and the same chain under each from a fixed seed, and compares every file they write byte
for byte; it runs in CI on every push. See `tests/backend_parity/`.

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
the data-source letters and the storage backend, and it is empty when you invoke cmake yourself.
