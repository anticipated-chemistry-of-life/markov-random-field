# Metabolite inference

## Compiling the repo

The build is driven by [`pixi`](https://pixi.sh), which wraps cmake and takes care of the conda
environment for you. Install the one prerequisite once:

```bash
brew install pixi          # or: curl -fsSL https://pixi.sh/install.sh | bash
```

Then clone the repo and create the build environment:

```bash
git clone https://github.com/anticipated-chemistry-of-life/markov-random-field
cd markov-random-field
pixi install
```

`pixi install` creates (or updates) the environment in `.pixi/envs/default` with cmake, ninja, a
compiler toolchain and the third-party libraries, as `pixi.toml` lists them. Every task runs inside
that environment, so there is nothing to activate by hand; `pixi shell` opens a shell in it.

`coretools` and `stattools` are no longer git submodules: cmake checks them out via `FetchContent`
on the first configure.

### Building and running

```bash
pixi run build             # debug build of ./acol
pixi run build release     # release build
pixi run run               # build and run ./acol
pixi run test              # build and run the unit tests
pixi run test release      # ... in release mode
```

The build mode and the data-source letters may be given in either order, and either may be left
out. They are typed task arguments, so everything meant for the executable goes after `--`:

```bash
pixi run run release -- --out results/acol --numThreads all
```

is the same as running `./acol --out results/acol --numThreads all` from a release build.

### Choosing the data sources

Which sources of information get compiled in is a compile-time decision. Pass one of the letter
sets, before or after the build mode:

| letter | cmake option       | data source             |
| ------ | ------------------ | ----------------------- |
| `l`    | `-DLOTUS=ON`       | LOTUS data              |
| `s`    | `-DSIMPLE_DATA=ON` | simple error model data |
| `m`    | `-DUSE_MS_DATA=ON` | mass spec data          |

The default is `ls`. The task takes the seven sets `l`, `m`, `s`, `ls`, `lm`, `sm` and `lsm`, and
names the rest — a set is spelled in l, s, m order, so `ls` is a set and `sl` is a typo. At least
one of `l` and `s` is required: with neither, nothing informs `Y`, and cmake stops before it
compiles anything.

```bash
pixi run build l           # debug, LOTUS only
pixi run build lsm         # debug, all three
pixi run run release lsm -- --out results/acol --numThreads all
```

Each combination gets its own build directory (`build/<mode>-<letters>`, e.g. `build/release-ls`),
so switching back and forth does not trigger a rebuild.

### Choosing the storage

Which storage backs the field, and which backs the node state, are two aliases in
`src/storages/storage_backend.h`. The two are chosen independently, and changing one is an edit to
one line. The build system takes no part.

| storage | field             | node state        |
| ------- | ----------------- | ----------------- |
| sparse  | `TStorageYSparse` | `TStorageZSparse` |
| dense   | `TStorageYDense`  | `TStorageZDense`  |

A sparse storage holds its cells in a hash map keyed by the linear index, so its memory tracks the
number of ones rather than the size of the container space. A dense one holds an array over the
whole space. Both hold the same cell, and a cell the storage has no entry for reads as state 0.

A stored cell costs a map node and a bucket slot, which is some tens of bytes where the cell itself
is one or two. Sparse therefore wins on memory well below **one one in twenty cells**, and not
merely below one in two. Choose it on the fill of the container and not on its size.

The interface an alias has to satisfy is the pair of concepts in
`src/storages/storage_concepts.h`, checked with `static_assert` rather than through virtual calls,
so nothing on a storage access path pays for the choice.

#### The posterior field is thinned

A field cell packs its posterior counter into the same 16-bit word as its state, which leaves the
counter 15 bits. A chain of `n` iterations is therefore counted one iteration in
`ceil(n / 32767)`, and that factor also decides which iterations get a trace line. Each run reports
it to its log file.

Both backends hold that cell, so both thin a chain identically. **A dense run before that was true
counted one iteration in `ceil(n / 65535)`**, because the dense counter had a 16th bit of its own.
A dense chain longer than 32767 iterations therefore writes a posterior field at half the previous
resolution, and traces of half the previous length. The numbers a run produces are otherwise
unchanged.

Both defaults are **dense** for now, which is one of the two pairings CI gates.
`docs/adr/0006-each-storage-brings-its-own-traversal.md` argues for a sparse field against a dense
node state, on fill rather than size; that pairing is one line away when the runs need it. The
header records which pairings CI gates.

An external define overrides either alias:

```bash
cmake --preset debug -DCMAKE_CXX_FLAGS="${CXXFLAGS:-} -DACOL_FIELD_STORAGE=TStorageYDense"
```

That is how `pixi run parity` builds two binaries from one source tree. It gates two of the four
pairings, sparse against sparse and dense against dense. It runs the same simulation and the same
chain under each from a fixed seed, then compares every file they write byte for byte. It runs in
CI on every push. See `tests/backend_parity/`.

Other tasks: `pixi run configure` (configure only), `pixi run bin` / `pixi run dir` (print the
binary or the build directory path), `pixi run clean`, `pixi run distclean`. `pixi task list` lists
them all, and `pixi shell` opens a shell inside the environment.

Every task is a cmake invocation in `pixi.toml`; there is no wrapper script. Which compiler a
build uses is decided in `cmake/toolchain.cmake`, which the presets name, so a plain
`cmake --preset` picks the same one.

### Using cmake directly

The pixi tasks are a convenience wrapper; the presets in `CMakePresets.json` work on their own.
Enter the environment first, since the presets do not do it for you:

```bash
pixi shell
cmake --preset debug -DLOTUS=ON -DSIMPLE_DATA=ON
cmake --build build/debug
```

The presets put their output in `build/<preset>$ACOL_FLAG_SUFFIX`; the tasks set
`ACOL_FLAG_SUFFIX` to the data-source letters, and it is empty when you invoke cmake yourself.
`cmake/toolchain.cmake` still chooses the compiler, because `CMakePresets.json` names it.
