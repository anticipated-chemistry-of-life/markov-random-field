# The dense-versus-sparse gate

Two binaries, built from the same sources and differing only in which storage the aliases in
`src/storages/storage_backend.h` select, have to produce byte-identical output from the same seed.
`run.sh` is the check; `just parity` runs it, and so does CI on every push
(`.github/workflows/backend-parity.yml`).

The comparison is between two *builds* because the storage is a compile-time choice. That is the
point: there is no dispatch layer that exists only for testing, so the gate exercises exactly what
ships.

## Which pairs it gates

The field and the node state choose their storage independently, so there are four pairings. The
gate builds two: sparse field against sparse node state, and dense against dense. Those two between
them exercise both storages. `src/storages/storage_backend.h` records which pairs are gated and
ADR-0006 gives the argument. The default build is dense against dense, so it is one of the two.

`run.sh` drives cmake itself rather than going through `just`, because nothing in the build system
chooses a storage any more. It passes `-DACOL_FIELD_STORAGE` and `-DACOL_NODE_STATE_STORAGE` on the
compiler command line, which is how an external define overrides an alias.

## The fixture

Four files, small enough that both chains run in well under a second:

| file            | shape                                                                   |
| --------------- | ----------------------------------------------------------------------- |
| `species.txt`   | a balanced binary tree: 15 nodes, of which 8 leaves and 7 internal       |
| `molecules.txt` | one root, three internal children, two leaves each: 10 nodes, 6 leaves   |
| `*_papers.txt`  | a distinct paper count per leaf, so research effort varies by cell       |

(Internal nodes include the root, per `CONTEXT.md`.)

The two trees are deliberately different shapes, so that no container the run asks for is square:
the field is 8x6, the species node state 7x6 and the molecules node state 8x4. A square
container is the one shape in which a row-walk and a column-walk can agree by accident, and telling
those two apart is what the sparse window does on every open.

The parameters are not fitted to anything and no result of the run is interpreted -- only that the
two backends produce the same bytes. What matters about the fixture is that it is varied enough for
a divergence to show, and small enough to run on every push.

## What it compares

Every file the runs write, except `*.log` -- which carries a fresh ntfy topic UUID and wall-clock
timings, and so differs between two runs of the *same* binary.

- `simulate` writes the field and both node states in full, the LOTUS and simple-error data
  drawn from them, and the per-iteration traces.
- `infer` writes the parameter traces, the field and node-state traces, the joint density and
  the posterior field.

Both backends run in the same working directory, one after the other, with the outputs moved aside
in between: `acol.parameters` echoes the command line with every path resolved against the working
directory, so two runs in two directories would differ there and nowhere else.

Everything compared lives under the two `--out` prefixes, and the script fails if a run leaves
anything beside them -- an output the comparison would not look at. That is why
`--write_branch_lengths` is not passed: `write_branch_length_grid` names its file from a literal
rather than from the prefix. The grid is a pure function of `n_bins` anyway, so it carries nothing
a backend could get wrong.

## Where it stops

Above 32767 iterations the two fields thin their posterior counters differently -- the sparse
counter shares its 16-bit word with the state bit, the dense one does not -- and the thinning factor
is also what decides which iterations get a trace line. The two would then write traces of different
lengths, by design rather than by regression. `run.sh` refuses a chain that long rather than let the
gate fail for a reason it is not testing.

One thread for the two backend chains. What that leaves out is the multi-batch commit of the
update's deferred inserts, which one thread never produces; `StorageEquivalence` in
`tests/TStorageConformance_Tests.cpp` covers that path instead.

## The thread-count check

A third run, at the foot of `run.sh`: the dense binary simulates the same chain again at four
threads (`ACOL_PARITY_THREADS`), and every file it writes has to match the one-thread run byte for
byte. `acol.parameters` is left out, because it echoes the command line and the command line is
where the two runs differ on purpose.

A cell's uniform is hashed from its position (ADR-0007), so a chain that draws nothing else gives
one answer however many threads it runs on. `simulate` is that chain: `IsSimulation` compiles the
alpha and nu moves out of the clique loop, and those moves hold the last draws still taken from the
thread-local generator.

`infer` is therefore **not** checked this way, and it would fail if it were. That is the half of
"reproducible at any thread count" that does not hold yet, and ADR-0007 records what a repair
costs.

## What it caught

The first thing this gate found, on its first CI run, was a real defect -- worth recording, because
it is the shape of the problem it exists to catch.

`TLotus::calculate_log_likelihood_of_L` splits its sum in two: the cells a cursor
yields, term by term through an accumulator, and every other cell folded into one closed-form
product. The cursor used to yield the cells the field *stored* -- which is a property of the
backend, not of the field: the sparse matrix holds the cells it was given, ones and zeros alike,
and the dense field holds all of them. So the two backends split the same sum differently and
reached answers about one ULP apart. That difference goes into a Metropolis ratio, which is compared
against a uniform draw; a few iterations later one chain accepted a move the other rejected, and
from there they were simply different chains.

The fix was to give the cursor a contract that does not mention the backend: it yields the cells
that are **one**. Every caller wanted exactly that anyway, a stored zero contributes precisely what
the closed-form term contributes for it, and the two backends now walk the same cells in the same
order. See `TStorageYMatrix::OnesCursor`.

The lesson generalises: anything that lets "what this backend happens to store" reach an
arithmetic result will diverge, and printed precision will hide it right up until it doesn't.
