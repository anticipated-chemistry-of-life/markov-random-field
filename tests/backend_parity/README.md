# The dense-versus-sparse gate

Two binaries, built from the same sources and differing only in which storage the alias in
`src/storages/storage_backend.h` selects, have to produce byte-identical output from the same seed.
`run.sh` is the check; `just parity` runs it, and so does CI on every push
(`.github/workflows/backend-parity.yml`).

The comparison is between two *builds* because the backend is a compile-time choice. That is the
point: there is no dispatch layer that exists only for testing, so the gate exercises exactly what
ships.

## The fixture

Four files, small enough that both chains run in well under a second:

| file            | shape                                                                   |
| --------------- | ----------------------------------------------------------------------- |
| `species.txt`   | a balanced binary tree: 15 nodes, of which 8 leaves and 7 internal       |
| `molecules.txt` | one root, three internal children, two leaves each: 10 nodes, 6 leaves   |
| `*_papers.txt`  | a distinct paper count per leaf, so research effort varies by cell       |

(Internal nodes include the root, per `CONTEXT.md`.)

The two trees are deliberately different shapes, so that no container the run asks for is square:
the field is 8x6, the species internal state 7x6 and the molecules internal state 8x4. A square
container is the one shape in which a row-walk and a column-walk can agree by accident, and telling
those two apart is what the sparse implementation's `fill_current_state` does for a living.

The parameters are not fitted to anything and no result of the run is interpreted -- only that the
two backends produce the same bytes. What matters about the fixture is that it is varied enough for
a divergence to show, and small enough to run on every push.

## What it compares

Every file the runs write, except `*.log` -- which carries a fresh ntfy topic UUID and wall-clock
timings, and so differs between two runs of the *same* binary.

- `simulate` writes the field and both internal states in full, the LOTUS and simple-error data
  drawn from them, and the per-iteration traces.
- `infer` writes the parameter traces, the field and internal-state traces, the joint density and
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

One thread, always. The field sweep is parallel and draws from a shared random generator, so a run
is not reproducible across thread counts even under one backend (issue #38). What that leaves out
is the multi-batch commit of the sweep's deferred inserts, which one thread never produces;
`StorageEquivalence` in `tests/TStorageConformance_Tests.cpp` covers that path instead.

Byte-identical is not bit-identical arithmetic. Where a likelihood is summed by merge-joining the
stored cells (`TLotus::_calculate_log_likelihood_of_L_no_collapsing`, `count_disagreements`), the
sparse field visits the cells it holds and folds the rest into one closed-form term, while the
dense field visits every cell. The two compute the same quantity by different summation orders, so
they can differ in the last bits of a double and still print the same digits. The gate compares
what is written, which is the thing that has to agree; a divergence that only appears in the last
bits is not something it can see, and one that grows into the printed digits is.
