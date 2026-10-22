# Each data source brings its own pass through the block update

The block update draws the field and both tree fields at every leaf pair. At each cell it multiplies in what every compiled-in data source says about that cell. How a data source took part was spread over about ten places, and none of them was the data source. The block model held one `#if` block per data source per hook. The three kernels took LOTUS, the simple error model and mass spectrometry as three named arguments. The data sources exposed members that existed only so that the block model could drive them from outside. Nothing tested what a real data source tallies during a real block update.

So each data source now brings its own pass. A data source hands a pass out when a block update begins (`begin_block_update()`), and takes it back when the update ends (`end_block_update(pass)`), to install the statistic the pass tallied. The block update reaches every pass through one set of hooks, combined by one fold.

## The six hooks

`DataSourcePass` (`src/data_sources/data_source_pass.h`) names them, in the order the block update calls them:

| Hook | When | Contract |
|---|---|---|
| begin pass | once, serial, before the parallel region | Ready anything that cannot be readied concurrently. |
| begin row | once per species leaf, in the parallel region | Given the species leaf and the row's first field cell. Returns the row's own state. |
| cell factor | once per cell | Given the row state and a molecule leaf of the row. Returns `{P(obs \| Y=0), P(obs \| Y=1)}`. May advance row-local state, such as a cursor. |
| record | once per cell, after its cell factor | Given the row state, the molecule leaf and the drawn field state only. |
| end row | once per row, in the parallel region | Folds the row's tally into the pass-wide one. |
| end pass | once, serial, after the parallel region | Finishes the pass-wide tally. Gets the field read-only, so a debug build can recount the statistic. |

Within a row, cells arrive in ascending molecule-leaf order, each once, with the cell factor before the record. This holds across tiles, because a row's state survives from one tile to the next. A debug build of the fold asserts it, so a traversal change that breaks it fails at the first wrong cell. Begin row, cell factor and record are `const`: many rows run at once, and the only state they may write is the row's own.

There is no per-column hook. No data source needs one for correctness. LOTUS could use one to settle its molecule class once per column, but that is a performance change, to be measured on its own.

## The pass is plain data

A pass holds what its data source observes the field with, and the source's current parameter values. The simple error model's pass (`simple_error_model::TBlockPass`, aliased `TSimpleErrorModel::TBlockPass`) holds D and the error rate. So a test builds a pass without a tree or stattools. The pass is the one seam through which a data source's statistic changes. Nothing outside the pass installs a disagreement count.

## The fold and its fixed order

`TBlockSources` (`src/data_sources/block_sources.h`) lists the compiled-in passes. It holds the only data-source `#ifdef`s the fold depends on; `TDataModel` names the data source that hands out each pass. `TBlockData` (`src/field/TBlockData.h`) holds one pass per entry, begins and ends them, fans the row hooks out, and multiplies the cell factors into the cell data factor. It is a compile-time fold over a tuple. There is no virtual call and no type-erased callable: the cell factor runs once per leaf pair, about 4.6e9 times per iteration on the full LOTUS data. `TBlockModel` is the two trees, unchanged, plus a `TBlockData`. Until LOTUS and mass spectrometry bring passes, it also keeps their per-cell factors and LOTUS's present-cell bookkeeping.

The kernels (`block_weights`, `block_probabilities`, `draw_block`) take one cell data factor instead of three named pairs.

The chain does not move. The kernel used to form `lotus * simple_error * mass_spec`, left to right, with `{1, 1}` for a source the build left out. The fold multiplies into a running product that starts at 1.0, in the fixed order LOTUS, the simple error model, mass spectrometry. `1.0 * x` is exact, so the product is the same double in every letter set. Floating-point multiplication is not associative, so the order is fixed. A data source not yet migrated to a pass keeps its old factor, multiplied in at its own position. LOTUS's factor is the running product the fold starts from, and mass spectrometry's multiplies the fold's result.

## Mass spectrometry takes part without satisfying `DataSource`

`DataSource` (`src/data_sources/data_source.h`) is a source's setup, likelihood, simulation and notifier stats. Mass spectrometry does not satisfy it, because it carries latent state of its own, the assignments, and they need a Metropolis update of their own every iteration. Its part in the block update is still per-cell shaped: given the current assignments, it answers `{P(runs | Y=0), P(runs | Y=1)}` for a cell. So it can satisfy `DataSourcePass` without satisfying `DataSource`. The two concepts are separate for that reason, and neither `DataSource` nor `ACOL_DATA_SOURCES` changes. Mass spectrometry tallies nothing: its likelihood reads the field itself.

## Thread-count invariance belongs to each data source

The old accumulator kept one slot per row for every data source, and merged them in row order after the parallel region. Its comment argued for that: every source counted integers, but a slot per row kept thread-count invariance from being a property the next source had to remember.

That reverses here. Each data source keeps its own tally and its own rule that the tally comes out the same at any thread count. The simple error model adds each row's disagreement count into its pass atomically, as LOTUS already adds its present cells (ADR-0015). Both are integer sums, so the order the rows end in does not matter. A test holds every data source to the rule: the `BlockData` suite runs the real fold through the real traversal at 1 thread and at 4, under every backend pairing, and checks the committed statistic against a fresh count of the field the update left.

The link counters keep their per-row slots and row-order merge. They are the block update's own tally, and the traversal returns them (`block_update::run`).

## Considered options

**Per-row slots for every data source, owned by the block update.** This is what the old accumulator did. It writes the reasoning about thread-count invariance once, but it is a shared layer whose fields had to be kept in step with the data sources: a per-cell outcome, a per-row outcome and an accumulator, each with one `#if` block per source. It also does not fit LOTUS, whose row tally is a row of the present-cell counts, so LOTUS was already the exception. Rejected.

**Virtual dispatch or type-erased callables over the passes.** This avoids a template fold, but the cell factor runs once per leaf pair, and an indirect call there keeps the compiler from inlining the fold back into the kernel. Rejected.

**A per-column hook.** No data source needs one for correctness. Left for a measured performance change.

**Keeping the single-cell field conditional.** ADR-0010 kept `field_math::prob_field_cell_is_one` with no production caller, as arithmetic a future repair might call. It took two named data pairs and left mass spectrometry out, so it was a second version of the formula that had already drifted. The block kernels hold the one formula for what the data does to a cell, and their tests cover what its tests did. It is deleted, with its four tests.

## Evidence

To be filled in by the last change of this migration: the byte-identity run against `main`, and the block update's CPU time before and after.
