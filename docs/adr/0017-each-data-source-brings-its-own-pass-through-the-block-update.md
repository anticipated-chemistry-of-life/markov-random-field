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

A pass holds what its data source observes the field with, and the source's current parameter values. The simple error model's pass (`simple_error_model::TBlockPass`, aliased `TSimpleErrorModel::TBlockPass`) holds D and the error rate. LOTUS's pass (`lotus_math::TBlockPass`, aliased `TLotus::TBlockPass`) holds the records, the paper-count classes and the reporting model, and the present cells it counts. Mass spectrometry's pass (`mass_spec_math::TBlockPass`, aliased `TMSMSData::TBlockPass`) holds the runs of every species leaf, and the filter and contamination probabilities read once when it is handed out. So a test builds a pass without a tree or stattools. The pass is the one seam through which a data source's statistic changes. Nothing outside a pass installs a disagreement count or present cells.

LOTUS's pass asks its record cursor about a cell once, in the cell factor, and keeps the answer in the row state for the record. Which kind a present cell is counted as -- its molecule class, and whether it is reported -- is decided in the pass alone.

## The fold and its fixed order

`TBlockSources` (`src/data_sources/block_sources.h`) lists the compiled-in passes. It holds the only data-source `#ifdef`s the fold depends on; `TDataModel` names the data source that hands out each pass. `TBlockData` (`src/field/TBlockData.h`) holds one pass per entry, begins and ends them, fans the row hooks out, and multiplies the cell factors into the cell data factor. It is a compile-time fold over a tuple. There is no virtual call and no type-erased callable: the cell factor runs once per leaf pair, about 4.6e9 times per iteration on the full LOTUS data. `TBlockModel` is the two trees, unchanged, plus a `TBlockData`. It names no data source.

The kernels (`block_weights`, `block_probabilities`, `draw_block`) take one cell data factor instead of three named pairs.

The chain does not move. The kernel used to form `lotus * simple_error * mass_spec`, left to right, with `{1, 1}` for a source the build left out. The fold multiplies into a running product that starts at 1.0, in the fixed order LOTUS, the simple error model, mass spectrometry. `1.0 * x` is exact, so the product is the same double in every letter set. Floating-point multiplication is not associative, so the order is fixed. While the sources moved onto passes one at a time, a source not yet migrated kept its old factor, multiplied in at its own position. The filter and contamination probabilities do not move during a block update, so mass spectrometry's pass reading them once gives the doubles a read per cell gave.

## Mass spectrometry takes part without satisfying `DataSource`

`DataSource` (`src/data_sources/data_source.h`) is a source's setup, likelihood, simulation and notifier stats. Mass spectrometry does not satisfy it, because it carries latent state of its own, the assignments, and they need a Metropolis update of their own every iteration. Its part in the block update is still per-cell shaped: given the current assignments, it answers `{P(runs | Y=0), P(runs | Y=1)}` for a cell, and `{1, 1}` for a species without runs. So it satisfies `DataSourcePass` without satisfying `DataSource`. The two concepts are separate for that reason, and neither `DataSource` nor `ACOL_DATA_SOURCES` changes. Mass spectrometry tallies nothing: its likelihood reads the field itself, so its record, end-row and end-pass hooks are empty. Its `end_block_update` only takes back the buffer the pass read the filter probabilities into, so the next update refills it in place.

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

Both checks compare `main` before the first change of this migration (f95dad5) with `main` after the last one (316c97a). They run on the balanced 255/511-node dataset under model validation (`model_validation/s_balanced_255_m_balanced_511`, 128 × 256 leaf pairs). It has simulated LOTUS records and simple error model data. Every run is 200 iterations at `--fixedSeed 42`. The dataset has no mass spectrometry file. So both sides read the same throwaway one, which is not committed. It has 4 species leaves, 3 runs each, 4 features per run, and 5 molecule-leaf candidates per feature with calibration scores.

**Byte identity.** Both binaries are `release-lsm` builds. Every run writes to the same output prefix, and the script moves its files aside. So even the files that echo the output path compare. The comparison leaves out `acol.log`, because it carries wall-clock timings.

| threads | assignment update of mass spectrometry | before against before | before against after |
|---|---|---|---|
| 1 | as shipped | | 20 of 20 files identical |
| 4 | as shipped | 2 of 20 files identical | 2 of 20 files identical |
| 4 | on one thread, both sides | 20 of 20 files identical | 20 of 20 files identical |

The runs write every trace they can: the parameter traces of the model and of both trees, the field trace and both node-state traces. The 20 files are those traces and the link-counter trace. They are also the joint density, the field posterior, both tree field posteriors, the posterior modes, the state posteriors, the mean and variance file, the config and the parameters.

At 4 threads, `main` before the migration differs from itself. `TMSMSData::update_all_MS_assignments` draws from a per-thread generator. So the 4-thread comparison runs that assignment update on one thread on both sides. A throwaway one-line change to its `num_threads` does this. The block update still runs on 4 threads. Those 4-thread runs also equal the 1-thread runs, except `acol.parameters`, which echoes `--numThreads`.

**CPU time.** Both binaries are `profile-lsm` builds. Each ran the same chain at 4 threads, without the field and node-state traces. The runs alternate before and after, five times each, at `CPUPROFILE_FREQUENCY=1000`. Each profile is symbolized with its own binary into a proto before any comparison. The block update's cumulative time is the profile focused on `^block_update::run$|^TMarkovField::_update_block$`:

| run | before | after |
|---|---|---|
| 1 | 0.56 s | 0.47 s |
| 2 | 0.67 s | 0.28 s |
| 3 | 0.70 s | 0.47 s |
| 4 | 0.68 s | 0.40 s |
| 5 | 0.73 s | 0.47 s |
| mean | 0.67 s | 0.42 s |

The block update is no slower. `pprof -diff_base` on run 3 puts most of the drop on two costs the old block model paid in release builds. The first is LOTUS's present-cell recount in `TMarkovField::_update_block`. It was a `DEBUG_ASSERT`, and coretools evaluates that argument in every build (64 ms). The second is the lookup of LOTUS through `TDataModel::get_lotus` at every cell (67 ms). The recount now sits in LOTUS's end-pass hook under `#ifndef NDEBUG`. The fold leaves no frame of its own. `TBlockData` and the cell factor and record of every pass are inlined into `block_update::update_band`.

These numbers have three limits. At 1 kHz the profiler loses signals while many threads are busy, so the absolute times undercount. Both sides ran under the same setting. At the default 100 Hz, three runs each gave 0.89, 0.42 and 0.41 s before, and 0.40, 0.49 and 0.46 s after. That is about 45 samples a run, too few to tell the two apart. The field of this dataset stays in cache. So the profile shows that the fold inlines back into the kernel. It does not show how the change behaves under the memory pressure of the full LOTUS data.
