# The LOTUS likelihood is read off counts of the field's present cells

Gamma and the LOTUS error rate are moved by Metropolis, and each move needs the likelihood of the whole of L at the proposed value. `TLotus::calculate_log_likelihood_of_L` computed it by merge-joining the field's ones with L's records: one pass over the field per move, on one thread. On the LOTUS data the field is 32,347 × 141,618 cells held densely, so a pass reads about 9 GB, and an iteration made three of them (two gammas and the error rate). A wall-clock profile put them at about a quarter of an iteration.

Neither parameter moves the field. And a cell's term depends on very little:

```
P(L = 1 | x = 1) = f_s(count of the species leaf) * f_m(count of the molecule leaf)      f_i(c) = 1 - exp(-gamma_i * log(c + 1))
P(L = 1 | x = 0) = error rate
```

So the likelihood is a function of the parameters and of how many present cells there are of each kind — (species leaf's paper count, molecule leaf's paper count, reported or not) — plus how many absent cells are reported and how many are not, which follow from the totals. The LOTUS trees have 465 and 459 distinct paper counts, so there are about 214,000 kinds where there are 4.6 billion cells.

`lotus_math::TPresentCells` holds those counts, indexed by each leaf's paper-count class (`TPaperCountClasses`), and `lotus_math::log_likelihood_of_L` reads the likelihood off them: two logarithms per kind. A move costs a few milliseconds and no pass over the field.

The counts change only when the field does, and the field changes only in the block update — which already visits every cell and knows whether LOTUS reports it. So the block update counts them as it draws: a row counts its present cells by molecule class (`TBlockModel::record`), and adds its row of counts to the field's when it ends (`TDataUpdateAccumulator::add_lotus_row`). The field's counts are counted from scratch once, when the data sources take their starting values (`lotus_math::count_present_cells`), and a debug build checks after every block update that what the rows counted is what a count from scratch finds.

This is the move ADR-0005 made for omega and the simple error model makes for its error rate: replace a pass over the cells by counts the block update keeps.

## Determinism

The counts are integers, so they come out the same whichever order the rows add them in. The rows add into one shared set of counts with atomic adds, rather than into a slot per row as `TDataUpdateAccumulator` does for everything else: a row's share is a row of the counts, and a slot per species leaf would be a copy of them per species leaf. The likelihood is then summed over the kinds in class order, so it is a function of the counts and the parameters alone — the same under either backend and at any thread count, which the parity gate checks.

The block update no longer sums the LOTUS likelihood per cell, so the per-row `TSumLogProbability` that existed to make that sum schedule-independent is gone.

## What changes for a chain

The likelihood is the same quantity summed in another order, so it differs from the old one in its last bits, and a Metropolis decision that lands exactly on that difference can go the other way. Over 5 iterations on `inference_lotus/` (15 decisions) the two builds accepted and rejected the same moves and wrote identical traces, link counters, node-state traces and posterior summaries.

## Considered options

**Parallelising the pass.** Splitting the merge join across threads, with a partial sum per fixed slice summed in slice order, keeps the sum deterministic and divides its time by up to the thread count. It still reads 9 GB per move, so memory bandwidth bounds it; this was not measured.

**Counting lazily, once per iteration, from the field.** The first move after a block update counts the field and the rest reuse the counts. One pass per iteration instead of three, and no change to the block update. It was rejected because the block update already holds every fact the counts need, at the moment it decides them.

**Separating the reported present cells.** `log(f_s f_m) = log f_s + log f_m`, so the reported present cells need only per-leaf counts. The unreported ones do not separate — `log(1 - f_s f_m)` — so the class-pair table is needed anyway, and one table for both keeps one code path.

## Evidence

The chain on `inference_lotus/` (`--iterations 5 --numBurnin 0 --thinning 1 --fixedSeed 42 --numThreads 32`), alternating the two builds on a machine at a load of 15–22, stopped once the chain ended: 2 min 29 s / 2 min 29 s before, 2 min 4 s / 2 min 3 s after. About 57 s of each is the chain start, which runs once, so an iteration went from about 18.4 s to about 13.3 s.
