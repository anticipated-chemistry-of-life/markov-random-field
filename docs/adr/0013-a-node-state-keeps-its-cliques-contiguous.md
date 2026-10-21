# A node state keeps its own cliques contiguous

A tree's walks read its node state one clique at a time: the Gibbs walk, the density, and the likelihoods of alpha, nu and the branch lengths each step through every node of a clique. A clique is the run of cells along the tree's own dimension (ADR-0011), and the dense node state stored every run row-major. That suits the molecule tree, whose dimension is the last, and does not suit the species tree: consecutive nodes of a species clique sat a whole row apart, 141,618 bytes on the LOTUS data, so every node read loaded a cache line to use one byte of it. The profile of a LOTUS run put 38% of CPU time in `cell_is_one`, 98% of it under `TCliqueCells::is_one`.

The dense node state now takes the dimension to keep contiguous when it is sized, and `TTree` sizes its own with its own dimension contiguous. Each clique of either tree is one run of adjacent cells.

**Where a cell sits is not what it is called.** A cell is still named by its row-major linear index under every storage: the uniform it draws is hashed from that index (ADR-0007), a file lists it, and a deferred insert carries it. The storage converts a name to a *position* when it reads or writes. The storage interface grew the way in for a caller that walks positions itself — `position_of`, `is_one_at`, `locate_at` — and `TCliqueCells` uses it: it computes the run's first position and its step once, and steps positions from there, so a walk pays no conversion per node. The sparse storage keys its map by the linear index, so a position there is that index and the layout argument is checked and otherwise ignored.

Because a cell's name does not move, the chain does not move. The layout is an implementation detail of the dense storage, not of the model.

## The block update pays for it, and tiling pays part of it back

The block update draws leaf pairs `(species leaf, molecule leaf)` and reads both node states at each. In the leaf block the two node states are transposes of each other, so no loop order reads both contiguously: a field row walks the molecule node state along a row and now walks the species node state across a stride. Drawn a whole row at a time, the next row came back to cache lines long evicted.

The block update now gives a thread a *band* of 64 species rows and draws it a *tile* of 128 columns at a time, row after row. The band's 64 species leaves at one molecule leaf are 64 adjacent bytes of the species node state — one cache line — and every row of the band uses it before it is evicted. Each row still walks its own columns in ascending order, settles what it settles once, and hands over its tally and the model's bookkeeping once, when the band ends; the band only interleaves the rows. No leaf pair reads a cell another writes, so the interleaving draws what the rows drawn one after another would.

Tiling recovers the leaf cells and not the parents. A leaf pair also reads its species leaf's parent, and on the LOTUS tree the species leaves are in an order unrelated to the tree: along the 32,347 leaves in file order the parent changes 32,313 times, and a band of 64 leaves has 62 distinct parents spread over about 50 cache lines of a column. So the parent read is a cache line of its own per leaf pair, whatever the band and tile. In production this leaves the block update where it was — neither the gain nor the loss of the layout — rather than paying it back.

## Considered options

**Armadillo, or any column-major matrix.** Measured on the `bench/dense-vs-armadillo` branch: addressing an `arma::Mat` at the same subscripts makes the species walk fast and the block update and the molecule walk slow, and an iteration took twice as long as with the row-major storage (22.5 s against 11.4 s). Column-major is not the fix; the fix is that each tree's own dimension is contiguous, which no single fixed order gives both trees. Armadillo has no 1-bit type, so it saves no memory either.

**Transposing the species node state's shape instead of its layout** — sizing it molecule leaves by species nodes. That changes the subscripts and with them the linear indices, so every uniform, every file and the dense–sparse parity would move with it. Keeping the shape and changing only where a cell sits leaves all of those where they were.

**Leaving the block update row by row.** On the synthetic, LOTUS-shaped trees, whose leaves are sorted by parent, the untiled block update went from 3.1 s to 5.7 s per iteration once the species node state was transposed, eating half of what the walks gained.

## Evidence

`bench/dense_storage_bench.cpp` on the full LOTUS trees, 32 threads, unloaded machine, before tiling: species walk 4.8 s → 1.9 s, species density 2.1 s → 0.86 s, block update 3.1 s → 5.7 s, per iteration. With tiling and the production block update (`own` against `own_clique`, alternating in one process, on a machine at a load of ~140 on 64 cores): block update 8.8 / 7.7 s → 8.3 / 9.7 s, species walk 14.4 / 14.3 s → 4.7 / 5.8 s, species density 7.4 / 9.2 s → 2.3 / 2.9 s, iteration 36.8 / 38.6 s → 19.1 / 23.5 s. Every variant ended in the same state.

The production binary before (`de510c6`) and after, profiled with gperftools over a complete chain of 3 iterations on `inference_lotus/` (`--numBurnin 0 --thinning 1 --fixedSeed 42 --numThreads 32`), CPU seconds summed over threads: clique traversals 823 → 303, block update 298 → 306, the whole run, reading and writing included, 1978 → 1170. The chain took 5 min 13 s → 3 min 32 s. On a machine that loaded, wall time is a poor measure — two earlier pairs of 5-iteration runs came out 7 min 39 s / 7 min 58 s against 7 min 24 s / 7 min 10 s — and CPU time is the one that held still.

Same chain: with 5 iterations, every output file of the two binaries is byte-identical, the 43 GB `acol_Y_posterior.tsv` included; only the output path recorded in `acol.parameters` differs.

## Consequences

The dense node state's memory is no longer in linear-index order. Code that walks a node state's memory directly — `cells()` — walks positions, and must ask the storage which cell it is holding. Code that reads by linear index pays a division per read under the transposed layout; a hot path reads by subscript or by position instead (`tree_field_posterior` reads by subscript for this reason).

`BAND_ROWS` and `TILE_COLUMNS` are tuned to a cache line and a core's L2, not to the model. They change the order in which leaf pairs are drawn and nothing else.

The species leaf order is what keeps the block update from gaining. Ordering the species leaves so that siblings are adjacent and their parents follow them would make a band's parents a few cache lines instead of fifty. That moves the leaf order away from file order (ADR-0004) and with it every linear index and every uniform, so it is a different chain, and it is not part of this decision.
