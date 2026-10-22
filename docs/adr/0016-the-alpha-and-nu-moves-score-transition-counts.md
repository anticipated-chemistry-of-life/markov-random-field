# The alpha and nu moves score a clique's transition counts

The alpha and nu moves score every clique under three transition grids: the installed grid, alpha's candidate and nu's candidate. `tree_update::detail::clique_LL` computed each score node by node with a `TSumLogProbability`. The alpha pass scored two grids in one walk over the clique's nodes. The nu move then made a second walk over every clique to score its own grid. On `inference_lotus/` that is 6.1e9 species-tree cells and 4.6e9 molecule-tree cells per walk. A CPU profile put `clique_LL` at 11% of a run. Most of that was the exponent tracking in `TDynamicBinnedSumLog::add`, whose multiplications form one long dependency chain.

A node's term depends on three things: the bin its branch sits in, its parent's state and its own state. A root's term depends on its state alone. So a clique's log-likelihood under any grid is

```
sum over (bin, from, to) of  n(bin, from, to) * log P_bin(from -> to)  +  sum over s of  r(s) * log pi(s)
```

That is `4 * n_bins + 2` integers per clique, 42 at the default 10 bins.

Both moves score one configuration. The walk runs in the alpha pass and not again, and the bins move only in the branch-length move, after nu. So the alpha pass counts each clique once after the walk (`transition_counts::count_transitions`). It reads both of alpha's scores off the counts (`transition_counts::log_likelihood`). The nu move reads its score off the same counts and opens no clique view. `tree_update::run` holds the counts (`TTransitionCountTable`, one block per clique) for one tree update only. The next block update moves the leaves.

A kind of term that no node falls in adds nothing, and its probability is not read. Zero times the log of a zero probability is NaN.

## Determinism

The counts are integers. They come out the same in any node order. A score adds the kinds in one fixed order: bin, from, to, then the roots. So a score depends on the counts and the grid alone. It is the same under either backend and at any thread count. The parity gate checks this.

## What changes for a chain

A score is the same quantity summed in another order. It differs from the old score in its last bits. A Metropolis decision that lands exactly on that difference can go the other way.

Over 5 iterations on `inference_lotus/`, the two builds made the same decisions. They wrote identical parameter traces, link counters and posterior summaries. Over 10 iterations they wrote identical parameter traces.

## Considered options

**Count inside the walk.** The walk visits every edge while it scores a node against its children. But a node's state is drawn after its children are read. Counting there needs a second loop over the children, or a per-node histogram that is cleared at every internal node. It also leaves `FixZ` without counts, because `FixZ` skips the walk. A separate count after the walk reads states that are still in cache.

**Log tables in the transition-grid table.** These would replace the products in the walk and in the branch-length ratio as well. That is a separate change, and it does not depend on this one.

## Evidence

The chain on `inference_lotus/` (`--numBurnin 0 --thinning 1 --fixedSeed 42 --numThreads 32`, LOTUS and mass spectrometry compiled in) alternated the two release builds on a machine at a load of about 8. Times are from the chain's progress lines, after the first iteration, which holds the chain start.

| run | iterations | before | after |
|---|---|---|---|
| 1 | 5 | 11.75 s per iteration | 10.25 s per iteration |
| 2 | 10 | 11.56 s per iteration | 10.00 s per iteration |
| 3 | 10 | 11.44 s per iteration | 10.00 s per iteration |

An iteration went from about 11.5 s to about 10.0 s.

The profile builds (`profile-lm`, 5 iterations, whole run) compared with `pprof -diff_base`:

| | before | after |
|---|---|---|
| alpha pass (walk, counts, two scores) | 187.2 s | 148.6 s |
| `clique_LL` / `count_transitions` | 88.8 s | 17.2 s |
| nu move | 35.2 s | 0.05 s |
| whole run | 819.3 s | 748.8 s |
