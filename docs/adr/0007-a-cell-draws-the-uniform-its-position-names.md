# A cell draws the uniform its position names

ADR-0006 records why the dense and the sparse backend traverse their storage differently on purpose, and it names the one thing that lets two traversals still produce one chain: a cell's uniform has to come from the cell, not from the order the cell was reached. This record is that stream.

## The defect

The field update and the node-state update are both parallel, and both drew from `coretools::instances::randomGenerator()`. That generator is `thread_local`, and only the thread that runs `TTask::initializeRandomGenerator` is ever seeded from `--fixedSeed`. A worker thread constructs its own on first use, from the wall clock.

So a run at more than one thread was not reproducible at all — not across machines, not across two runs on the same machine. At one thread it was, which is why every byte-identical claim in the repository is qualified that way, and why `tests/backend_parity/run.sh` pins `--numThreads 1`.

Two consequences follow, and the second is the one that forces the timing. A bug that only appears under parallelism could not be reproduced, so it could not be bisected. And the traversal work of ADR-0006 could not land at all: once the two backends walk the storage in different orders, a shared generator hands the same cell different uniforms, and the gate that says the two backends agree could never be green again.

## The decision

A cell's uniform is a hash of the seed, the stream, the tree, the iteration and the cell's linear index. `TCellUniforms`, in `src/random/`.

The class is a key and a counter. Seed, stream, tree and iteration mix down to one 64-bit key when the object is built, once per update; the cell's linear index is the counter, and one more mix turns key and counter into the cell's uniform. The mixing function is the SplitMix64 finaliser and the counter steps by SplitMix64's own increment, so the cells of one stream walk a SplitMix64 sequence — a generator whose statistical properties are published, rather than a hash invented here.

Nothing in the derivation counts the draws that came before, so a cell gets the same number whichever thread reaches it, and whichever order the update visits the container in. That is the whole property, and everything else in this record follows from it.

## Why a stream label

Two containers hold a cell at linear index 5, and two draws visit the same container.

The field and each tree's node state index their cells in spaces of their own, so without a label the field's cell 5 and the species node state's cell 5 would share a number at every iteration — a correlation between two variables the model says are conditionally independent. The tree is a second component for the same reason: the two trees' node states are two containers.

The draw a simulated chain starts from carries its own label as well. Without it, the state a cell is given before the chain starts and the state the chain's first update gives it would be one draw rather than two.

## Why the two-state draw stopped using the odds-ratio helper

`sample` used `coretools::TAcceptOddsRatio::accept`, which draws its own uniform. Taking a supplied one would need a second overload in coretools, and coretools is checked out by `FetchContent` rather than committed here — an edit there is an edit to somebody else's tree that this repository would lose on the next fetch.

So the draw is written locally, as one comparison of the uniform against `1 / (1 + exp(-log odds))`. That is the probability of state 1, and it costs one exponential where the helper's lookup table often avoided a logarithm.

It is the faster of the two, and by more than the table saves. The helper's cost is dominated not by the logarithm it avoids but by the draw it makes: a Mersenne twister step and a uniform distribution. Twenty million draws in a release build, on one core:

| draw                                                  | per draw |
| ----------------------------------------------------- | -------- |
| the generator plus the odds-ratio table                | 22.2 ns  |
| the position hash plus one exponential                 | 4.0 ns   |

Measured in a tight loop, so it flatters both sides equally and neither carries the update's cache traffic. The point it settles is only that this exchange is not a cost, which is the question a reader would otherwise have to ask.

## What this record does not cover

The alpha, nu and branch-length moves still draw from the thread-local generator.

They run inside the same parallel loop over cliques, through `stattools::TParameter::propose` and `acceptOrReject`, which own their draws the way the odds-ratio helper did. Those are cliques rather than cells, so no linear cell index names them, and the number of uniforms one move consumes depends on the value it is scoring — `evalLogH` draws nothing when the log Hastings ratio is positive. A position-hashed stream cannot be threaded through an interface that does not take one.

The mass-spec assignment moves are the same shape of problem, in `TMSMSData::update_all_MS_assignments`: a parallel region whose moves draw a variable number of uniforms each. They are compiled out unless `USE_MS_DATA` is on, and no run constructs the data yet, so they cost nothing today.

The consequence is worth stating plainly rather than leaving to be discovered. Measured on the parity fixture at 200 iterations, with LOTUS and the simple error model compiled in and mass spec off:

| run        | one thread, twice | four threads, twice | one against four |
| ---------- | ----------------- | ------------------- | ---------------- |
| `simulate` | identical         | identical           | identical        |
| `infer`    | identical         | **differ**          | **differ**       |

`simulate` is reproducible at any thread count because `IsSimulation` compiles the parameter moves out of the clique loop, which leaves only cell draws inside it. `infer` is not, and two runs at four threads differ from each other as well — so for that build the residual sits on the parameter moves, and nowhere in the field or the node state. `run.sh` now runs the `simulate` half of that table on every push, so the green half is gated rather than remembered. A difference in a parameter reaches every cell on the next iteration, which is why it takes the whole run with it. `run.sh` therefore keeps `--numThreads 1`.

Most of the serial draws stay on the generator, and that is not an oversight. Initial values, the branch-length proposals, the pair sampler, and the LOTUS and simple-error data drawn from a simulated field all run on one thread, outside every parallel region, in an order the code fixes. A generator is the right tool there, and the measured table above shows those outputs already identical at four threads.

Two serial draws did move: the field and the node state a simulated chain starts from. They draw one uniform per cell, which is exactly what a stream serves, and `sample` no longer draws one for them. Hashing them costs nothing and takes the ordering assumption out of two loops that a later change could reasonably want to parallelise or reorder.

Fixing the parameter moves is a change to the shape of the clique loop, not to this stream. The proposals would move to a serial pass before the loop, the parallel pass would score every combination the two moves can reach, and the accept-or-reject calls would run serially after it — with the branch-length term, which reads the grid those moves install, moving with them. That is a separate piece of work, and the loop it rewrites is the one #54 rewrites anyway.

## Considered options

**Seed each thread's generator deterministically.** Rejected. `schedule(dynamic)` decides which thread takes which clique, so even identically seeded threads consume their draws in an order that changes between runs. It also does not survive ADR-0006: two traversals still visit cells in different orders within one thread.

**One shared generator behind a lock.** Rejected. It serialises the hot loop and it still fixes nothing: the order threads reach the lock in is the order that varies.

**Reseed the thread-local generator per cell from a position hash.** Rejected. Seeding a Mersenne twister initialises 624 words, which is far more work than the draw it enables, and `setSeed` adds the thread number to the seed it is given — the very dependency being removed.

**Draw a block of uniforms serially and hand them out.** Rejected. It keeps a run reproducible only at a fixed traversal: the block is indexed by visit order, so ADR-0006's second traversal indexes it differently. It also makes the update carry a buffer the size of the container.

## Consequences

**Every chain moves.** The uniforms are different numbers, so no trace from before this change compares with one from after it. Nothing about the target distribution changes, so a run is as valid as it was; it is a different sample from the same posterior.

**The parity gate is now load-bearing in a new way.** ADR-0006 says the gate is the only thing asserting that two traversals compute one chain. This stream is what makes that possible, so anything that moves a cell update back onto a running generator takes the gate down with it.

**`sample` needs a uniform from every caller.** That is deliberate. A caller that cannot name the cell it is drawing for is a caller that has no business drawing, and the compiler now says so.

**The seed has to be read outside a parallel region.** `run_seed` reads the thread-local generator, and a worker thread's copy carries a seed of its own. Every construction of a stream therefore sits before the region that uses it, which is also where it belongs: the key is built once per update, not once per cell.
