# Each storage brings its own traversal, and the sampler keeps one kernel

ADR-0005 changed the model. This record changes the shape of the code that runs it, and it exists because that shape looks wrong from the outside: the storage seam promised one sampler over one concept, and what it delivers is **two traversals**. A reader who finds the dense path indexing a vector while the sparse path materialises a window will want to know whether that is a design or an accident.

It is a design. The dense and the sparse storage want opposite things from a traversal. Dense wants to index its state vector: a read is an address computation, and forcing it through a cache adds a copy, an existence test and a deferred-write branch it never takes. Sparse wants to walk a row once and answer every question from the result: point lookups cost a search, so paying that search per cell is the difference between usable and not. One traversal cannot serve both. The current arrangement — a current-state class that fills from both storages and hands out a dense array — serves the sparse case and taxes the dense one, which is why the dense path carries machinery whose only purpose is to be ignored.

So the traversal moves next to the storage, and the arithmetic stays shared.

## The window

Every storage opens a **window**: a strided view over itself, given a start, a count and a stride. The dense window holds a base and a stride, and a read or a write through it compiles to an indexing of the state vector — no buffer, no fill, no copy. The sparse window materialises on open, exactly as the current clique fill does, writes in place where it already holds the cell, buffers the inserts, and flushes them when it closes.

One abstraction covers both shapes the sampler needs. The block update walks a row: stride 1, one leaf pair after another. A clique walks a column of its tree's node state: stride equal to the leaf count of the other dimension. The existing clique fill is already parameterised exactly this way, which is why this is a generalisation of something that works rather than a new mechanism.

What crosses into the shared kernel is the probability arithmetic and nothing else: the read scalars and one uniform go in, the drawn states and the sufficient-statistic bucket come out. That boundary is where it is because the incremental counter arithmetic is the one thing that must not be written twice. A sign error in one traversal's bucket delta is invisible to every unit test and surfaces only as a wrong error probability after a long chain.

## The readback contract

A sparse window must show its own buffered write to a later read on the same window.

This is not a convenience. A node-state walk goes in post-order, so it reaches a parent after its children and reads the states they were just given. On dense storage those are simply there. On sparse storage a write to a cell the matrix does not yet hold is deferred out of the parallel region, and a naive read would return the old value — so the two backends would compute different chains inside a single update, which is exactly the divergence the byte-identical gate exists to prevent.

`TCliqueWalkStates` was introduced to paper over this: a vector, seeded from the storages, carrying the states the walk had assigned. Once the window owns the readback, that class is solving a problem that no longer exists, and it is deleted. Two mechanisms for one invariant is how the invariant eventually gets stated two different ways.

## The window's lifetime

A clique's window stays open across the alpha, nu and branch-length moves that follow its walk, because those moves have to see the states the walk assigned.

This is longer than the lifetime a reader would assume from the block update, where a window is opened for one row and closed at the end of it. It is written down here because it is the one place where a window is not scoped to the loop that created it, and because closing it early would be a silent correctness bug rather than a compile error: the moves would read the storage, and on the sparse backend they would read stale values for every deferred cell.

## Selecting a backend

The field and the node state are chosen **independently**, by editing one alias each, in one header. The cmake option that chose one backend for both, and the per-backend build directories it drove, are removed.

All four combinations compile. The one the science wants is a sparse field against a dense node state, and the reason is **fill, not size**. Under the AND a field cell is one only where both trees agree, so the field is sparser than either tree field (ADR-0005) and sparse storage pays. A node state, by contrast, sits near its clique's stationary rate — for a mid-range alpha, roughly half ones — where sparse storage costs more per cell than dense and saves nothing.

Size points the other way, which is worth stating so that nobody reaches for the opposite pairing on a memory argument. The node state's own dimension spans every node, not every leaf, so with `n_nodes` near `2 * n_leaves` each node state holds about twice the field's cells. The large object is the one being stored densely; it is still the right call, because a dense byte beats a sparse entry whenever the cell is as likely to be one as not. That pairing could not be expressed at all under one shared option.

Two of the four are gated in continuous integration — dense against dense, and sparse against sparse. Each storage's traversal is exercised by one of those two, so a mixed build cannot diverge in a way the pure builds do not. The other two compile and are covered in-process by the block update's tests, which are instantiated over all four. The header records which pairs are gated, so that the ungated ones are not read as unsupported.

## Why the gate now depends on the random stream

Two traversals that visit cells in different orders can only produce identical output if a cell's uniform depends on the cell and not on the order it was reached. Under a shared mutable generator they would draw in different sequences, and the byte-identical gate could never be green again.

So the position-hashed stream — a uniform derived by hashing the seed, the iteration and the linear cell index — stops being only a reproducibility feature and becomes the load-bearing wall under the gate. It has to land **before** the traversal work, not after, and anything that moves the sampler back onto a shared generator takes the gate down with it.

The gate matters more after this change than before it, not less. Before, the two backends ran the same traversal over different containers, and the ways they could disagree were few and local. Now they run different code, and the gate is the only thing asserting that different code computes one chain.

## A note on the closed-form log link probabilities

Relocated here from `TFieldMath.h`, which now keeps the result and points at this section.

`log_prob_for_bucket` returns `{ log(1 - P_k), log P_k }` from closed forms chosen so that nothing cancels. Taking logs of `prob_for_bucket` would be shorter, and was what the code did first. It is worse in both directions, and not where one would guess.

`1 - P_k` is *exact* whenever `P_k >= 0.5` by Sterbenz's lemma, so `log1p` buys nothing for the both-ones bucket. The loss is at the other end, where `P_0 = omega^2` is tiny and `1 - P_0` rounds to 1. Separately, `log P_2` formed as `log((1 - w) * (1 - w))` inherits the squaring's rounding — 4e-11 relative at `omega = 1e-6` — where ADR-0005's observation that `log P_k` is *affine* in `k` gives it exactly.

So the logs come from the affine form, and the complements from expansions that subtract nothing near-equal:

```
1 - P_0 = (1 - w)(1 + w)
1 - P_1 = 1 - w + w^2
1 - P_2 = w(2 - w)
```

## Considered options

**One traversal on the shared concept, materialising a clique.** This is what the code does today, and what the storage spec originally named as the upgrade path. Rejected because dense can only satisfy it by copying its own vector into a buffer, which is the cost this work exists to remove. The abstraction would be shaped entirely by the sparse case and paid for by the dense one.

**One templated update branching internally on the storage.** Rejected because it puts the sparse deferral machinery in front of every dense reader, and because with independent selection the branching is on a *pair* of storages, so it multiplies rather than adds.

**Keeping the current-state class for the clique walk only.** Rejected because the window's readback contract already provides what the class provides. Keeping both leaves one invariant with two owners.

**Landing the traversal work on the trunk as a changes-no-numbers commit,** as the storage phase before it did. Not available: on the trunk the sampler still draws from a shared generator, so changing the traversal changes the draw order and therefore changes numbers. It can only be a no-numbers change after the position-hashed stream, and that stream belongs to the model branch. The two-phase split is recorded as closed rather than kept open as a plan the branch has overtaken.

## Consequences

**The two paths can drift, and no type prevents it.** This is the real cost. Two traversals over one kernel means a change to one loop that should have been made to both is a divergence the compiler will not catch. The gate is the whole defence, and the gate is only possible because of the random stream, which is why that dependency is written down above rather than left to be rediscovered.

**`is_stored` and the clique fill leave the shared interface.** The first becomes a question only the sparse window asks itself; the second becomes the sparse window's constructor. Both were on the interface only because the sampler had to know which storage it was talking to. It no longer does.

**The benchmark for the clique fill is retargeted rather than deleted.** It measures speed, which no correctness seam covers, and window opening is the whole of the sparse path's cost model.

**This record supersedes nothing.** ADR-0005 stands unchanged; it decides the model, and this decides the code that runs it. ADR-0004's canonical ordering is what makes a leaf pair land at the same `(row, column)` in the field and in either node state, so the windows the block update opens need no index conversion — this record depends on it rather than revising it.
