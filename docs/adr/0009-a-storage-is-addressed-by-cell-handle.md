# A storage is addressed by cell handle, and no window owns a lifetime

ADR-0006 argued that the dense and the sparse storage want opposite things from a traversal, and it made the **window** the abstraction that lets each bring its own. The argument rested on one premise: a sparse point lookup costs a search of a line, so a run of cells has to amortise one. **That premise is gone.** Both sparse storages hold their cells in a hash map keyed by the linear index, a lookup is a hash, and the two backends traverse alike.

So a storage is addressed one cell at a time, and the window is deleted. **This record supersedes ADR-0006 outright.** It also carries the hash-map decision, because the two are one decision: the map is what makes a single cell addressable, and a single addressable cell is what makes the window unnecessary.

## The decision

**A storage hands out a handle.** An updater has four questions to ask before it writes a cell — what state is the cell in, does the storage hold it, what is its linear index, and where is it — and `locate` answers all four at once. `IsOneResult` carries the answers, and its `in_container` flag is true exactly when its pointer is not null. `src/storages/cell_handle.h`.

The handle carries a **raw pointer** and not an iterator. A hash map iterates over key-value pairs, and a caller that dereferenced one would have to know which half of the pair holds the state. A pointer says "the cell is here" under every storage.

**`write_or_defer` is the whole write, and it is one branch.** A cell the storage holds is written where it lies. A cell the storage does not hold, written to a one, waits in a list of deferred inserts. A cell the storage does not hold, written to a zero, is left out: it already reads as zero, so storing it would store a cell to say what its absence says. The helper takes the handle and not the storage, so one body serves every backend and it is a template only over the cell the handle points at.

`open_window`, the window concept, the two window types, `close`, `take_buffered_inserts`, the buffered-write protocol, the two exit paths and the row-versus-column walk a window had to decide on every open are all deleted. So are `LocatableStorage` and `write_state_if_held`: every storage answers `locate` now, so `locate` joins `BinaryStorage` and there is no second concept to select on.

## The hash map, which is the same decision

The sparse field and the sparse node state kept their cells in a sorted-vector matrix. That matrix keeps every cell **twice**, once in its row and once in its column. Two things follow, and the first is why ADR-0006 could not have been written differently at the time.

**It has no single cell to point at.** A handle would have to name one of the two copies, which makes the other stale. This is why the window came first and the handle second: the container the window was built over could not hand one out.

**And a point lookup searched a line.** That is ADR-0006's premise verbatim, and it was true of that container rather than of sparse storage as such.

The two copies could also part company, and the measurement of that is what settles the matter. `cleanUp` drops a cell whose value is the default one, and it leaves a line of one entry alone. The window healed the asymmetry, because it read one copy and wrote both. A write through `TSparseMatrix::update` writes each copy where it finds it, and the asymmetry then grows until `is_one` answers from whichever copy is shorter. **That sent the parity gate apart at the thirty-fifth iteration.** So the window was load-bearing for the matrix, and for nothing else — which is a reason to replace the matrix, not a reason to keep the window.

**A hash map keyed by the linear index holds the cells instead.** A position the map does not hold reads as state 0, so `is_one` is total over the container space and memory tracks the number of ones rather than the size of the space. A stored cell carries a state and not merely a key, so a cell that goes to zero is a **write** and not an erase — which is what lets an update inside a parallel region write it in place. `remove_zeros` reclaims such a cell between iterations, where nothing else is running.

**One array body serves each side.** `TSparseStorage` is the sparse mirror of `TDenseStorage`, and both are templates over the cell they hold, so the map, the handle, the index arithmetic and the ones cursor are each written once. Everything in both is written against `cell_is_one` and `write_state`, so neither body knows whether it holds a bare state or a state packed with a counter.

**The fill argument is ADR-0006's, sharpened, and it points the same way.** A stored cell costs a map node and a bucket slot, which is some tens of bytes where the cell itself is one or two. So the sparse form wins on memory well below **one one in twenty cells**, and not merely below one in two. ADR-0006's rule — choose on the fill of the container and not on its size — stands; the threshold a reader would have inferred from it was wrong by an order of magnitude. The pairing the science wants is unchanged: a sparse field against a dense node state, because under the AND the field is sparser than either tree field (ADR-0005) while a node state sits near its clique's stationary rate.

## The two invariants a window's lifetime enforced structurally

A window had a lifetime, and two properties were true because of that lifetime rather than because anyone stated them. A handle is a pointer with no lifetime of its own, so both now have to be written down and checked.

### A handle does not survive a restructure

The pointer is valid until the next write that can restructure the storage. **A bulk insert or a `remove_zeros` between the `locate` and the write invalidates it.** The convention is stated in `storages/cell_handle.h`, and `initialize_dimensions` clears the handle flag with the cells, because dropping every cell invalidates every handle.

A window could not state this convention, because there was no pointer for a caller to keep: the window *was* the lifetime, and it ended in one of two ways it controlled itself. That is the property being given up, and it is the real cost of this record.

The conformance suite states the convention and shows it is load-bearing rather than decorative. Across each of the two restructuring calls, a handle taken before answers something that is no longer true, and a handle taken after is right again. That a stale handle gets *used* is what the convention forbids, and no test can assert that.

### The ones cursor's cache makes a const call mutate, so it is not thread-safe

A hash map has no order of its own, and the merge joins that read a storage need one. So the cursor walks a sorted vector of the ones, in ascending linear-index order. The vector is `mutable`, and `ones_cursor() const` rebuilds it whenever it is stale.

Three rules follow from that one line, and all three used to be somebody else's problem:

- **A const call writes.** Two threads must not call `ones_cursor()` on one array at once. No caller does, because the merge joins run outside the parallel regions.
- **One cursor at a time per array.** A second `ones_cursor()` may rebuild the very vector the first one walks. A merge join takes one cursor from each of *two* storages, which is what every caller does.
- **A handle defeats the staleness flag.** A write through a handle does not pass through the array, so from the first `locate` the array can no longer see when its ones changed. The dirty flag alone cannot carry this: a cursor taken between the `locate` and the write would clear the flag, and the write would then leave a cache the array believes in.

So a storage that has ever handed out a handle rebuilds the cache on **every** call. That is what the field needs, because the field is written through handles every iteration. An immutable observation locates nothing and still sorts once, which is what the plain dirty flag is for — the LOTUS records sort on their first cursor and never again.

The flag recording that a handle was handed out is sticky: it is set once and never cleared. It is **atomic**, because an update locates a cell from every thread, and a plain `bool` written from two threads is a data race however harmless the value looks. The store is skipped once the flag is set, so after the first one every thread only reads a line nothing writes again, and a relaxed order is enough: what the flag guards is read after the parallel region has joined. `TStickyFlag`, `src/storages/TSparse.h`.

Every mutator sets the staleness flag, `remove_zeros` included, and that costs a sort that finds the same ones. It is kept because every mutator setting the flag is **one** rule. An exemption per mutator is a rule per mutator, and the next mutator added gets it wrong.

## The readback the window owed

ADR-0006's readback contract was the third thing a window's lifetime provided, and unlike the two above it had an owner ready to take it.

A node-state walk goes in post-order, so it reaches a parent after its children and reads the states they were just given. On a sparse node state a write to an absent cell is deferred, and a naive read would return the old value — so the two backends would compute different chains inside one update, which is exactly what the parity gate exists to prevent.

`write_or_defer` answers whether the write landed in the storage, and the caller keeps what it could not place. `TCliqueView` is where the clique's walk keeps it. A caller that reads a cell back before the deferred list is committed needs that answer; every other caller ignores it.

## What ADR-0006 decided that this record keeps

Four of ADR-0006's decisions had nothing to do with the traversal, and they carry over unchanged. They are restated here because a reader sent from a code comment to a superseded record needs to know which half of it to believe.

**No write inside a parallel region may insert.** ADR-0006's argument is untouched: a clique owns its column and shares every row, and splitting the work the other way only swaps which of the two is shared, so no partition of a matrix keeps both private. Under a hash map the reason has the same shape — an insert rehashes the map while other threads read it — and the mechanism is the same. Each caller hands out a list of linear indices, and one bulk insert runs after the region. A dense caller hands out an empty list and the dense bulk insert does nothing, so one loop body serves both backends.

**The field and the node state are each selected by an alias of their own**, by editing one line in one header. All four pairings compile, and continuous integration gates two of them -- sparse against sparse, and dense against dense, which is the default. Those two between them exercise both storages, and the header records the table, so an ungated pair reads as untested rather than as unsupported. Two storages make no choice: the simple error model data follows the field, because the two are read cell for cell against each other, and the LOTUS records are pinned sparse whatever the build selects. Each would otherwise be a third pairing for the gate to cover.

**The arithmetic stays in one shared kernel.** Still the rule, and cheaper to keep now: with one traversal there is no second loop for it to drift from. ADR-0006 named that drift as its own headline cost — two paths, and no type preventing a change to one from missing the other. **This record pays that cost off.**

**The parity gate still rests on ADR-0007's stream, and less tightly than it did.** ADR-0006 made the stream load-bearing under the gate: two traversals visit cells in different orders, and a shared generator hands the same cell different uniforms. One traversal removes that particular dependency, because the gate pins `--numThreads 1` and both backends now visit one order. The stream stays for ADR-0007's own reasons -- a run at more than one thread is otherwise not reproducible at all -- and the gate is still the only thing asserting that the two backends compute one chain.

One more part of ADR-0006 is neither superseded nor restated. Its section on the closed-form log link probabilities is a numerical-analysis result about `log_prob_for_bucket`, and it was never about the traversal. `TFieldMath.h` points at that section, and it stays there.

The code cites ADR-0006 in about a dozen comments. Those citations are left alone rather than swept: ADR-0006 carries a banner that names this record, so a reader who follows one lands here, and a rewrite of a dozen comments is a large diff over a pointer that already works.

## The dense field's counter width, and the thinning factor

A field cell packs its posterior counter into the same 16-bit word as its state, which leaves the counter **15 bits**. The dense field used to keep its counters in an array of their own, a full `uint16_t` each, so it had a sixteenth bit the sparse cell did not have. The two backends therefore thinned one chain differently: `ceil(n / 65535)` against `ceil(n / 32767)`.

**Both hold the packed cell now, and both count one iteration in `ceil(n / 32767)`.**

One consequence is wanted, and it is why the change was made. The resolution is the **cell's** and no longer the backend's, so a posterior field written by one backend compares with a posterior field written by the other, and the parity gate no longer refuses a long chain — the two write traces of the same length whatever the chain. The thinning factor also decides which iterations get a trace line, and it is the resolution of the posterior field, of both tree field posteriors and of the field's trace at once, so one number had to describe both backends.

One consequence is paid. **A dense run of more than 32 767 iterations writes a posterior field at half the previous resolution, and half as many trace lines.** Nothing else a run produces changes. Every run reports its thinning factor to its log file, and `README.md` records what changed, so a reader comparing a new posterior field against an old one is told why the line counts differ.

The dense field gives up a bit it could have kept. That is deliberate: a resolution that depends on which alias the build selected is a number no reader can quote.

## Considered options

**Give the sorted-vector matrix a cell to point at.** Rejected. It holds every cell twice, so a handle has to name one copy and leave the other stale — and the two already drift under `cleanUp`, which is the parity failure above. Repairing the drift *and* adding a handle is more work than replacing the container, and it leaves a container whose point lookup is still a search of a line.

**Keep the window for the sparse backend and address the dense one by cell.** Rejected. That is ADR-0006's two paths with the reason for them removed, so it keeps the drift ADR-0006 named as its own headline cost and buys nothing back. It also keeps the sparse deferral machinery in front of a caller that would rather not know.

**Keep the window over the hash map, as a batching device.** Rejected. With a lookup at hash cost there is no search left to amortise, so a window buys the amortisation of nothing. It still costs a lifetime, two exit paths, and a walk direction to decide on every open.

**Erase a cell that goes to zero, rather than storing a zero.** Rejected. An erase restructures the map while other threads read it, so a one-to-zero transition would have to be deferred as well — and that is the common transition, where a zero-to-one on an absent cell is the rare one. Storing the zero keeps every write inside a parallel region in place except the rare one, and `remove_zeros` reclaims between iterations, where nothing else runs.

## Consequences

**A handle is a raw pointer, and the compiler does not check its lifetime.** This is what the window absorbed and what a reader is now trusted with. The convention is one sentence, the conformance suite shows it is load-bearing, and an invalidated handle is the thing to suspect first when a sparse run and a dense run part company.

**One traversal, so ADR-0006's drift is gone.** The two backends run the same loop over different containers. A change to that loop cannot be made to one path and missed on the other, because there is one path.

**The benchmark is retargeted rather than deleted.** It measured window opening, which was the whole of the sparse path's cost model. It measures the point lookups that replace it, which is what a hash-map backing has to beat.

**The parity gate's fixture rationale changed with the traversal.** Its non-square shapes now catch a wrong clique-to-cell mapping rather than a wrong row-versus-column walk, which is the same bug class one layer up.

**This record supersedes ADR-0006 and nothing else.** ADR-0004's canonical ordering is still what makes a leaf pair land at the same `(row, column)` in the field and in either node state, so a caller reaching a cell needs no index conversion. ADR-0005 decides the model and stands untouched. ADR-0007's stream stands, for the reasons that record gives rather than for the one ADR-0006 gave it. ADR-0008 changes who draws the leaf layer, which is a caller of this interface and not a part of it.

**CONTEXT.md keeps *Window* as a retired term.** It names what ADR-0006 built and this record retired, so the word can still be resolved where it survives in an old comment or an old commit.
