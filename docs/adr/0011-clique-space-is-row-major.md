# Clique space has one convention, and it is row-major

A clique of a tree carries a leaf in every dimension but that tree's own, which carries a 0 (CONTEXT.md, "Clique"). Numbering the cliques means fixing a correspondence between a clique number and that subscript, and the code held two of them. `TTree::_clique_index` read a clique number as a **row-major** subscript, through `coretools::getSubscriptsAsArray`. `TTree::transition_grid_of_cell` walked a hand-rolled **column-major** stride. The two are inverses of each other and they agreed, so nothing failed.

They agreed by coincidence. Clique space has exactly one dimension above 1 when there are two trees — the owning tree's carries a 1 — and with one non-unit dimension every linearisation lands on the same number. The header said so, in a comment that ended "A third tree would need one convention here." That comment was the only thing holding the invariant, and a comment cannot be asked at build time.

`TCliqueSpace` (`src/tree/clique/TCliqueSpace.h`) now holds the correspondence in both directions, row-major. The column-major loop is deleted rather than ported.

Row-major wins because the rest of the codebase already reads that way. `TDenseStorage` and `TSparseStorage` both linearise through `coretools::getLinearIndex` and `getSubscriptsAsArray`, so the field, both node states and the observed data all agree on what a subscript means. The column-major stride had one author, one call path, and no second example in `src/`. Choosing it would have made clique space the one space in the model that reads backwards.

## Considered options

**Keeping both and asserting they agree.** This states the coincidence instead of removing it, and the assertion would pass on every configuration that can be built today. It buys nothing a comment did not already buy.

**Column-major everywhere.** It would require changing all four storages and `node_state_shape.h`, which is a large diff to make the model read against its own grain.

**No template parameter, and the divergence recorded as prose.** `IndexArray` is `std::array<size_t, 2>`, so a clique space that stores one cannot be handed three dimensions, and the discriminating case becomes a worked example in this record rather than a test. It was rejected because the reason for the change then stops being executable: a future edit that reintroduced a second convention would break nothing.

## Evidence

The two conventions were compared directly over the leaf counts `{7, 1, 3}` with the middle dimension owned. At two dimensions they disagree on **0 of 7** cells. At three they disagree on **18 of 21**. The comment's claim was exact.

`just parity` was run at `cacf9ca`, before `TTree` moved onto the module, and its output kept. It was run again after the move and again after `TBlockModel` was handed the whole leaf pair. All 58 non-log files — both backends, simulate and infer — are byte-identical across the three runs. The refactor moves no chain.

Note what that comparison is and is not. The parity gate's own check gives two backends the same work and requires the same bytes, so it cannot see an error symmetric across both. The before-and-after comparison above is a different one: the same backend against its own earlier output, at a fixed seed and one thread. It has no such blind spot, and it is the check a behaviour-preserving refactor of this code should be held to.

## Consequences

`TCliqueSpace` carries the dimension count as a template parameter, `NumDim`, defaulting to `NUMBER_OF_TREES`. **This is a test affordance and not a seam.** Production instantiates it once. Nothing varies across the parameter at run time, and no second instantiation ships. It exists so that `TCliqueSpace_Tests` can build the three-dimension case that tells the two conventions apart — the case the shipping configuration cannot produce.

It should not be read as a half-built N-tree feature and finished, and it should not be deleted as an unused abstraction. The model is two-tree in the places that matter: `block_update` carries a `static_assert`, the eight-state table names two tree fields, and `IndexArray` is two long. A third tree is a change to those, and `TCliqueSpace` is deliberately not one of the places that would have to be argued about when it happens.

The rule "one adapter means a hypothetical seam" — which ADR-0009 used to delete the window — does not reach this. A seam is a place behaviour varies. `NumDim` is a type parameter serving a test, and the test is the reason the module exists at all.
