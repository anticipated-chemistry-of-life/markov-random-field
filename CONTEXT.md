# Metabolite Inference

A Markov random field over which molecules occur in which species. Two phylogenies — one of species, one of molecules — jointly constrain a latent binary presence matrix, which is observed only indirectly through noisy literature and assay data.

## The trees

**Tree**:
One of the two phylogenies the field is defined over, one of species and one of molecules. A tree may have more than one root; each root's subtree is drawn independently, with the root's own state taken from the stationary distribution rather than from a parent.
_Avoid_: forest, taxonomy, hierarchy

**Node**:
One vertex of a tree.
_Avoid_: vertex, taxon

**Leaf**:
A node with no children. Leaves are the only nodes the observed data may reference: the field, LOTUS records and simple error model data are all indexed in leaf space.
_Avoid_: tip, terminal

**Root**:
A node with no parent. A tree may have several. A root has no branch, and therefore no branch length; its state comes from the stationary distribution instead of from a parent.
_Avoid_: ancestor, origin

**Internal node**:
A node with at least one child, roots included. The distinction is about tree structure only: every node carries a node state, leaves included.
_Avoid_: ancestral node, non-leaf

**Branch**:
The edge from a node to its parent, identified by that child node. Every node except a root has exactly one, so a tree has `n_nodes - n_roots` branches — the count the branch-length budget is built from.
_Avoid_: edge, leaves and internal nodes without roots

## The latent field

**Field**:
The latent binary matrix recording, for every (species leaf, molecule leaf) pair, whether that molecule occurs in that species. Written `Y`. It is the *reconciled* field: a noisy AND of the two tree fields, and not a variable the two trees share. See ADR-0005.
_Avoid_: presence matrix, occurrence matrix, Y-space, shared field (except when naming the model ADR-0005 retired)

**Node state**:
The latent binary state of a node of one tree, replicated across every clique of that tree. Written `Z`. Every node carries one, leaves included: the species tree's node state is `n_nodes(species) x n_leaves(molecules)` and the molecule tree's is `n_leaves(species) x n_nodes(molecules)`.
_Avoid_: internal state, ancestral state, hidden state

**Tree field**:
The leaf block of one tree's node state — that tree's own view of the leaf-level field, before the two are reconciled. Written `Z_s` and `Z_m`, the leaf block of the node state written with the same letter. There is one per tree, and together they are what the field reconciles. A tree field and the field are addressed at the same `(row, column)` for a given leaf pair, so the correspondence between them is the identity rather than a conversion. See ADR-0005.
_Avoid_: per-tree field, own field, leaf state, Z at the leaves

**Error probability**:
The probability that a tree field cell is corrupted before the two tree fields are reconciled into the field. Written `omega`. One scalar, shared by both trees, constrained to the open interval `(0, 0.5)`. It is estimated, under an exponential prior truncated to that interval, and its Metropolis move reads the link counters rather than the cells. See ADR-0005.
_Avoid_: noise rate, flip probability, error rate (too easily confused with the simple error model's misreport probability)

**Bucket**:
The number of tree fields in state 1 at one leaf pair, so 0, 1 or 2. The link table depends on the two tree field states only through it, which is what collapses four cells to three probabilities. See ADR-0005.
_Avoid_: class, category, sum of the tree fields

**Link counters**:
The link's sufficient statistic: six integers, `n(bucket, field state)`, counting the leaf pairs of the whole field. The link's whole likelihood is a function of them and the error probability, so the error probability's move costs the same whatever the size of the field. The block update tallies them as it goes. Traced to `<prefix>_link_counters_trace.txt`, which is what the AND diagnostic reads. See ADR-0005.
_Avoid_: sufficient statistics, the six counts, contingency table

**Clique**:
A set of nodes that vary along exactly one tree's dimension while every other dimension is fixed at a leaf. Cliques belong to a tree: a species-tree clique is identified by a *molecule* leaf, and vice versa.
_Avoid_: slice, column, replicate

**Alpha**:
The stationary probability that a node in a given clique is in state 1. One value per clique.
_Avoid_: prevalence, base rate, pi

**Nu**:
The switching rate of the two-state continuous-time process running along a tree's branches. One value per clique, carried in log space.
_Avoid_: rate, mu, lambda

**Transition grid**:
One clique's two-state process discretised onto the bin grid: one transition matrix per bin, plus the stationary distribution its roots are drawn from. Built from an alpha, a nu and a bin grid, and immutable — a Metropolis proposal builds a second grid rather than mutating the first, so there is no "try" state. `TTransitionGrid`, `src/tree/branch/`.
_Avoid_: transition matrices, lambda matrices, clique process, try matrix

**Neutral dimension**:
A tree dimension whose parameters are pinned so that every transition matrix row is exactly (0.5, 0.5), making that tree's node state an independent coin flip at every node. A neutral dimension carries no phylogenetic signal. It does still *influence* the field under ADR-0005 — its tree field is half ones, so the field's expected density is exactly half the active tree's corrupted rate — where under the old shared-field model its factor cancelled from every field conditional instead.
_Avoid_: disabled tree, ignored dimension, flat tree

**Field normalising constant**:
Whatever a field distribution would have to be divided by to make it a proper density. Written `C`. Under this model it is identically 1, at every parameter value: the field has exactly one conditional density, `p(Y | Z_s, Z_m, omega)`, and it sums to 1 over fields by construction, so there is nothing left to normalise. The term survives only to name what ADR-0002 diagnosed and ADR-0005 removed — under the old shared-field model `C` was the sum, over every possible field, of the product of the two trees' field likelihoods, and it moved with both trees' parameters, biasing them toward small nu. See ADR-0005.
_Avoid_: partition function, Z (that is the node state), evidence

## The chain

**Update**:
One full pass over a set of variables. The block update visits every leaf pair; a tree's node-state update visits every internal node of every clique of that tree. An *iteration* is one turn of the whole chain, and holds several updates in a fixed order: block-update every leaf pair, update each tree's internal node state, then the parameters.
_Avoid_: sweep, pass, scan

**Chain start**:
The configuration a chain holds before its first update. Both tree fields start at one wherever a LOTUS record exists, and zero elsewhere. The field starts matching them. Each tree then initialises every internal node from its children, in one forward pass. Every state is a mode and not a draw. Under the AND a record is strong evidence that both tree fields are one at that cell, so the chain starts near the posterior mode. With no record anywhere the start is all zeros. `leaf_layer_start` in `src/field/` starts the leaf layer, and `TTree::initialize_Z_from_children` the nodes above it. See ADR-0005.
_Avoid_: initial values, seed, guess, warm-up

**Block update**:
The joint draw over the field and both tree fields at one leaf pair, taken from all eight combinations at once rather than one variable at a time. Exact, not an approximation. It is what escapes the state the AND makes metastable: with a small error probability a field cell at one pins both tree fields to one, and single-variable draws can only escape through the field. See ADR-0005.
_Avoid_: Y update, field update, joint draw, eight-state sweep

**Cell uniform**:
The one uniform a cell's update draws, derived by hashing the seed, the stream, the tree, the iteration and the cell's linear index instead of taken from a running generator. Two cells, two iterations, two containers and two seeds share one only by chance, and the number a cell gets does not move when the thread count changes or when an update visits the cells in another order. That last property is what lets the dense and the sparse backend traverse their storage differently and still run one chain. `TCellUniforms`, `src/random/`. See ADR-0007.
_Avoid_: random number, uniform variate, the cell's random draw

**Window**:
The strided view a storage opens over itself, given a start, a count and a stride. An update reads and writes its cells through a window rather than through a cache of its own, so each storage brings the traversal that suits it: the dense window indexes the state vector, and the sparse window walks its line once and buffers the inserts it cannot make in place. A window shows its own write to a later read on the same window. It ends once: it either commits its buffered inserts or hands them to the caller, which is the only exit open to a window inside a parallel region. `TDenseWindow` and `TSparseWindow`, `src/storages/`. See ADR-0006.
_Avoid_: slice, view, buffer, current state

## Branch lengths

**Bin**:
The discrete index, in `0 … n_bins-1`, standing for a branch's length. Branch lengths are never continuous in the model; they are only ever bins.
_Avoid_: branch index, discretised length

**Grid branch length**:
The continuous length a bin represents, taken at the bin's centre.
_Avoid_: bin midpoint, branch length

**Branch-length budget**:
The total of a tree's bins, which is fixed at `n_branches · n_bins / 2` and conserved for the lifetime of a chain. Equivalent to requiring the mean grid branch length to be exactly 1. A set of branch lengths that misses the budget is unreachable, not merely improbable.
_Avoid_: branch length sum, normalisation constraint

**Bin grid**:
One tree's bin↔length correspondence: the bin width, the grid branch length each bin stands for, the branch-length budget their sum must hit, and the ±1 step that conserves it. A pure function of `n_bins` — it knows nothing of the tree's topology, of the parameters, or of the random generator, which is what lets it be tested without running a chain. `TBinGrid`, `src/tree/branch/`.
_Avoid_: branch-length grid, binning, discretisation

## Observations

**LOTUS record**:
An observed (species, molecule) occurrence reported in the literature. Written `L`. Absence of a record is uninformative wherever research effort is low.
_Avoid_: citation, literature record, observation

**Research effort**:
How thoroughly a (species, molecule) pair has been looked for. Determines whether a missing LOTUS record means "absent" or "unstudied". One factor per tree, multiplied together: `1 - exp(-gamma_i · log(count_i + 1))`, where `count_i` is the paper count of the leaf the pair occupies in dimension `i`. Note there is **one gamma per tree**, not one overall — the independent reference simulates with a single scalar for both trees, which is a special case, not the model.
_Avoid_: coverage, sampling effort, detection probability

**Paper count**:
The number of publications covering a leaf. The sole input to research effort.
_Avoid_: occurrence counter, citation count

**Simple error model data**:
An observation of the field in which every cell is independently misreported with a fixed probability. Written `D`. Unlike a LOTUS record, it is dense: every cell is observed.
_Avoid_: noisy Y, simple data, flipped field
