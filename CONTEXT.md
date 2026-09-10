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

**Tree field posterior**:
How often each cell of one tree field was a one, over the counted iterations of a chain. One file per tree, `<prefix>_<tree>_tree_field_posterior.txt`, beside the field's own. The field's own cannot stand in for it: the field's expected density is the product of the two corrupted rates, so it says nothing about how the rate splits between the trees. The two tree fields moving in opposite directions along that ridge is what shows up first. A node state carries no counter, so this one stands beside it and counts the leaf block alone. `--write_tree_field_posteriors` decides whether a run holds one: it is a counter per leaf pair per tree, which a run that chose the sparse field chose not to pay for the field itself. `TTreeFieldPosterior`, `src/field/`. See ADR-0005.
_Avoid_: Z posterior, leaf posterior, per-tree posterior

**Error probability**:
The probability that a tree field cell is corrupted before the two tree fields are reconciled into the field. Written `omega`. One scalar, shared by both trees, constrained to the open interval `(0, 0.5)`. It is estimated, under an exponential prior truncated to that interval, and its Metropolis move reads the link counters rather than the cells. See ADR-0005.
_Avoid_: noise rate, flip probability, error rate (too easily confused with the simple error model's misreport probability)

**Bucket**:
The number of tree fields in state 1 at one leaf pair, so 0, 1 or 2. The link table depends on the two tree field states only through it, which is what collapses four cells to three probabilities. See ADR-0005.
_Avoid_: class, category, sum of the tree fields

**Link counters**:
The link's sufficient statistic: six integers, `n(bucket, field state)`, counting the leaf pairs of the whole field. The link's whole likelihood is a function of them and the error probability, so the error probability's move costs the same whatever the size of the field. The field update retallies them over every leaf pair as it goes. Traced to `<prefix>_link_counters_trace.txt`, which is what the AND diagnostic reads. See ADR-0005.
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
One full pass over a set of variables. A tree's node-state update visits every node of every clique of that tree, leaves included; the field update visits every leaf pair. An *iteration* is one turn of the whole chain, and holds several updates in a fixed order: the species tree's node state, the molecule tree's node state, the field, then the parameters.
_Avoid_: sweep, pass, scan

**Chain start**:
The configuration a chain holds before its first update. Both tree fields start at one wherever a LOTUS record exists, and zero elsewhere. The field starts matching them. Each tree then initialises every internal node from its children, in one forward pass. Every state is a mode and not a draw. Under the AND a record is strong evidence that both tree fields are one at that cell, so the chain starts near the posterior mode. With no record anywhere the start is all zeros. `leaf_layer_start` in `src/field/` starts the leaf layer, and `TTree::initialize_Z_from_children` the nodes above it. See ADR-0005.
_Avoid_: initial values, seed, guess, warm-up

**Field update**:
The field's own pass over its cells. It visits every leaf pair and draws that cell from the two tree field cells at that pair, and from the data that observes the field. It retallies the six link counters as it goes. It is the last state update of an iteration, so the counters describe the configuration the error probability then proposes against. The tree fields are not its to draw: each is drawn by its own tree, as the leaf block of that tree's node state. `field_update::run`, `src/field/`.
_Avoid_: Y update, Y sweep, leaf pair update

**Block update**:
The joint draw over the field and both tree fields at one leaf pair, taken from all eight combinations at once rather than one variable at a time. The term survives only to name what ADR-0005 built and ADR-0008 retired. Each tree now draws its own leaf states, and the field has an update of its own, so the three variables move one at a time. ADR-0005 built the block to escape the state the AND makes metastable: with a small error probability a field cell at one pins both tree fields to one, and single-variable draws can only escape through the field. ADR-0008 records what dropping it costs. See ADR-0008.
_Avoid_: joint draw, eight-state sweep

**Joint density**:
The log density of the whole model at one configuration: `log p(Z_s | theta_s) + log p(Z_m | theta_m) + log p(Y | Z_s, Z_m, omega) + log p(L, D | Y)`. Every factor is a proper conditional density, so the sum is one too — which the sum of the two trees' likelihoods it replaces was not (ADR-0002). It is the no-drift instrument: one number an iteration, and a chain that drifts moves it. A tree's own factor scores each node once, against its parent or against the stationary distribution, so each branch is counted once. Each factor takes a column of `<prefix>_joint_density.txt`. `--write_joint_log_prob_density` decides whether the file is written at all, because the two tree factors cost a pass over every node. See ADR-0005.
_Avoid_: likelihood, posterior, complete joint density

**Cell uniform**:
The one uniform a cell's update draws, derived by hashing the seed, the stream, the tree, the iteration and the cell's linear index instead of taken from a running generator. Two cells, two iterations, two containers and two seeds share one only by chance, and the number a cell gets does not move when the thread count changes or when an update visits the cells in another order. That last property is what lets the dense and the sparse backend traverse their storage differently and still run one chain. `TCellUniforms`, `src/random/`. See ADR-0007.
_Avoid_: random number, uniform variate, the cell's random draw

**Window**:
The strided view a storage opened over itself, given a start, a count and a stride. The term survives only to name what ADR-0006 built and ADR-0009 retired. A storage is addressed one cell at a time now: it hands an updater a pointer to the cell, and the updater writes the cell where it lies or defers the one insert it cannot make in place. ADR-0006 built the window so that each storage could bring the traversal that suits it, because a point lookup on the sorted-vector matrix cost a search of a line. Both sparse storages hold their cells in a hash map now, so a lookup is a hash and the two backends traverse alike. ADR-0009 records the two invariants the window's lifetime used to enforce for free. See ADR-0009.
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
