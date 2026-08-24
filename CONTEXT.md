# Metabolite Inference

A Markov random field over which molecules occur in which species. Two phylogenies — one of species, one of molecules — jointly constrain a latent binary presence matrix, which is observed only indirectly through noisy literature and assay data.

## The latent field

**Field**:
The latent binary matrix recording, for every (species leaf, molecule leaf) pair, whether that molecule occurs in that species. Written `Y`.
_Avoid_: presence matrix, occurrence matrix, Y-space

**Internal state**:
The latent binary state of a non-leaf node of one tree, replicated across every clique of that tree. Written `Z`.
_Avoid_: ancestral state, hidden state

**Clique**:
A set of nodes that vary along exactly one tree's dimension while every other dimension is fixed at a leaf. Cliques belong to a tree: a species-tree clique is identified by a *molecule* leaf, and vice versa.
_Avoid_: slice, column, replicate

**Alpha**:
The stationary probability that a node in a given clique is in state 1. One value per clique.
_Avoid_: prevalence, base rate, pi

**Nu**:
The switching rate of the two-state continuous-time process running along a tree's branches. One value per clique, carried in log space.
_Avoid_: rate, mu, lambda

**Neutral dimension**:
A tree dimension whose parameters are pinned so that every transition matrix row is exactly (0.5, 0.5), making that tree contribute a constant factor to every field conditional. A neutral dimension cannot influence the field.
_Avoid_: disabled tree, ignored dimension, flat tree

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

## Observations

**LOTUS record**:
An observed (species, molecule) occurrence reported in the literature. Written `L`. Absence of a record is uninformative wherever research effort is low.
_Avoid_: citation, literature record, observation

**Research effort**:
How thoroughly a (species, molecule) pair has been looked for, derived from how many papers cover each of the two leaves. Determines whether a missing LOTUS record means "absent" or "unstudied".
_Avoid_: coverage, sampling effort, detection probability

**Paper count**:
The number of publications covering a leaf. The sole input to research effort.
_Avoid_: occurrence counter, citation count

**Simple error model data**:
An observation of the field in which every cell is independently misreported with a fixed probability. Written `D`. Unlike a LOTUS record, it is dense: every cell is observed.
_Avoid_: noisy Y, simple data, flipped field
