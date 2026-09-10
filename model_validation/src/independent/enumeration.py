"""A model small enough to enumerate, to assert that the joint sums to one.

ADR-0005 replaces the product of two tree likelihoods with a directed
factorisation:

    p(Z_s | theta_s) * p(Z_m | theta_m) * p(Y | Z_s, Z_m, omega)

Every factor is a proper conditional density, so the joint sums to one at every
parameter value. That is an argument. This module is the evidence. It shrinks
both trees until every `(Z_s, Z_m, Y)` can be listed. It then adds the joint up
over all of them and hands back the total.

The same enumeration computed ADR-0002's `C` and reported how far it moved with
the parameters. The new model leaves nothing there to measure, so it asserts
instead.

Read the total for what it is. The link hands each leaf pair one Bernoulli
probability, so summing the field out of a *correct* factorisation is an
identity. The total therefore reads two things: that each tree's node-state
density is proper, and that no variable is scored twice. The second is ADR-0002's
defect. That is why the sum runs over the product space, and not over one tree at
a time.

The total cannot see the shape of the link, the error probability or the
orientation of a leaf block. The field rate and the reference draw pin those, and
`test_independent` holds both.

Both trees stay active. ADR-0005 says why a neutral tree would hide the defect
this model repairs.

The maths comes from ADR-0005. It shares no code with `src/field/`, which is the
implementation this harness exists to disagree with.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import link as L
from .field import N_BINS, grid_branch_lengths, transition_matrices
from .indexing import TreeIndex

MAX_ENUMERABLE_STATES = 2**22
"""How many configurations one enumeration may list.

A tree grows as `2 ** (nodes * cliques)`, so a handful of leaves reaches this
bound. The sum over every configuration is the intractable object. Shrink the
trees rather than raise the bound.
"""


@dataclass(frozen=True)
class TreeProcess:
    """One tree's process: its topology, its branch bins, and its parameters."""

    index: TreeIndex
    bins: np.ndarray
    """One bin per branch, in branch order."""

    alphas: np.ndarray
    """One per clique, and a clique of this tree is named by a leaf of the other."""

    nus: np.ndarray
    """One per clique. The nu itself, not its logarithm."""

    def __post_init__(self) -> None:
        if len(self.bins) != self.index.n_branches:
            raise ValueError(
                f"Expected {self.index.n_branches} branch bins, got {len(self.bins)}."
            )
        if len(self.alphas) != len(self.nus):
            raise ValueError(
                f"{len(self.alphas)} alphas against {len(self.nus)} nus; a clique "
                "carries one of each."
            )

    @property
    def n_cliques(self) -> int:
        return len(self.alphas)


@dataclass(frozen=True)
class TreeFieldDistribution:
    """Every tree field of one tree, and the density of the node state it came from."""

    tree_fields: np.ndarray
    """`(n_node_states, n_species_leaves, n_molecule_leaves)`."""

    probability: np.ndarray
    """`p(Z | theta)`, one per node state."""


@dataclass(frozen=True)
class Enumeration:
    """Every `(Z_s, Z_m, Y)` of one model, with the joint density of each."""

    fields: np.ndarray
    """`(n_fields, n_species_leaves, n_molecule_leaves)`, every field."""

    species_probability: np.ndarray
    """`p(Z_s | theta_s)`, one per species node state."""

    molecule_probability: np.ndarray
    """`p(Z_m | theta_m)`, one per molecule node state."""

    table: np.ndarray
    """The joint, indexed `[species node state, molecule node state, field]`."""

    @property
    def mass(self) -> float:
        """The total the model claims is one."""
        return float(self.table.sum())

    @property
    def field_rate(self) -> np.ndarray:
        """`P(Y = 1)` at every leaf pair, marginalising both node states."""
        return np.tensordot(self.table.sum(axis=(0, 1)), self.fields, axes=(0, 0))


def check_budget(n_states: int, subject: str) -> int:
    """Refuse an enumeration that does not fit. `subject` opens the message."""
    if n_states > MAX_ENUMERABLE_STATES:
        raise ValueError(
            f"{subject} needs {n_states} configurations, past the "
            f"{MAX_ENUMERABLE_STATES} this enumeration lists. The sum is "
            "intractable beyond a few leaves. Shrink the trees rather than raise "
            "the bound."
        )
    return n_states


def enumerate_binary_arrays(shape: tuple[int, ...]) -> np.ndarray:
    """Every binary array of `shape`, as `(2 ** size, *shape)`.

    Bit `k` of the row index is element `k` of the flattened array, so a field's
    cell `(i, j)` is bit `i * n_molecule_leaves + j` -- the C++'s row-major linear
    index, with the last dimension varying fastest.
    """
    size = int(np.prod(shape))
    n_states = check_budget(1 << size, f"An array of shape {tuple(shape)}")
    bits = (np.arange(n_states, dtype=np.int64)[:, None] >> np.arange(size)) & 1
    return bits.reshape((n_states, *shape)).astype(bool)


def code_of(array: np.ndarray) -> int:
    """The row `enumerate_binary_arrays` gives this array. The inverse of it."""
    flat = np.asarray(array, dtype=np.int64).reshape(-1)
    return int((flat * (1 << np.arange(len(flat)))).sum())


def node_state_distribution(
    process: TreeProcess, n_bins: int = N_BINS
) -> tuple[np.ndarray, np.ndarray]:
    """Every node state of one tree, and `p(Z | theta)` for each.

    The density is the one the draw follows. Each root takes the stationary
    probability of its own state. Every other node takes the transition
    probability given its parent. Cliques multiply, because a tree runs one
    independent process per clique.
    """
    states = enumerate_binary_arrays((process.index.n_nodes, process.n_cliques))
    matrices = transition_matrices(
        process.alphas, process.nus, grid_branch_lengths(n_bins)
    )
    cliques = np.arange(process.n_cliques)

    probability = np.ones(len(states))
    for node in range(process.index.n_nodes):
        child = states[:, node].astype(np.int64)
        parent = int(process.index.parent[node])
        if parent < 0:
            factor = np.where(child == 1, process.alphas, 1.0 - process.alphas)
        else:
            bin_ix = int(process.bins[process.index.branch_of_node[node]])
            factor = matrices[
                cliques, bin_ix, states[:, parent].astype(np.int64), child
            ]
        probability *= factor.prod(axis=1)

    return states, probability


def species_tree_fields(
    process: TreeProcess, n_bins: int = N_BINS
) -> TreeFieldDistribution:
    """The species tree's leaf block, which is already `[species, molecule]`."""
    states, probability = node_state_distribution(process, n_bins)
    return TreeFieldDistribution(states[:, process.index.leaves, :], probability)


def molecule_tree_fields(
    process: TreeProcess, n_bins: int = N_BINS
) -> TreeFieldDistribution:
    """The molecule tree's leaf block, transposed into `[species, molecule]`.

    The molecule tree draws its own nodes down the rows. A tree field and the
    field share their subscripts (ADR-0005), so its leaf block turns over.
    """
    states, probability = node_state_distribution(process, n_bins)
    return TreeFieldDistribution(
        states[:, process.index.leaves, :].transpose(0, 2, 1), probability
    )


def joint(
    species: TreeProcess,
    molecules: TreeProcess,
    omega: float,
    n_bins: int = N_BINS,
) -> Enumeration:
    """Enumerate `p(Z_s) p(Z_m) p(Y | Z_s, Z_m, omega)` over every configuration."""
    if species.n_cliques != molecules.index.n_leaves:
        raise ValueError(
            f"The species tree has {species.n_cliques} cliques against "
            f"{molecules.index.n_leaves} molecule leaves; a species clique is named "
            "by a molecule leaf."
        )
    if molecules.n_cliques != species.index.n_leaves:
        raise ValueError(
            f"The molecule tree has {molecules.n_cliques} cliques against "
            f"{species.index.n_leaves} species leaves; a molecule clique is named "
            "by a species leaf."
        )

    return joint_over_tree_fields(
        species_tree_fields(species, n_bins),
        molecule_tree_fields(molecules, n_bins),
        omega,
    )


def joint_over_tree_fields(
    species: TreeFieldDistribution,
    molecules: TreeFieldDistribution,
    omega: float,
) -> Enumeration:
    """The joint over every `(Z_s, Z_m, Y)`, from the two tree fields.

    Both densities enter as they are given, so a caller can deform one and watch
    the total leave one. That is what makes the assertion a test. A sum that only
    ever returns one proves nothing.
    """
    if species.tree_fields.shape[1:] != molecules.tree_fields.shape[1:]:
        raise ValueError(
            f"The two tree fields are {species.tree_fields.shape[1:]} and "
            f"{molecules.tree_fields.shape[1:]}; both are indexed "
            "[species_leaf, molecule_leaf] and must agree."
        )

    fields = enumerate_binary_arrays(species.tree_fields.shape[1:])
    n_species = len(species.probability)
    n_molecules = len(molecules.probability)
    n_fields = len(fields)
    check_budget(n_species * n_molecules * n_fields, "The joint")

    flat_fields = fields.reshape(n_fields, -1)
    table = np.empty((n_species, n_molecules, n_fields))
    for ix in range(n_species):
        # One species node state against every molecule node state at once.
        p_one = L.prob_field_is_one(
            species.tree_fields[ix], molecules.tree_fields, omega
        ).reshape(n_molecules, 1, -1)
        cells = np.where(flat_fields[None, :, :], p_one, 1.0 - p_one)
        table[ix] = (
            species.probability[ix]
            * molecules.probability[:, None]
            * cells.prod(axis=2)
        )

    return Enumeration(
        fields=fields,
        species_probability=species.probability,
        molecule_probability=molecules.probability,
        table=table,
    )
