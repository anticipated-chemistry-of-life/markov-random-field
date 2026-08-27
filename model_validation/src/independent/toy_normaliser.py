"""A field small enough to enumerate exactly, to isolate the normalising constant.

The field's joint density is the product of the two trees' likelihoods
(`TMarkovField::_calculate_complete_joint_density`, and the same product drives
every acceptance ratio):

    U(Y, Z_s, Z_m | theta) = p_s(Y, Z_s | theta_s) * p_m(Y, Z_m | theta_m)

Each factor on its own is a normalised tree likelihood. Their product is not,
because the two share Y:

    sum over everything of U = sum_Y p_s(Y | theta_s) * p_m(Y | theta_m) = C(theta)

and `C` depends on both parameter sets. A sampler that treats `U` as a density
targets `p_correct(theta) * C(theta)` rather than `p_correct(theta)`, so it is
pulled toward whatever maximises `C`.

This module computes `C` exactly for two balanced trees of equal depth. Depth 2
gives seven nodes each and a 4x4 field with 2^16 states — the doubly balanced
configuration in miniature, and small enough that no approximation enters
anywhere. Depth 3 would need 2^64 and is out of reach, which is the whole
difficulty: `C` is exactly the intractable object.

The neutral case that `ADR-0001` relies on falls out as a special case: when the
molecules rows are (0.5, 0.5), `p_m(Y)` is the same for every `Y`, so `C` factors
into a constant and the bias vanishes. That is *why* the independent-field
harness passes, and why it cannot see this.
"""

from __future__ import annotations

import functools
import itertools

import numpy as np

from .field import N_BINS, transition_matrix

MAX_ENUMERABLE_DEPTH = 2


@functools.lru_cache(maxsize=None)
def balanced_tree(depth: int) -> tuple[np.ndarray, tuple[int, ...]]:
    """A balanced binary tree: `(parent, leaves)`, node 0 the root.

    Node `i` has children `2i + 1` and `2i + 2`, so the leaves are the last
    `2 ** depth` nodes and the ordering is topological.
    """
    n_nodes = 2 ** (depth + 1) - 1
    parent = np.array([-1] + [(i - 1) // 2 for i in range(1, n_nodes)])
    leaves = tuple(range(2**depth - 1, n_nodes))
    return parent, leaves


def leaf_pattern_probabilities(
    alpha: float, nu: float, depth: int, bin_index: int = 5, n_bins: int = N_BINS
) -> np.ndarray:
    """`p(leaf pattern)` for one clique, marginalising the internal nodes.

    Returns `2 ** n_leaves` probabilities summing to 1, indexed by the pattern
    code: bit `k` of the code is the state of the `k`-th leaf.

    Every branch sits in the same bin, whose grid length at the default bin 5 and
    ten bins is exactly 1.0.
    """
    parent, leaves = balanced_tree(depth)
    internals = [node for node in range(len(parent)) if node not in set(leaves)]

    matrix = transition_matrix(alpha, nu, (bin_index + 0.5) * 2.0 / (n_bins + 1.0))
    stationary = (1.0 - alpha, alpha)

    out = np.zeros(2 ** len(leaves))
    for assignment in itertools.product((0, 1), repeat=len(internals)):
        state = dict(zip(internals, assignment))
        weight = stationary[state[0]]
        for node in internals[1:]:
            weight *= matrix[state[parent[node]], state[node]]

        # Each leaf contributes independently given its parent, so the pattern
        # distribution for this internal assignment is an outer product rather
        # than another 2^n_leaves loop.
        block = np.array([weight])
        for leaf in leaves:
            row = matrix[state[parent[leaf]]]
            block = np.concatenate([block * row[0], block * row[1]])
        out += block

    return out


@functools.lru_cache(maxsize=None)
def pattern_codes(depth: int) -> tuple[np.ndarray, np.ndarray]:
    """Row and column pattern codes for every field of two depth-`depth` trees.

    Cell `(i, j)` of the field is bit `i * n_leaves + j`, matching the C++'s
    row-major linear index with the last dimension varying fastest.
    """
    if depth > MAX_ENUMERABLE_DEPTH:
        raise ValueError(
            f"Depth {depth} needs 2^{4 ** depth} fields. Enumeration stops at "
            f"depth {MAX_ENUMERABLE_DEPTH}; C is intractable beyond it, which is "
            "the point."
        )

    n_leaves = 2**depth
    n_cells = n_leaves * n_leaves
    fields = np.arange(2**n_cells, dtype=np.int64)
    bits = ((fields[:, None] >> np.arange(n_cells)[None, :]) & 1)
    bits = bits.reshape(-1, n_leaves, n_leaves)

    powers = 2 ** np.arange(n_leaves)
    rows = (bits * powers[None, None, :]).sum(axis=2)
    cols = (bits * powers[None, :, None]).sum(axis=1)
    return rows, cols


def field_factors(
    species: np.ndarray, molecules: np.ndarray, depth: int
) -> tuple[np.ndarray, np.ndarray]:
    """`p_s(Y)` and `p_m(Y)` for every field, from the two pattern distributions.

    A species clique is one column of the field and a molecules clique one row,
    so each factor is a product over the cliques of that tree.
    """
    rows, cols = pattern_codes(depth)
    return np.prod(species[cols], axis=1), np.prod(molecules[rows], axis=1)


def normalising_constant(
    species: np.ndarray, molecules: np.ndarray, depth: int
) -> float:
    """`C(theta) = sum_Y p_s(Y) * p_m(Y)`, the constant the sampler omits."""
    p_species, p_molecules = field_factors(species, molecules, depth)
    return float(np.dot(p_species, p_molecules))


def field_distribution(
    species: np.ndarray, molecules: np.ndarray, depth: int
) -> np.ndarray:
    """The properly normalised field distribution, `p_s(Y) p_m(Y) / C`."""
    p_species, p_molecules = field_factors(species, molecules, depth)
    joint = p_species * p_molecules
    return joint / joint.sum()


def expected_log_likelihood_profile(
    truth: np.ndarray,
    log_nu_grid: np.ndarray,
    alpha_species: float,
    molecules: np.ndarray,
    depth: int,
    bin_index: int = 5,
    n_bins: int = N_BINS,
) -> tuple[np.ndarray, np.ndarray]:
    """Profile both objectives in the species `log nu`, under the true field law.

    Returns `(targeted, correct)`, each the expectation under `truth` of one
    per-field objective:

    - `targeted` is what the C++ maximises: `log p_s(Y)`, with `log p_m(Y)`
      dropped as a constant in `theta_s`.
    - `correct` subtracts `log C(theta_s, theta_m)`, which does depend on
      `theta_s` whenever the molecules dimension is not neutral.

    Using the expectation rather than a sample removes Monte Carlo noise: by
    Gibbs' inequality `correct` is maximised at the truth, so any gap between the
    two argmaxes is the bias, not sampling error.
    """
    targeted = np.empty(len(log_nu_grid))
    correct = np.empty(len(log_nu_grid))

    for ix, log_nu in enumerate(log_nu_grid):
        species = leaf_pattern_probabilities(
            alpha_species, float(np.exp(log_nu)), depth, bin_index, n_bins
        )
        p_species, p_molecules = field_factors(species, molecules, depth)

        log_p_species = np.log(np.maximum(p_species, 1e-300))
        targeted[ix] = float(np.dot(truth, log_p_species))
        correct[ix] = targeted[ix] - np.log(np.dot(p_species, p_molecules))

    return targeted, correct


def run_chain(
    rng: np.random.Generator,
    *,
    correct_normaliser: bool,
    true_log_nu_species: float,
    log_nu_molecules: float,
    alpha_species: float,
    alpha_molecules: float,
    depth: int,
    n_iterations: int,
    n_replicates: int = 200,
    proposal_sd: float = 0.15,
    prior_sd: float = 10.0,
    bin_index: int = 5,
    n_bins: int = N_BINS,
) -> np.ndarray:
    """Metropolis on the species `log nu`, over fields observed without error.

    `n_replicates` independent fields are drawn once from the correctly
    normalised law at the true parameters and then held fixed, which is the
    best case for inference: the field is known exactly, so nothing here depends
    on the data model, on Gibbs mixing, or on how well `Y` is recovered. Many
    replicates rather than one only sharpen the posterior; they do not change
    where it sits.

    With `correct_normaliser` false the acceptance ratio is the C++'s — the
    species likelihood ratio plus the prior ratio. With it true the ratio also
    carries `-log C(theta)`. Nothing else differs between the two, so a gap
    between the traces is attributable to that term alone.

    The prior on `log nu` is a Normal centred at zero and wide enough that it
    cannot itself account for a drift.
    """
    molecules = leaf_pattern_probabilities(
        alpha_molecules, float(np.exp(log_nu_molecules)), depth, bin_index, n_bins
    )
    truth_species = leaf_pattern_probabilities(
        alpha_species, float(np.exp(true_log_nu_species)), depth, bin_index, n_bins
    )
    truth = field_distribution(truth_species, molecules, depth)
    fields = rng.choice(len(truth), size=n_replicates, p=truth)

    def log_likelihood(log_nu: float) -> float:
        species = leaf_pattern_probabilities(
            alpha_species, float(np.exp(log_nu)), depth, bin_index, n_bins
        )
        p_species, p_molecules = field_factors(species, molecules, depth)
        total = float(np.sum(np.log(np.maximum(p_species[fields], 1e-300))))
        if correct_normaliser:
            total -= n_replicates * float(np.log(np.dot(p_species, p_molecules)))
        return total

    def log_prior(log_nu: float) -> float:
        return -0.5 * (log_nu / prior_sd) ** 2

    log_nu = true_log_nu_species
    current = log_likelihood(log_nu) + log_prior(log_nu)
    trace = np.empty(n_iterations)

    for iteration in range(n_iterations):
        proposed = log_nu + rng.normal(0.0, proposal_sd)
        candidate = log_likelihood(proposed) + log_prior(proposed)
        if np.log(rng.random()) < candidate - current:
            log_nu, current = proposed, candidate
        trace[iteration] = log_nu

    return trace
