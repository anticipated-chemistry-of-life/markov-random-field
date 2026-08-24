"""The two-state tree process, derived from the maths rather than from the C++.

The process on a branch of length `t` is a two-state continuous-time Markov chain
with stationary distribution `(1 - alpha, alpha)` and switching rate `nu`:

    Lambda = [[  -alpha * nu   ,     alpha * nu   ],
              [ (1-alpha) * nu , (alpha - 1) * nu ]]

whose matrix exponential has the closed form implemented in `transition_matrix`.

The C++ never evaluates that exponential at `t` directly: it builds bin 0 as
`expm(Lambda * Delta / 2)` and then walks up the grid by repeated multiplication,
`P_k = P_{k-1} @ expm(Lambda * Delta)` (`src/TClique.h:103`). This module
evaluates the closed form at `t = Delta * (k + 0.5)` instead. That divergence is
deliberate: replicating the recursion would reproduce any error accumulating in
it rather than exposing it.
"""

from __future__ import annotations

import numpy as np

# ProgramOptions::BRANCH_LENGTHS_BINS, src/cli.h:28
N_BINS = 10

# Above this rate the C++ short-circuits to the stationary distribution
# (TMatrices::set_lambda, src/TClique.h:161). This module never replicates that
# approximation; callers assert which side of it they are on.
STATIONARY_NU_THRESHOLD = 25.0


def delta(n_bins: int = N_BINS) -> float:
    """Bin width. `_initialize_grid_branch_lengths`, tree_branches.cpp:17."""
    return 2.0 / (n_bins + 1.0)


def grid_branch_lengths(n_bins: int = N_BINS) -> np.ndarray:
    """The continuous length each bin stands for: its centre, `Delta * (k + 0.5)`."""
    return delta(n_bins) * (np.arange(n_bins) + 0.5)


def transition_matrix(alpha: float, nu: float, t: float) -> np.ndarray:
    """`P[i, j] = P(child = j | parent = i)` after time `t`."""
    return transition_matrices(np.array([alpha]), np.array([nu]), np.array([t]))[0, 0]


def transition_matrices(
    alphas: np.ndarray, nus: np.ndarray, times: np.ndarray
) -> np.ndarray:
    """Transition matrices for every (clique, time) pair.

    Returns shape `(n_cliques, n_times, 2, 2)`, indexed `[c, k, parent, child]`.
    """
    alpha = np.asarray(alphas, dtype=float)[:, None]
    decay = np.exp(-np.asarray(nus, dtype=float)[:, None] * np.asarray(times)[None, :])

    out = np.empty(decay.shape + (2, 2), dtype=float)
    out[..., 0, 0] = 1.0 - alpha + alpha * decay
    out[..., 0, 1] = alpha * (1.0 - decay)
    out[..., 1, 0] = (1.0 - alpha) * (1.0 - decay)
    out[..., 1, 1] = alpha + (1.0 - alpha) * decay
    return out


def branch_length_budget(n_branches: int, n_bins: int = N_BINS) -> int:
    """The total of a tree's bins, which the MCMC conserves for life.

    `_bin_branch_lengths` forces this total at startup (tree_branches.cpp:191),
    and every proposal thereafter moves +1 on one branch and -1 on another
    (tree_branches.cpp:46-47). Bins summing to anything else are unreachable.
    """
    return n_branches * n_bins // 2


def repair_to_budget(
    rng: np.random.Generator, bins: np.ndarray, n_bins: int = N_BINS
) -> np.ndarray:
    """Walk `bins` onto the branch-length budget by random +-1 steps."""
    bins = np.array(bins, dtype=np.int64)
    goal = branch_length_budget(len(bins), n_bins)
    total = int(bins.sum())

    while total != goal:
        ix = int(rng.integers(0, len(bins)))
        if total < goal:
            if bins[ix] >= n_bins - 1:
                continue
            step = 1
        else:
            if bins[ix] <= 0:
                continue
            step = -1
        bins[ix] += step
        total += step

    return bins


def sample_binned_branch_lengths(
    rng: np.random.Generator, n_branches: int, n_bins: int = N_BINS
) -> np.ndarray:
    """Draw bins uniformly on `0 .. n_bins-1`, then repair onto the budget."""
    return repair_to_budget(rng, rng.integers(0, n_bins, size=n_branches), n_bins)


def bin_from_length(length: float, n_bins: int = N_BINS) -> int:
    """Replica of `TTree::_get_bin_branch_length` (tree_branches.cpp:28)."""
    if length <= 0.0:
        return 0
    return min(int(length / delta(n_bins)), n_bins - 1)


def bins_from_tree_lengths(lengths: np.ndarray, n_bins: int = N_BINS) -> np.ndarray:
    """Replica of the C++ read path: normalise to mean 1, then bin.

    Mirrors `_bin_branch_lengths_from_tree` (tree_branches.cpp:212) *without* the
    trailing repair, so tests can check that writing grid centres round-trips
    back to the bins that produced them.
    """
    lengths = np.asarray(lengths, dtype=float)
    positive = lengths[lengths > 0.0]
    normalised = lengths / positive.mean()
    return np.array([bin_from_length(x, n_bins) for x in normalised], dtype=np.int64)


def sample_states(
    rng: np.random.Generator,
    tree,
    bins: np.ndarray,
    alphas: np.ndarray,
    nus: np.ndarray,
    n_bins: int = N_BINS,
) -> np.ndarray:
    """One pass down the tree, for every clique at once.

    Roots are drawn from the stationary distribution `(1 - alpha, alpha)`; every
    other node is drawn given its parent. Returns shape `(n_nodes, n_cliques)`.
    """
    alphas = np.asarray(alphas, dtype=float)
    nus = np.asarray(nus, dtype=float)
    n_cliques = len(alphas)
    if len(nus) != n_cliques:
        raise ValueError("alphas and nus must have one entry per clique.")
    if len(bins) != tree.n_branches:
        raise ValueError(
            f"Expected {tree.n_branches} branch lengths, got {len(bins)}."
        )

    # p_one[c, k, parent] = P(child = 1 | parent, clique c, bin k)
    matrices = transition_matrices(alphas, nus, grid_branch_lengths(n_bins))
    p_one = matrices[..., 1]

    states = np.empty((tree.n_nodes, n_cliques), dtype=bool)
    cliques = np.arange(n_cliques)

    # Node order is topological, so a single forward pass suffices.
    for node in range(tree.n_nodes):
        parent = tree.parent[node]
        if parent < 0:
            probability = alphas
        else:
            bin_ix = bins[tree.branch_of_node[node]]
            probability = p_one[cliques, bin_ix, states[parent].astype(np.int64)]
        states[node] = rng.random(n_cliques) < probability

    return states


def sibling_disagreement_probability(
    alpha: float, nu: float, bin_left: int, bin_right: int, n_bins: int = N_BINS
) -> float:
    """`P(left leaf != right leaf)` for two siblings, marginalising the parent.

    The parent of a leaf is an internal node, whose marginal is stationary, so
    this is a closed-form prediction that both implementations must reproduce.
    """
    grid = grid_branch_lengths(n_bins)
    left = transition_matrix(alpha, nu, grid[bin_left])
    right = transition_matrix(alpha, nu, grid[bin_right])
    stationary = np.array([1.0 - alpha, alpha])
    return float(
        sum(
            stationary[s] * (left[s, 0] * right[s, 1] + left[s, 1] * right[s, 0])
            for s in (0, 1)
        )
    )
