"""Verification of the independent reference implementation.

Without these, an unverified Python would be judging an unverified C++, and a
disagreement between them would say nothing about which side is wrong.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.linalg import expm

from src.independent import field as F
from src.independent import scenario
from src.independent.data import research_effort, simulate_simple_error
from src.independent.indexing import build_tree_index
from src.tree import Tree, TreeType

ALPHAS = [0.05, 0.3, 0.5, 0.72, 0.95]
NUS = [0.05, 0.6, 3.0]
TIMES = [0.05, 0.5, 1.0, 1.7]


def edges_of(tree: Tree) -> list[tuple[str, str]]:
    frame = tree.to_dataframe()
    return list(zip(frame["child"].astype(str), frame["parent"].astype(str)))


# --------------------------------------------------------------------------
# The transition matrix
# --------------------------------------------------------------------------


@pytest.mark.parametrize("alpha", ALPHAS)
@pytest.mark.parametrize("nu", NUS)
@pytest.mark.parametrize("t", TIMES)
def test_rows_are_distributions(alpha, nu, t):
    matrix = F.transition_matrix(alpha, nu, t)
    assert np.allclose(matrix.sum(axis=1), 1.0)
    assert (matrix >= 0.0).all()


@pytest.mark.parametrize("alpha", ALPHAS)
@pytest.mark.parametrize("nu", NUS)
def test_zero_time_is_identity(alpha, nu):
    assert np.allclose(F.transition_matrix(alpha, nu, 0.0), np.eye(2))


@pytest.mark.parametrize("alpha", ALPHAS)
@pytest.mark.parametrize("nu", NUS)
def test_infinite_time_is_stationary(alpha, nu):
    matrix = F.transition_matrix(alpha, nu, 1e6)
    assert np.allclose(matrix, np.array([[1 - alpha, alpha], [1 - alpha, alpha]]))


@pytest.mark.parametrize("alpha", ALPHAS)
@pytest.mark.parametrize("nu", NUS)
@pytest.mark.parametrize("t", TIMES)
def test_stationary_distribution_is_preserved(alpha, nu, t):
    stationary = np.array([1 - alpha, alpha])
    assert np.allclose(stationary @ F.transition_matrix(alpha, nu, t), stationary)


@pytest.mark.parametrize("alpha", ALPHAS)
@pytest.mark.parametrize("nu", NUS)
@pytest.mark.parametrize("t", TIMES)
def test_closed_form_matches_matrix_exponential(alpha, nu, t):
    """Two independent routes to the same number: analytic and numerical."""
    generator = np.array(
        [[-alpha * nu, alpha * nu], [(1 - alpha) * nu, (alpha - 1) * nu]]
    )
    assert np.allclose(F.transition_matrix(alpha, nu, t), expm(generator * t))


def test_generator_rows_sum_to_zero():
    alpha, nu = 0.3, 0.7
    generator = np.array(
        [[-alpha * nu, alpha * nu], [(1 - alpha) * nu, (alpha - 1) * nu]]
    )
    assert np.allclose(generator.sum(axis=1), 0.0)


# --------------------------------------------------------------------------
# Branch length bins
# --------------------------------------------------------------------------


def test_grid_centres_have_mean_one_at_budget():
    """The budget is exactly the constraint 'mean grid branch length is 1'."""
    n_branches, n_bins = 254, F.N_BINS
    budget = F.branch_length_budget(n_branches, n_bins)
    mean_bin = budget / n_branches
    assert np.isclose(F.delta(n_bins) * (mean_bin + 0.5), 1.0)


@pytest.mark.parametrize("n_branches", [2, 7, 64, 254])
def test_repair_reaches_the_budget(n_branches):
    rng = np.random.default_rng(0)
    bins = F.repair_to_budget(rng, rng.integers(0, F.N_BINS, size=n_branches))
    assert bins.sum() == F.branch_length_budget(n_branches)
    assert bins.min() >= 0
    assert bins.max() <= F.N_BINS - 1


def test_repair_is_a_noop_when_already_on_budget():
    rng = np.random.default_rng(1)
    bins = np.full(100, F.N_BINS // 2, dtype=np.int64)
    assert np.array_equal(F.repair_to_budget(rng, bins), bins)


@pytest.mark.parametrize("seed", range(5))
def test_sampled_bins_are_on_budget(seed):
    rng = np.random.default_rng(seed)
    bins = F.sample_binned_branch_lengths(rng, 254)
    assert bins.sum() == F.branch_length_budget(254)


@pytest.mark.parametrize("k", range(F.N_BINS))
def test_bin_of_its_own_grid_centre(k):
    assert F.bin_from_length(F.grid_branch_lengths()[k]) == k


@pytest.mark.parametrize("seed", range(5))
def test_grid_centres_round_trip_through_the_cpp_read_path(seed):
    """Writing grid centres recovers the bins that produced them.

    Only true because bins on budget have mean grid length exactly 1, so the
    C++'s normalise-to-mean-1 step is a no-op.
    """
    rng = np.random.default_rng(seed)
    bins = F.sample_binned_branch_lengths(rng, 254)
    lengths = F.grid_branch_lengths()[bins]
    assert np.array_equal(F.bins_from_tree_lengths(lengths), bins)


def test_flat_lengths_bin_to_the_middle():
    """A flat tree file starts every branch at the budget's mean bin."""
    bins = F.bins_from_tree_lengths(np.full(254, 0.2))
    assert set(np.unique(bins)) == {F.N_BINS // 2}
    assert bins.sum() == F.branch_length_budget(254)


# --------------------------------------------------------------------------
# Node ordering
# --------------------------------------------------------------------------


@pytest.mark.parametrize("n_nodes", [7, 15, 255])
def test_ordering_matches_the_existing_replica(n_nodes):
    tree = Tree(n_nodes, TreeType.balanced, "species")
    index = build_tree_index(edges_of(tree))
    assert index.leaf_names() == tree.get_leaf_names_in_cpp_order()
    assert index.branch_names() == tree.get_non_root_node_names_in_cpp_order()


def test_internal_nodes_include_the_root():
    tree = Tree(255, TreeType.balanced, "species")
    index = build_tree_index(edges_of(tree))
    assert index.n_leaves == 128
    assert index.n_internals == 127
    assert index.n_branches == 254
    roots = [n for n in range(index.n_nodes) if index.parent[n] < 0]
    assert len(roots) == 1
    assert roots[0] in set(index.internals.tolist())


def test_every_parent_precedes_its_child():
    index = build_tree_index(edges_of(Tree(255, TreeType.balanced, "species")))
    for node in range(index.n_nodes):
        assert index.parent[node] < node


# --------------------------------------------------------------------------
# The sampler
# --------------------------------------------------------------------------


@pytest.mark.parametrize("alpha", [0.2, 0.5, 0.8])
def test_leaf_marginal_converges_to_alpha(alpha):
    """Every node's marginal is stationary, so leaves average to alpha."""
    tree = Tree(15, TreeType.balanced, "species")
    index = build_tree_index(edges_of(tree))
    rng = np.random.default_rng(20)

    n_cliques = 4000
    states = F.sample_states(
        rng,
        index,
        F.sample_binned_branch_lengths(rng, index.n_branches),
        np.full(n_cliques, alpha),
        np.full(n_cliques, 0.6),
    )
    assert states[index.leaves].mean() == pytest.approx(alpha, abs=0.03)


def test_a_frozen_process_copies_the_root():
    """As nu goes to zero, no branch ever switches."""
    index = build_tree_index(edges_of(Tree(15, TreeType.balanced, "species")))
    rng = np.random.default_rng(3)
    states = F.sample_states(
        rng,
        index,
        F.sample_binned_branch_lengths(rng, index.n_branches),
        np.full(200, 0.5),
        np.full(200, 1e-9),
    )
    assert (states == states[0]).all()


def test_a_fast_process_forgets_the_parent():
    """Past the stationary threshold, children are independent of their parent."""
    index = build_tree_index(edges_of(Tree(15, TreeType.balanced, "species")))
    rng = np.random.default_rng(4)
    alpha = 0.5
    states = F.sample_states(
        rng,
        index,
        F.sample_binned_branch_lengths(rng, index.n_branches),
        np.full(4000, alpha),
        np.full(4000, 1e4),
    )
    children = [n for n in range(index.n_nodes) if index.parent[n] >= 0]
    agreement = np.mean([states[n] == states[index.parent[n]] for n in children])
    assert agreement == pytest.approx(0.5, abs=0.02)


def test_sampler_rejects_a_branch_length_mismatch():
    index = build_tree_index(edges_of(Tree(15, TreeType.balanced, "species")))
    rng = np.random.default_rng(5)
    with pytest.raises(ValueError, match="branch lengths"):
        F.sample_states(rng, index, np.zeros(3, dtype=np.int64), [0.5], [0.6])


@pytest.mark.parametrize("bin_left", [0, 4, 9])
@pytest.mark.parametrize("bin_right", [0, 7])
def test_sibling_disagreement_matches_simulation(bin_left, bin_right):
    """The Z-free statistic's analytic prediction is checked against sampling."""
    alpha, nu = 0.35, 0.8
    index = build_tree_index([("leaf_l", "root"), ("leaf_r", "root")])
    bins = np.empty(2, dtype=np.int64)
    bins[index.branch_of_node[index.names.index("leaf_l")]] = bin_left
    bins[index.branch_of_node[index.names.index("leaf_r")]] = bin_right

    n_cliques = 200_000
    states = F.sample_states(
        np.random.default_rng(7),
        index,
        bins,
        np.full(n_cliques, alpha),
        np.full(n_cliques, nu),
    )
    observed = np.mean(
        states[index.names.index("leaf_l")] != states[index.names.index("leaf_r")]
    )
    predicted = F.sibling_disagreement_probability(alpha, nu, bin_left, bin_right)
    assert observed == pytest.approx(predicted, abs=0.005)


# --------------------------------------------------------------------------
# Neutrality of a pinned dimension
# --------------------------------------------------------------------------


def test_neutral_parameters_give_exactly_uninformative_rows():
    """The assumption the whole validation rests on (ADR-0001)."""
    nu = np.exp(5.0)
    assert nu > F.STATIONARY_NU_THRESHOLD
    matrix = F.transition_matrix(0.5, nu, F.grid_branch_lengths()[0])
    assert np.allclose(matrix, 0.5)


def test_neutral_rows_are_identical_across_every_bin():
    nu = np.exp(5.0)
    matrices = F.transition_matrices(
        np.array([0.5]), np.array([nu]), F.grid_branch_lengths()
    )
    assert np.allclose(matrices, 0.5)


# --------------------------------------------------------------------------
# Observation models
# --------------------------------------------------------------------------


def test_initial_value_filenames_keep_their_reader_marker():
    """stattools picks a reader by filename, not by content.

    A name/value file is only matched up by parameter name when its filename
    contains one of these markers (TReadInitialValues.h:133); otherwise the C++
    rejects it. Renaming these files without keeping a marker breaks every run
    script, so fail here rather than in a C++ stack trace.
    """
    for filename in (scenario.PINNED_MOLECULES, scenario.SIMULATE_PARAMETERS):
        assert any(m in filename for m in scenario.INITIAL_VALUE_MARKERS), filename


def test_neutral_pinning_covers_every_molecules_parameter():
    """Leaving one free would let a meaningless chain wander into the traces."""
    pinned = {flag.lstrip("-").split(".")[0]
              for flag in scenario._PIN_MOLECULES.split()
              if flag.startswith("--molecules_")}
    assert pinned == {
        "molecules_alpha", "molecules_log_nu", "molecules_mean_log_nu",
        "molecules_var_log_nu", "molecules_branch_lengths",
    }


def test_research_effort_uses_log_paper_counts():
    """The C++ log-transforms paper counts when reading them (TTree.h:517)."""
    papers = np.array([0, 1, 3, 9])
    gamma = 1.1
    expected = 1.0 - np.exp(-gamma * np.log(papers + 1.0))
    assert np.allclose(research_effort(papers, np.array([3]), gamma)[:, 0],
                       expected * (1.0 - np.exp(-gamma * np.log(4.0))))


def test_a_leaf_with_no_papers_is_never_reported():
    """log(0 + 1) = 0, so an unstudied leaf has zero research effort."""
    assert research_effort(np.array([0]), np.array([5]), gamma=1.1)[0, 0] == 0.0


def test_research_effort_is_a_probability():
    effort = research_effort(np.arange(1, 9), np.arange(1, 6), gamma=1.1)
    assert effort.shape == (8, 5)
    assert ((effort > 0.0) & (effort < 1.0)).all()


def test_research_effort_rises_with_papers():
    effort = research_effort(np.arange(1, 9), np.array([4]), gamma=1.1)[:, 0]
    assert (np.diff(effort) > 0).all()


@pytest.mark.parametrize("epsilon", [0.0, 0.05, 0.5, 1.0])
def test_simple_error_flips_at_the_stated_rate(epsilon):
    rng = np.random.default_rng(11)
    truth = rng.random((400, 400)) < 0.3
    observed = simulate_simple_error(rng, truth, epsilon)
    assert np.mean(observed != truth) == pytest.approx(epsilon, abs=0.01)
