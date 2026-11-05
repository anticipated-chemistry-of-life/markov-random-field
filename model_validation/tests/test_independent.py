"""Verification of the independent reference implementation.

Without these, an unverified Python would be judging an unverified C++, and a
disagreement between them would say nothing about which side is wrong.
"""

from __future__ import annotations

import pathlib
import re
import tempfile

import numpy as np
import pandas as pd
import pytest
from scipy.linalg import expm

from src.independent import enumeration as E
from src.independent import field as F
from src.independent import io
from src.independent import link as L
from src.independent import mass_spec
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
def test_the_three_blocks_partition_the_nodes(n_nodes):
    """Leaves, then internal non-root nodes, then roots -- and nothing else."""
    index = build_tree_index(edges_of(Tree(n_nodes, TreeType.balanced, "species")))
    n = index.n_nodes
    n_roots = int(np.sum(index.parent < 0))
    has_child = set(int(p) for p in index.parent if p >= 0)

    for node in range(n):
        is_leaf = node not in has_child
        assert is_leaf == (node < index.n_leaves), f"leaf block at {node}"
        assert (index.parent[node] < 0) == (node >= n - n_roots), f"root at {node}"
    for node in range(index.n_leaves, n - n_roots):
        assert index.parent[node] >= 0 and node in has_child


@pytest.mark.parametrize("n_nodes", [7, 15, 255])
def test_index_spaces_are_arithmetic(n_nodes):
    index = build_tree_index(edges_of(Tree(n_nodes, TreeType.balanced, "species")))
    n_roots = int(np.sum(index.parent < 0))
    assert index.leaves.tolist() == list(range(index.n_leaves))
    assert index.internals.tolist() == list(range(index.n_leaves, index.n_nodes))
    assert index.branch_nodes.tolist() == list(range(index.n_nodes - n_roots))
    n_branches = index.n_nodes - n_roots
    assert index.branch_of_node.tolist() == list(range(n_branches)) + [-1] * n_roots


def test_leaves_are_numbered_depth_first_with_siblings_adjacent():
    """The leaf order the C++ uses (ADR-0014), on a tree where it differs from file
    order: the root's leaf children `a` and `c` come together, before the clade
    below `x`. `PhylogenyOrder.leaves_are_depth_first_with_siblings_adjacent` in
    tests/TPhylogeny_Tests.cpp pins the same tree, so the two implementations
    cannot drift apart unnoticed."""
    edges = [("a", "r"), ("x", "r"), ("b", "x"), ("c", "r"), ("d", "x")]
    index = build_tree_index(edges)
    assert index.leaf_names() == ["a", "c", "b", "d"]
    assert index.internal_names() == ["x", "r"]


def test_internal_nodes_include_the_root():
    tree = Tree(255, TreeType.balanced, "species")
    index = build_tree_index(edges_of(tree))
    assert index.n_leaves == 128
    assert index.n_internals == 127
    assert index.n_branches == 254
    roots = [n for n in range(index.n_nodes) if index.parent[n] < 0]
    assert len(roots) == 1
    assert roots[0] in set(index.internals.tolist())


def test_a_tree_file_gives_each_row_its_own_child_s_branch_length():
    """Rows are in edge order, per-branch arrays are in branch order.

    The two coincided for a balanced tree while node order was file order, and
    stopped coinciding under canonical order (ADR-0004) -- so writing the lengths
    positionally would hand every branch someone else's.
    """
    tree = Tree(15, TreeType.balanced, "species")
    edges = edges_of(tree)
    index = build_tree_index(edges)
    assert [child for child, _ in edges] != index.branch_names(), (
        "this fixture no longer distinguishes edge order from branch order"
    )

    lengths = np.arange(index.n_branches, dtype=float) + 1.0
    with tempfile.TemporaryDirectory() as tmp:
        path = pathlib.Path(tmp) / "species.txt"
        io.write_tree(path, edges, index, lengths)
        frame = pd.read_csv(path, sep="\t")

    expected = dict(zip(index.branch_names(), lengths))
    for child, length in zip(frame["child"].astype(str), frame["length"]):
        assert length == expected[child], f"{child} got the wrong branch length"


def test_the_true_branch_lengths_option_writes_the_true_branch_lengths():
    """--true_branch_lengths must start the chain *on* the truth, not on a
    permutation of it. The truth file is name-keyed, the tree file is not."""
    with tempfile.TemporaryDirectory() as tmp:
        out = pathlib.Path(tmp) / "scenario"
        scenario.build_scenario(
            out,
            scenario.ScenarioConfig(
                n_species_nodes=15, n_molecule_nodes=15, true_branch_lengths=True
            ),
        )
        tree_file = pd.read_csv(out / "species.txt", sep="\t")
        truth = pd.read_csv(out / "truth_species.txt", sep="\t")

    bin_of = {
        name.removeprefix("species_branch_lengths_"): int(value)
        for name, value in zip(truth["name"].astype(str), truth["value"])
        if name.startswith("species_branch_lengths_")
    }
    grid = F.grid_branch_lengths()
    for child, length in zip(tree_file["child"].astype(str), tree_file["length"]):
        assert length == pytest.approx(grid[bin_of[child]]), child


def test_every_child_precedes_its_parent():
    """The post-order guarantee, stated as the property the sampler relies on."""
    index = build_tree_index(edges_of(Tree(255, TreeType.balanced, "species")))
    for node in range(index.n_nodes):
        if index.parent[node] >= 0:
            assert index.parent[node] > node


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
#
# No scenario neutralises anything any more: the reference is exact with both
# trees active (ADR-0005). What a neutral dimension *is* is still worth pinning
# down, because ADR-0001's rung survives as a cheap regression check.
# --------------------------------------------------------------------------


def test_neutral_parameters_give_exactly_uninformative_rows():
    """What makes a dimension neutral: transition rows of exactly (0.5, 0.5)."""
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
# The link, and the field it draws
# --------------------------------------------------------------------------

OMEGAS = [1e-4, 0.01, 0.05, 0.2, 0.4999]


def _brute_force_link(z_s: bool, z_m: bool, omega: float) -> float:
    """`P(Y = 1 | Z_s, Z_m)` by summing over both corruption events.

    The closed form in `link` is a product of two independent factors. This is the
    definition it came from: corrupt each cell, then AND the results. Nothing here
    is derived from that product, so agreement is a check and not a restatement.
    """
    total = 0.0
    for corrupted_s in (False, True):
        for corrupted_m in (False, True):
            probability = (omega if corrupted_s else 1.0 - omega) * (
                omega if corrupted_m else 1.0 - omega
            )
            read_s = (not z_s) if corrupted_s else z_s
            read_m = (not z_m) if corrupted_m else z_m
            if read_s and read_m:
                total += probability
    return total


@pytest.mark.parametrize("omega", OMEGAS)
def test_the_link_is_an_and_over_two_independently_corrupted_reads(omega):
    for z_s in (False, True):
        for z_m in (False, True):
            expected = _brute_force_link(z_s, z_m, omega)
            got = float(L.prob_field_is_one(np.array([z_s]), np.array([z_m]), omega)[0])
            assert got == pytest.approx(expected), (z_s, z_m, omega)


@pytest.mark.parametrize("omega", OMEGAS)
def test_the_bucket_pools_the_two_mixed_cells(omega):
    """The table depends on the two tree fields only through their sum."""
    mixed = float(L.prob_field_is_one(np.array([True]), np.array([False]), omega)[0])
    other = float(L.prob_field_is_one(np.array([False]), np.array([True]), omega)[0])
    assert mixed == pytest.approx(other)
    assert float(L.prob_for_bucket(1, omega)) == pytest.approx(mixed)


@pytest.mark.parametrize("omega", OMEGAS)
def test_both_parameter_free_constraints_hold_at_every_error_probability(omega):
    """Three Bernoulli rates pinned by one parameter (ADR-0005, derivation 2)."""
    p_0, p_1, p_2 = L.prob_for_bucket(np.arange(3), omega)
    assert p_1**2 == pytest.approx(p_0 * p_2, abs=1e-15)
    assert np.sqrt(p_0) + np.sqrt(p_2) == pytest.approx(1.0)


@pytest.mark.parametrize("omega", [0.0, -0.1, 0.5, 0.9, 1.0])
def test_the_error_probability_must_lie_inside_the_open_interval(omega):
    """At 0 the link is deterministic; at 0.5 and above the tree fields are
    anti-correlated with the field. Both are statements about the model."""
    with pytest.raises(ValueError):
        L.check_error_probability(omega)


def test_the_field_is_drawn_at_the_rate_the_link_names():
    rng = np.random.default_rng(20260907)
    omega = 0.15
    z_s = np.array([[True, True], [False, False]])
    z_m = np.array([[True, False], [True, False]])

    draws = np.mean(
        [L.sample_field(rng, z_s, z_m, omega) for _ in range(20000)], axis=0
    )
    assert draws == pytest.approx(L.prob_field_is_one(z_s, z_m, omega), abs=0.01)


def test_the_counters_tally_every_cell_once():
    rng = np.random.default_rng(20260908)
    z_s = rng.random((17, 13)) < 0.4
    z_m = rng.random((17, 13)) < 0.6
    field = L.sample_field(rng, z_s, z_m, 0.1)

    counters = L.link_counters(z_s, z_m, field)
    assert counters.sum() == z_s.size
    for bucket in range(L.N_BUCKETS):
        holding = L.buckets(z_s, z_m) == bucket
        assert counters[bucket, 1] == int(field[holding].sum())
        assert counters[bucket, 0] == int((~field[holding]).sum())


def test_the_marginal_field_rate_is_the_product_of_the_two_adjusted_rates():
    """ADR-0005, derivation 3. The field's density says nothing about the split.

    Two very different pairs of alphas with the same product of adjusted rates
    give the same field density, which is the identifiability limit the rung
    ladder is meant to report rather than be surprised by.
    """
    rng = np.random.default_rng(20260909)
    omega = 0.1
    shape = (400, 400)

    def density(alpha_s: float, alpha_m: float) -> float:
        z_s = rng.random(shape) < alpha_s
        z_m = rng.random(shape) < alpha_m
        return float(L.sample_field(rng, z_s, z_m, omega).mean())

    for alpha_s, alpha_m in ((0.8, 0.3), (0.3, 0.8), (0.5, 0.5)):
        expected = float(
            L.adjusted_rate(alpha_s, omega) * L.adjusted_rate(alpha_m, omega)
        )
        assert density(alpha_s, alpha_m) == pytest.approx(expected, abs=0.005)


# --------------------------------------------------------------------------
# The scenario, with both trees active
# --------------------------------------------------------------------------


def _small_scenario(tmp: str, **overrides):
    out = pathlib.Path(tmp) / "scenario"
    config = scenario.ScenarioConfig(
        n_species_nodes=63, n_molecule_nodes=31, **overrides
    )
    return out, config, scenario.build_scenario(out, config)


def _indices_of(out: pathlib.Path):
    """The two tree indices, read back from the tree files the scenario wrote."""

    def index(name: str):
        frame = pd.read_csv(out / f"{name}.txt", sep="\t")
        return build_tree_index(
            list(zip(frame["child"].astype(str), frame["parent"].astype(str)))
        )

    return index("species"), index("molecules")


def test_the_scenario_draws_every_node_of_both_trees():
    """A node state spans every node of its own tree, leaves included (ADR-0005).

    The leaf rows are the half of it the link reads, and the writer used to leave
    them at zero.
    """
    with tempfile.TemporaryDirectory() as tmp:
        out, _, _ = _small_scenario(tmp)
        species, molecules = _indices_of(out)

        species_states = io.read_node_states(
            out / "simulated_Z_species.txt", molecules.n_leaves
        )
        molecule_states = io.read_node_states(
            out / "simulated_Z_molecules.txt", molecules.n_nodes
        )

    assert species_states.shape == (species.n_nodes, molecules.n_leaves)
    assert molecule_states.shape == (species.n_leaves, molecules.n_nodes)
    # Both leaf blocks carry states rather than the zeros the old writer left.
    assert species_states[species.leaves].any()
    assert molecule_states[:, molecules.leaves].any()


def test_the_field_is_the_and_of_the_two_tree_fields_it_was_drawn_from():
    """The written field, the written node states and the link agree.

    Every one of the three files is read back and the link's counters recomputed
    from them, so a transposed tree field or a mis-shaped node-state file shows up
    as a field density that the buckets cannot explain.
    """
    with tempfile.TemporaryDirectory() as tmp:
        out, config, meta = _small_scenario(tmp)
        species, molecules = _indices_of(out)
        field = io.read_field(out / "simulated_Y.txt", molecules.n_leaves)
        species_field = io.read_node_states(
            out / "simulated_Z_species.txt", molecules.n_leaves
        )[species.leaves]
        molecule_field = io.read_node_states(
            out / "simulated_Z_molecules.txt", molecules.n_nodes
        )[:, molecules.leaves]

    counters = L.link_counters(species_field, molecule_field, field)
    assert counters.tolist() == meta["link_counters"]
    assert field.mean() == pytest.approx(meta["field_ones_fraction"])

    # The rate at which the field reads 1 in each bucket, against the link's own
    # P_k. Bucket 0 is rare at a small omega, so only the buckets that hold cells
    # are judged.
    for bucket in range(L.N_BUCKETS):
        total = counters[bucket].sum()
        if total < 200:
            continue
        observed = counters[bucket, 1] / total
        predicted = float(L.prob_for_bucket(bucket, config.error_probability))
        assert observed == pytest.approx(predicted, abs=0.05), bucket


def test_the_simulate_parameters_carry_both_trees():
    """The replicate comparison runs the C++ under the reference's own draw, so
    every parameter of *both* trees has to be in the one file it is given."""
    with tempfile.TemporaryDirectory() as tmp:
        out, _, _ = _small_scenario(tmp)
        names = set(
            pd.read_csv(out / scenario.SIMULATE_PARAMETERS, sep="\t")["name"].astype(
                str
            )
        )

    for tree in ("species", "molecules"):
        for parameter in ("alpha", "log_nu", "branch_lengths"):
            assert any(n.startswith(f"{tree}_{parameter}_") for n in names), (
                f"{tree}_{parameter} is missing"
            )
        for scalar in ("mean_log_nu", "var_log_nu"):
            assert f"{tree}_{scalar}" in names


def test_no_run_script_neutralises_a_tree():
    """Neutralisation is retired. A rung that pinned one tree would reach the
    error probability through one tree where the model has two (ADR-0005), and
    the reference no longer needs it to be exact.

    The rung ladder itself is issue #44's; what this pins is that nothing the
    scenario writes still points at a neutralised molecules dimension.
    """
    with tempfile.TemporaryDirectory() as tmp:
        out, _, _ = _small_scenario(tmp)
        bodies = {path.name: path.read_text() for path in sorted(out.glob("*.sh"))}

    assert set(bodies) == {f"{name}.sh" for name, _, _, _ in scenario.RUNGS} | {
        "replicates.sh"
    }
    for name, body in bodies.items():
        assert "pinned_molecules" not in body, name
        if name == "replicates.sh":
            # The one script that legitimately hands the C++ the truth, because
            # both implementations have to run the same parameters.
            continue
        for parameter in ("alpha", "log_nu", "mean_log_nu", "var_log_nu"):
            assert f"--molecules_{parameter} " not in body, f"{name} pins {parameter}"


@pytest.mark.parametrize("omega", [0.05, 0.005])
def test_every_rung_infers_omega_under_a_prior_whose_mean_is_the_truth(omega):
    """The C++'s default prior on omega has mean 0.005. Left at that default
    against a truth of 0.05, the prior pulls omega's posterior down by about
    rate x variance, which goes from 0.0006 at rung 1 to about 0.04 at rung 3,
    and it takes alpha's coverage with it (#111)."""
    with tempfile.TemporaryDirectory() as tmp:
        out, _, _ = _small_scenario(tmp, error_probability=omega)
        bodies = {p.name: p.read_text() for p in sorted(out.glob("rung*.sh"))}

    for name, body in bodies.items():
        rates = re.findall(r"--error_probability_prior_rate (\S+)", body)
        assert len(rates) == 1, name
        # Truncation to (0, 0.5) moves the mean by less than 1e-4 at these rates.
        assert 1.0 / float(rates[0]) == pytest.approx(omega, rel=1e-12), name


@pytest.mark.parametrize("error_rate", [0.001, 0.01])
def test_every_rung_infers_epsilon_under_a_prior_whose_mean_is_the_truth(error_rate):
    """The C++'s default prior on epsilon is Beta(1, 1), whose mean 0.5 is far
    above a truth of 0.001. Once Y is free, LOTUS alone does not pin epsilon and
    omega apart, so that prior decides where the posterior sits (#113)."""
    with tempfile.TemporaryDirectory() as tmp:
        out, _, _ = _small_scenario(tmp, error_rate=error_rate)
        bodies = {p.name: p.read_text() for p in sorted(out.glob("rung*.sh"))}

    for name, body in bodies.items():
        priors = re.findall(r"--epsilon_prior_parameters (\S+),(\S+)", body)
        assert len(priors) == 1, name
        alpha, beta = (float(v) for v in priors[0])
        assert alpha / (alpha + beta) == pytest.approx(error_rate, rel=1e-12), name


@pytest.mark.parametrize("gamma", [1.1, 0.4])
def test_every_rung_infers_gamma_under_a_prior_whose_mean_is_the_truth(gamma):
    """Detection saturates once gamma is above about 10, so the likelihood is
    flat out there. Once Y is free, the chain walks into that plateau unless the
    prior closes it (#114). Gamma(shape, rate) has mean shape / rate."""
    with tempfile.TemporaryDirectory() as tmp:
        out, _, _ = _small_scenario(tmp, gamma=gamma)
        bodies = {p.name: p.read_text() for p in sorted(out.glob("rung*.sh"))}

    for name, body in bodies.items():
        priors = re.findall(r"--gamma_prior_parameters (\S+),(\S+)", body)
        assert len(priors) == 1, name
        shape, rate = (float(v) for v in priors[0])
        assert shape / rate == pytest.approx(gamma, rel=1e-12), name


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
    assert any(
        m in scenario.SIMULATE_PARAMETERS for m in scenario.INITIAL_VALUE_MARKERS
    ), scenario.SIMULATE_PARAMETERS


def test_research_effort_uses_log_paper_counts():
    """The C++ log-transforms paper counts when reading them (TTree.h:517)."""
    papers = np.array([0, 1, 3, 9])
    gamma = 1.1
    expected = 1.0 - np.exp(-gamma * np.log(papers + 1.0))
    assert np.allclose(
        research_effort(papers, np.array([3]), gamma)[:, 0],
        expected * (1.0 - np.exp(-gamma * np.log(4.0))),
    )


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


# ---------------------------------------------------------------------------
# The joint sums to one
# ---------------------------------------------------------------------------
#
# ADR-0005 argues that the joint sums to one by construction. These enumerate a
# model small enough to list every (Z_s, Z_m, Y) and check it. Both trees stay
# active, for the reason ADR-0005 gives.
#
# The total reads the two node-state densities, and it reads whether a variable
# is scored twice. It is blind to the shape of the link and to the orientation of
# a leaf block, which the field rate and the reference draw below pin instead.
# `enumeration` says why.

# Each entry is (species tree, molecule tree), as (node count, type). The shapes
# differ on purpose: a species tree with two roots, and a pair whose leaf counts
# disagree. Enumerability caps the trees at a handful of leaves, so the third
# pair is lopsided rather than large.
SMALL_MODELS = [
    ((3, TreeType.balanced), (3, TreeType.balanced)),
    ((4, TreeType.grass), (3, TreeType.balanced)),
    ((4, TreeType.star), (2, TreeType.grass)),
]


def _process(index, n_cliques: int, seed: int) -> E.TreeProcess:
    """One tree's process, with an alpha and a nu per clique."""
    rng = np.random.default_rng(seed)
    return E.TreeProcess(
        index=index,
        bins=F.sample_binned_branch_lengths(rng, index.n_branches),
        alphas=rng.uniform(0.2, 0.8, n_cliques),
        nus=np.exp(rng.normal(-0.5, 0.5, n_cliques)),
    )


def _both_trees(model, seed: int = 20260907):
    """Both trees of one model. A clique of each is named by a leaf of the other."""
    species_spec, molecule_spec = model
    species_index = build_tree_index(edges_of(Tree(*species_spec, "species")))
    molecule_index = build_tree_index(edges_of(Tree(*molecule_spec, "molecules")))
    return (
        _process(species_index, molecule_index.n_leaves, seed),
        _process(molecule_index, species_index.n_leaves, seed + 1),
    )


@pytest.mark.parametrize("model", SMALL_MODELS)
@pytest.mark.parametrize("omega", [0.02, 0.15, 0.4])
def test_the_joint_sums_to_one(model, omega):
    """ADR-0005's claim, enumerated rather than argued."""
    species, molecules = _both_trees(model)
    mass = E.joint(species, molecules, omega).mass
    assert mass == pytest.approx(1.0, abs=1e-12), f"the joint sums to {mass!r}"


@pytest.mark.parametrize("model", SMALL_MODELS)
def test_neither_tree_is_neutral_where_the_joint_is_summed(model):
    """The parameters the sum runs under, pinned. A neutral tree carries no
    phylogenetic signal, and the sum above must not be read under one."""
    for process in _both_trees(model):
        rows = F.transition_matrices(
            process.alphas, process.nus, F.grid_branch_lengths()
        )
        assert np.abs(rows - 0.5).max() > 0.05
        assert np.abs(process.alphas - 0.5).max() > 0.05


@pytest.mark.parametrize("model", SMALL_MODELS)
def test_each_node_state_is_a_distribution(model):
    """The two tree factors, each on its own. A transposed transition matrix
    leaves one of them, which is what the total then reports."""
    for process in _both_trees(model):
        _, probability = E.node_state_distribution(process)
        assert probability.sum() == pytest.approx(1.0)
        assert (probability > 0.0).all()


@pytest.mark.parametrize("model", SMALL_MODELS)
def test_summing_the_field_out_returns_the_two_node_states(model):
    """The link is a conditional distribution at every pair of tree fields."""
    species, molecules = _both_trees(model)
    enumeration = E.joint(species, molecules, 0.1)
    assert enumeration.table.sum(axis=2) == pytest.approx(
        np.outer(enumeration.species_probability, enumeration.molecule_probability)
    )


def test_a_factor_counted_twice_leaves_the_sum():
    """The defect class ADR-0002 records is a variable scored twice, and the sum
    is what sees it. A check that cannot fail is not a check."""
    species, molecules = _both_trees(SMALL_MODELS[0])
    tree_fields = E.species_tree_fields(species)

    # One species leaf, scored a second time under its own stationary rate.
    alpha = species.alphas[0]
    cell = tree_fields.tree_fields[:, 0, 0]
    twice = tree_fields.probability * np.where(cell, alpha, 1.0 - alpha)

    deformed = E.joint_over_tree_fields(
        E.TreeFieldDistribution(tree_fields.tree_fields, twice),
        E.molecule_tree_fields(molecules),
        0.1,
    )
    assert deformed.mass == pytest.approx(alpha**2 + (1.0 - alpha) ** 2)
    assert deformed.mass < 0.99


def test_the_enumerated_field_rate_is_the_product_of_the_two_adjusted_rates():
    """ADR-0005, derivation 3, read off the enumeration rather than a simulation.

    A species clique is named by a molecule leaf and a molecule clique by a
    species leaf, so the two alphas enter the cell at (species leaf, molecule
    leaf) from opposite sides. This is one of the two tests that pin the
    orientation of a leaf block; the sum itself cannot.
    """
    species, molecules = _both_trees(SMALL_MODELS[1])
    omega = 0.12
    assert E.joint(species, molecules, omega).field_rate == pytest.approx(
        np.outer(
            L.adjusted_rate(molecules.alphas, omega),
            L.adjusted_rate(species.alphas, omega),
        )
    )


def test_the_enumerated_field_matches_the_reference_draw():
    """The enumeration and the harness's own draw describe the same model.

    Without this the joint could sum to one and still be a statement about a
    model that nothing simulates.
    """
    species, molecules = _both_trees(SMALL_MODELS[0])
    omega, n_draws = 0.1, 20_000
    enumeration = E.joint(species, molecules, omega)

    rng = np.random.default_rng(4)
    counts = np.zeros(len(enumeration.fields))
    for _ in range(n_draws):
        z_s = F.sample_states(
            rng, species.index, species.bins, species.alphas, species.nus
        )
        z_m = F.sample_states(
            rng, molecules.index, molecules.bins, molecules.alphas, molecules.nus
        )
        field = L.sample_field(
            rng, z_s[species.index.leaves], z_m[molecules.index.leaves].T, omega
        )
        counts[E.code_of(field)] += 1

    assert counts / n_draws == pytest.approx(
        enumeration.table.sum(axis=(0, 1)), abs=0.01
    )


@pytest.mark.parametrize("shape", [(2, 2), (3, 1), (1, 4)])
def test_a_code_names_the_array_the_enumeration_lists(shape):
    """`code_of` inverts `enumerate_binary_arrays`, so the draw above and the
    enumeration index a field the same way."""
    arrays = E.enumerate_binary_arrays(shape)
    assert [E.code_of(array) for array in arrays] == list(range(len(arrays)))


def test_a_clique_that_names_no_leaf_of_the_other_tree_is_refused():
    species, molecules = _both_trees(SMALL_MODELS[0])
    one_clique = E.TreeProcess(
        index=species.index,
        bins=species.bins,
        alphas=species.alphas[:1],
        nus=species.nus[:1],
    )
    with pytest.raises(ValueError, match="named by a molecule leaf"):
        E.joint(one_clique, molecules, 0.1)


def test_a_model_beyond_enumeration_is_refused():
    with pytest.raises(ValueError, match="intractable"):
        E.enumerate_binary_arrays((8, 4))


# --------------------------------------------------------------------------
# The LOTUS+MS scenarios
# --------------------------------------------------------------------------


def _ms_scenario(tmp: str, layout: str = "dense", **overrides):
    """A small seeded scenario with mass-spec data.

    The molecule tree needs at least as many leaves as the longest candidate list,
    so it is larger than `_small_scenario`'s.
    """
    out = pathlib.Path(tmp) / "scenario"
    config = scenario.ScenarioConfig(
        n_species_nodes=63,
        n_molecule_nodes=63,
        mass_spec=mass_spec.LAYOUTS[layout],
        **overrides,
    )
    return out, config, scenario.build_scenario(out, config)


def _ms_files(out: pathlib.Path):
    data = pd.read_csv(out / mass_spec.DATA_FILE, sep="\t")
    truth = pd.read_csv(out / mass_spec.ASSIGNMENT_FILE, sep="\t")
    return data, truth


FEATURE = ["species", "extract", "feature"]


def test_every_molecule_feature_lists_its_true_molecule():
    with tempfile.TemporaryDirectory() as tmp:
        out, _, _ = _ms_scenario(tmp)
        data, truth = _ms_files(out)

    # Every feature that is not unknown has exactly one true row, and that row is
    # one of the structures the data file lists for it.
    unknown = truth[truth["molecule"] == mass_spec.UNKNOWN]
    molecule_features = unknown[unknown["fraction"] == 0][FEATURE]
    assert len(molecule_features) > 0
    assigned = truth[
        (truth["fraction"] == 1) & (truth["molecule"] != mass_spec.UNKNOWN)
    ].merge(molecule_features, on=FEATURE)
    assert len(assigned) == len(molecule_features)
    listed = data.merge(assigned, on=FEATURE + ["molecule"])
    assert len(listed) == len(assigned)


def test_no_molecule_explains_two_features_of_one_run():
    with tempfile.TemporaryDirectory() as tmp:
        out, _, _ = _ms_scenario(tmp)
        _, truth = _ms_files(out)

    assigned = truth[(truth["fraction"] == 1) & (truth["molecule"] != mass_spec.UNKNOWN)]
    assert not assigned.duplicated(["species", "extract", "molecule"]).any()


def test_each_feature_has_exactly_one_explanation_and_unknown_features_list_only_wrong_rows():
    with tempfile.TemporaryDirectory() as tmp:
        out, _, _ = _ms_scenario(tmp)
        _, truth = _ms_files(out)

    per_feature = truth.groupby(FEATURE)["fraction"].sum()
    assert (per_feature == 1).all()
    unknown = truth[(truth["molecule"] == mass_spec.UNKNOWN) & (truth["fraction"] == 1)]
    assert len(unknown) > 0
    rows = truth.merge(unknown[FEATURE], on=FEATURE)
    wrong = rows[rows["molecule"] != mass_spec.UNKNOWN]
    assert len(wrong) > 0
    assert (wrong["fraction"] == 0).all()


def test_the_number_of_proposed_structures_is_between_1_and_20_for_both_kinds():
    with tempfile.TemporaryDirectory() as tmp:
        out, _, _ = _ms_scenario(tmp)
        data, truth = _ms_files(out)

    n_f = data.groupby(FEATURE).size().rename("n_f").reset_index()
    is_unknown = (
        truth[truth["molecule"] == mass_spec.UNKNOWN][FEATURE + ["fraction"]]
        .rename(columns={"fraction": "unknown"})
    )
    merged = n_f.merge(is_unknown, on=FEATURE)
    for kind in (0, 1):
        sizes = merged[merged["unknown"] == kind]["n_f"]
        assert len(sizes) > 0
        assert sizes.min() >= 1
        assert sizes.max() <= 20
    assert (merged["n_f"] == 20).any() and (merged["n_f"] == 1).any()


def test_rho_is_lambda_over_the_molecule_leaves_and_lambda_follows_the_true_field():
    """lambda is recomputed here from the written field, filter probabilities and
    kappa, none of which the simulator hands back."""
    with tempfile.TemporaryDirectory() as tmp:
        out, _, _ = _ms_scenario(tmp)
        species, molecules = _indices_of(out)
        field = io.read_field(out / "simulated_Y.txt", molecules.n_leaves)
        values = pd.read_csv(out / scenario.SIMULATE_PARAMETERS, sep="\t")
        values = dict(zip(values["name"].astype(str), values["value"].astype(float)))
        data, _ = _ms_files(out)

    index = np.array([values[f"filter_proba_{n}"] for n in molecules.leaf_names()])
    kappa = values["contamination_proba"]
    assert kappa == 0.01
    assert ((index >= 0) & (index <= 254) & (index == index.astype(int))).all()

    # Dense: every species leaf has runs, so every run's expectation counts.
    expected = np.where(field, index / 255.0, kappa).sum(axis=1).mean()
    assert values["ms_lambda"] == pytest.approx(expected)
    assert values["ms_unknown_rate"] == pytest.approx(
        values["ms_lambda"] / molecules.n_leaves
    )


def test_the_true_assignment_is_keyed_by_the_features_of_the_data_file():
    with tempfile.TemporaryDirectory() as tmp:
        out, _, _ = _ms_scenario(tmp)
        data, truth = _ms_files(out)

    data_rows = set(map(tuple, data[FEATURE + ["molecule"]].to_numpy()))
    truth_rows = set(map(tuple, truth[FEATURE + ["molecule"]].to_numpy()))
    unknown_rows = {r for r in truth_rows if r[3] == mass_spec.UNKNOWN}
    assert truth_rows - unknown_rows == data_rows
    assert {r[:3] for r in unknown_rows} == {r[:3] for r in data_rows}


def test_about_half_of_all_features_are_unknown():
    with tempfile.TemporaryDirectory() as tmp:
        out, _, _ = _ms_scenario(tmp)
        _, truth = _ms_files(out)

    unknown = truth[truth["molecule"] == mass_spec.UNKNOWN]["fraction"]
    assert unknown.mean() == pytest.approx(0.5, abs=0.05)


def test_the_dense_layout_has_two_runs_per_species_leaf():
    with tempfile.TemporaryDirectory() as tmp:
        out, _, _ = _ms_scenario(tmp, "dense")
        species, _ = _indices_of(out)
        data, _ = _ms_files(out)

    runs = data[["species", "extract"]].drop_duplicates()
    assert set(runs["species"]) == set(species.leaf_names())
    assert (runs.groupby("species").size() == 2).all()


def test_the_sparse_layout_has_eight_species_with_one_run_each():
    with tempfile.TemporaryDirectory() as tmp:
        out, _, _ = _ms_scenario(tmp, "sparse")
        data, _ = _ms_files(out)

    runs = data[["species", "extract"]].drop_duplicates()
    assert runs["species"].nunique() == 8
    assert (runs.groupby("species").size() == 1).all()


def test_the_ms_rungs_hand_the_binary_the_data_and_the_true_rate_and_nothing_else_of_the_truth():
    with tempfile.TemporaryDirectory() as tmp:
        out, _, meta = _ms_scenario(tmp)
        bodies = {p.name: p.read_text() for p in sorted(out.glob("rung*.sh"))}

    assert set(bodies) == {f"{name}.sh" for name, _, _, _ in scenario.MS_RUNGS}
    rho = meta["mass_spec"]["unknown_rate"]
    for name, body in bodies.items():
        assert f"--ms_data {mass_spec.DATA_FILE}" in body, name
        assert f"--ms_unknown_rate {rho!r}" in body, name
        # The scenarios have no outside structures, so tau is 0 in every rung.
        assert "--ms_outside_rate 0 " in body, name
        # Only the data file is read: the assignment, the filter probabilities and
        # kappa are inferred in every rung.
        for truth in (mass_spec.ASSIGNMENT_FILE, scenario.SIMULATE_PARAMETERS):
            assert truth not in body, name
        for parameter in (mass_spec.FILTER_PARAMETER, mass_spec.CONTAMINATION_PARAMETER):
            assert f"--{parameter}" not in body, name
        assert 'FLAGS="lm"' in body, name


def test_a_scenario_without_mass_spec_writes_no_mass_spec_files_or_flags():
    with tempfile.TemporaryDirectory() as tmp:
        out, _, _ = _small_scenario(tmp)
        names = {p.name for p in out.iterdir()}
        bodies = [p.read_text() for p in out.glob("*.sh")]

    assert mass_spec.DATA_FILE not in names
    assert mass_spec.ASSIGNMENT_FILE not in names
    assert not any("--ms_data" in body for body in bodies)


def _with_options(**options):
    return {"ms_options": mass_spec.MassSpecOptions(**options)}


def test_the_longest_candidate_list_is_an_option():
    with tempfile.TemporaryDirectory() as tmp:
        out, _, _ = _ms_scenario(tmp, **_with_options(max_candidates=3))
        data, _ = _ms_files(out)

    lengths = data.groupby(FEATURE).size()
    assert lengths.max() == 3
    assert lengths.min() == 1


def test_the_written_incorrect_density_is_an_option():
    with tempfile.TemporaryDirectory() as tmp:
        out, _, _ = _ms_scenario(tmp, **_with_options(log_incorrect_density=-1.5))
        data, _ = _ms_files(out)

    assert (data["log_incorrect_density"] == -1.5).all()


def test_the_correct_score_distribution_is_an_option():
    # Beta(1, 1) is flat, so the log correct density is 0 at every score.
    with tempfile.TemporaryDirectory() as tmp:
        out, _, _ = _ms_scenario(tmp, **_with_options(correct_shape=1.0))
        data, _ = _ms_files(out)

    assert (data["log_correct_density"] == 0.0).all()


def test_the_incorrect_score_distribution_is_an_option():
    # Beta(1000, 1) puts a wrong row's log score within 1e-2 of 0 unless u < e**-10,
    # where the Beta(2, 1) log density is ln 2.
    with tempfile.TemporaryDirectory() as tmp:
        out, _, _ = _ms_scenario(
            tmp, **_with_options(correct_shape=2.0, incorrect_shape=1000.0)
        )
        data, truth = _ms_files(out)

    true_rows = truth[
        (truth["fraction"] == 1) & (truth["molecule"] != mass_spec.UNKNOWN)
    ][FEATURE + ["molecule"]]
    marked = data.merge(true_rows, on=FEATURE + ["molecule"], how="left", indicator=True)
    wrong = marked[marked["_merge"] == "left_only"]
    assert len(wrong) > 0
    assert np.allclose(wrong["log_correct_density"], np.log(2.0), atol=1e-2)
