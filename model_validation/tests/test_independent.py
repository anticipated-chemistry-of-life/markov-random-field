"""Verification of the independent reference implementation.

Without these, an unverified Python would be judging an unverified C++, and a
disagreement between them would say nothing about which side is wrong.
"""

from __future__ import annotations

import pathlib
import tempfile

import numpy as np
import pandas as pd
import pytest
from scipy.linalg import expm

from src.independent import field as F
from src.independent import io
from src.independent import link as L
from src.independent import scenario
from src.independent import toy_normaliser as TN
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


def test_the_leaf_block_keeps_file_order():
    """A reordered tree file must give a correspondingly reordered output, not an
    unrelated one -- so within a block, first-appearance order survives."""
    tree = Tree(255, TreeType.balanced, "species")
    edges = edges_of(tree)
    index = build_tree_index(edges)

    appearance: list[str] = []
    for child, parent in edges:
        for name in (parent, child):
            if name not in appearance:
                appearance.append(name)
    rank = {name: i for i, name in enumerate(appearance)}

    leaf_ranks = [rank[name] for name in index.leaf_names()]
    assert leaf_ranks == sorted(leaf_ranks)


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
# The normalising constant of the two-tree product
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("depth", [1, 2])
@pytest.mark.parametrize("alpha,nu", [(0.5, 0.6), (0.2, 2.0), (0.8, 0.1)])
def test_leaf_patterns_are_a_distribution(depth, alpha, nu):
    patterns = TN.leaf_pattern_probabilities(alpha, nu, depth)
    assert len(patterns) == 2 ** (2**depth)
    assert patterns.sum() == pytest.approx(1.0)
    assert (patterns > 0.0).all()


@pytest.mark.parametrize("depth", [1, 2])
def test_leaf_marginal_is_stationary(depth):
    """Marginalising all but one leaf must return the stationary probability."""
    alpha, n_leaves = 0.3, 2**depth
    patterns = TN.leaf_pattern_probabilities(alpha, 0.7, depth)
    codes = np.arange(len(patterns))
    for leaf in range(n_leaves):
        on = ((codes >> leaf) & 1).astype(bool)
        assert patterns[on].sum() == pytest.approx(alpha)


@pytest.mark.parametrize("depth", [1, 2])
def test_pattern_codes_agree_with_a_transpose(depth):
    """Row codes of a field are the column codes of its transpose."""
    rows, cols = TN.pattern_codes(depth)
    n_leaves = 2**depth
    bits = (
        (np.arange(len(rows))[:, None] >> np.arange(n_leaves * n_leaves)[None, :]) & 1
    ).reshape(-1, n_leaves, n_leaves)
    powers = 2 ** np.arange(n_leaves)
    assert (cols == (bits.transpose(0, 2, 1) * powers[None, None, :]).sum(2)).all()
    assert rows.shape == cols.shape == (len(rows), n_leaves)


@pytest.mark.parametrize("depth", [1, 2])
def test_neutral_molecules_make_the_constant_independent_of_species(depth):
    """ADR-0001's assumption, stated as an equation.

    With neutral molecules `C` must not move when the species parameters do —
    that is exactly why the independent-field harness is unbiased and also why
    it cannot see the effect this module isolates.
    """
    neutral = TN.leaf_pattern_probabilities(
        0.5, np.exp(F.STATIONARY_NU_THRESHOLD), depth
    )
    constants = [
        TN.normalising_constant(
            TN.leaf_pattern_probabilities(0.4, nu, depth), neutral, depth
        )
        for nu in (0.1, 0.5, 2.0, 10.0)
    ]
    assert constants == pytest.approx([constants[0]] * len(constants))


@pytest.mark.parametrize("depth", [1, 2])
def test_non_neutral_molecules_make_the_constant_move(depth):
    """The converse: without neutrality the omitted term is not a constant."""
    molecules = TN.leaf_pattern_probabilities(0.5, 0.2, depth)
    constants = [
        TN.normalising_constant(
            TN.leaf_pattern_probabilities(0.5, nu, depth), molecules, depth
        )
        for nu in (0.1, 0.5, 2.0, 10.0)
    ]
    assert max(constants) > min(constants) * 1.05


@pytest.mark.parametrize("depth", [1, 2])
def test_field_distribution_sums_to_one(depth):
    species = TN.leaf_pattern_probabilities(0.4, 0.9, depth)
    molecules = TN.leaf_pattern_probabilities(0.6, 0.3, depth)
    assert TN.field_distribution(species, molecules, depth).sum() == pytest.approx(1.0)


@pytest.mark.parametrize("depth", [1, 2])
def test_correct_objective_is_maximised_at_the_truth(depth):
    """Gibbs' inequality, which is what licenses reading the gap as bias."""
    true_log_nu = -0.5
    molecules = TN.leaf_pattern_probabilities(0.5, np.exp(-1.0), depth)
    truth = TN.field_distribution(
        TN.leaf_pattern_probabilities(0.5, np.exp(true_log_nu), depth), molecules, depth
    )
    grid = np.linspace(-4.0, 2.0, 121)
    _, correct = TN.expected_log_likelihood_profile(truth, grid, 0.5, molecules, depth)
    assert grid[int(np.argmax(correct))] == pytest.approx(true_log_nu, abs=0.05)


@pytest.mark.parametrize("depth", [1, 2])
def test_targeted_objective_is_biased_downward_off_neutrality(depth):
    """The C++'s objective peaks below the truth once molecules are active."""
    true_log_nu = -0.5
    molecules = TN.leaf_pattern_probabilities(0.5, np.exp(-2.0), depth)
    truth = TN.field_distribution(
        TN.leaf_pattern_probabilities(0.5, np.exp(true_log_nu), depth), molecules, depth
    )
    grid = np.linspace(-4.0, 2.0, 121)
    targeted, correct = TN.expected_log_likelihood_profile(
        truth, grid, 0.5, molecules, depth
    )
    assert grid[int(np.argmax(targeted))] < grid[int(np.argmax(correct))] - 0.2


def test_the_two_objectives_coincide_under_neutrality():
    """No bias to find when the molecules dimension is switched off."""
    depth = 2
    neutral = TN.leaf_pattern_probabilities(
        0.5, np.exp(F.STATIONARY_NU_THRESHOLD), depth
    )
    truth = TN.field_distribution(
        TN.leaf_pattern_probabilities(0.5, np.exp(-0.5), depth), neutral, depth
    )
    grid = np.linspace(-3.0, 1.0, 41)
    targeted, correct = TN.expected_log_likelihood_profile(
        truth, grid, 0.5, neutral, depth
    )
    assert np.argmax(targeted) == np.argmax(correct)


def test_correcting_the_chain_removes_the_drift():
    """End to end: the same chain, drifting and not drifting."""
    shared = dict(
        true_log_nu_species=-0.5,
        log_nu_molecules=-2.0,
        alpha_species=0.5,
        alpha_molecules=0.5,
        depth=2,
        n_iterations=600,
    )
    drifting = TN.run_chain(
        np.random.default_rng(7), correct_normaliser=False, **shared
    )
    corrected = TN.run_chain(
        np.random.default_rng(7), correct_normaliser=True, **shared
    )
    assert corrected[-300:].mean() == pytest.approx(-0.5, abs=0.2)
    assert drifting[-300:].mean() < corrected[-300:].mean() - 0.4


def test_depth_beyond_enumeration_is_refused():
    with pytest.raises(ValueError, match="intractable"):
        TN.pattern_codes(TN.MAX_ENUMERABLE_DEPTH + 1)
