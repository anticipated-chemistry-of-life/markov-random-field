"""Assemble a complete validation scenario: truth, observations, and run scripts.

Both trees are active. Each draws its own node state top-down over every one of
its nodes, the leaf block of that node state is that tree's **tree field**, and
the field is a noisy AND of the two (`link`). Neutralising the molecules
dimension is no longer the price of having an independent reference: under
ADR-0005 the reference is exact with both trees running, which is precisely what
ADR-0001 could not do.

The reference shares no code with the C++ binary. It re-derives the draw from the
maths, so a disagreement between the two is a bug in one of them.

The three rungs are run in the order written and abandoned at the first failure.
Each pins strictly less than the one before, so a failure localises the fault:

1. field *and* node states pinned — only the parameter updates move, so this is
   close to closed form and sets the empirical ceiling for the other two.
2. field pinned, node states inferred — adds the node-state Gibbs update.
3. nothing pinned, inferred from observations — the production path, run against
   the simple error model alone, then LOTUS alone, then both.
"""

from __future__ import annotations

import json
import pathlib
from dataclasses import dataclass

import numpy as np

from ..tree import Tree, TreeType
from . import field as F
from . import io
from . import link as L
from .data import sample_paper_counts, simulate_lotus, simulate_simple_error
from .indexing import TreeIndex, build_tree_index

DEFAULT_REPLICATES = 20

# The `simulated` substring in both names is required by stattools' filename-
# driven reader dispatch; see the comment in `build_scenario`. `test_independent`
# asserts it, so a rename fails in Python rather than mysteriously in C++.
INITIAL_VALUE_MARKERS = (
    "trace",
    "simulated",
    "meanVar",
    "statePosteriors",
    "posteriorMode",
)
SIMULATE_PARAMETERS = "simulated_parameters.txt"


@dataclass
class ScenarioConfig:
    seed: int = 42
    n_species_nodes: int = 255
    n_molecule_nodes: int = 255
    mean_log_nu: float = -0.5
    var_log_nu: float = 0.25
    error_probability: float = 0.05
    """Omega: the rate at which a tree field cell is corrupted before the AND."""
    epsilon: float = 0.05
    gamma: float = 1.1
    error_rate: float = 0.001
    iterations: int = 10_000
    burnin: int = 1_000
    n_burnin_rounds: int = 10
    true_branch_lengths: bool = False
    """Write grid centres into both tree files instead of flat lengths.

    Flat is the default: it makes the chain find the branch lengths rather than
    start on them, so 'never left the initialisation' cannot pass as recovery.
    """


@dataclass
class TreeTruth:
    """One tree's process, and the node state it drew under it."""

    index: TreeIndex
    alphas: np.ndarray
    """One per clique, and a clique of this tree is named by a leaf of the other."""
    log_nus: np.ndarray
    bins: np.ndarray
    states: np.ndarray
    """`(n_nodes, n_cliques)`, every node of this tree."""

    @property
    def tree_field(self) -> np.ndarray:
        """The leaf block: this tree's own view of the leaf-level field."""
        return self.states[self.index.leaves]


def _edges(tree: Tree) -> list[tuple[str, str]]:
    frame = tree.to_dataframe()
    return list(zip(frame["child"].astype(str), frame["parent"].astype(str)))


def _draw_tree(
    rng: np.random.Generator,
    index: TreeIndex,
    n_cliques: int,
    name: str,
    config: ScenarioConfig,
) -> TreeTruth:
    """Draw one tree's parameters, then its whole node state under them."""
    alphas = rng.beta(0.5, 0.5, size=n_cliques)
    log_nus = rng.normal(config.mean_log_nu, np.sqrt(config.var_log_nu), size=n_cliques)
    nus = np.exp(log_nus)
    if nus.max() >= F.STATIONARY_NU_THRESHOLD:
        raise ValueError(
            f"A {name} clique drew nu = {nus.max():.2f}, at or past the C++'s "
            f"stationary short-circuit ({F.STATIONARY_NU_THRESHOLD}). The reference "
            "does not replicate that approximation, so the two would diverge "
            "legitimately. Re-seed or lower var_log_nu."
        )

    bins = F.sample_binned_branch_lengths(rng, index.n_branches)
    return TreeTruth(
        index=index,
        alphas=alphas,
        log_nus=log_nus,
        bins=bins,
        states=F.sample_states(rng, index, bins, alphas, nus),
    )


def build_scenario(out: pathlib.Path, config: ScenarioConfig) -> dict:
    rng = np.random.default_rng(config.seed)
    out.mkdir(parents=True, exist_ok=True)
    omega = L.check_error_probability(config.error_probability)

    species_tree = Tree(config.n_species_nodes, TreeType.balanced, "species")
    molecule_tree = Tree(config.n_molecule_nodes, TreeType.balanced, "molecules")
    species_edges, molecule_edges = _edges(species_tree), _edges(molecule_tree)
    species = build_tree_index(species_edges)
    molecules = build_tree_index(molecule_edges)

    # -- the two trees, each drawn top-down over every one of its nodes ---------
    # A species-tree clique is named by a molecule leaf, and a molecule-tree
    # clique by a species leaf. So each tree's node state is
    # (its own nodes) x (the other tree's leaves).
    species_truth = _draw_tree(rng, species, molecules.n_leaves, "species", config)
    molecule_truth = _draw_tree(rng, molecules, species.n_leaves, "molecules", config)

    # Both tree fields are indexed [species_leaf, molecule_leaf]; the molecule
    # tree drew its own nodes down the rows, so its leaf block is transposed into
    # that shape. A tree field and the field share their subscripts (ADR-0005).
    species_field = species_truth.tree_field
    molecule_field = molecule_truth.tree_field.T

    # -- the field: a noisy AND of the two tree fields -------------------------
    latent_field = L.sample_field(rng, species_field, molecule_field, omega)

    # -- observations ---------------------------------------------------------
    species_papers = sample_paper_counts(rng, species.n_leaves)
    molecule_papers = sample_paper_counts(rng, molecules.n_leaves)
    simple_data = simulate_simple_error(rng, latent_field, config.epsilon)
    lotus = simulate_lotus(
        rng,
        latent_field,
        species_papers,
        molecule_papers,
        config.gamma,
        config.error_rate,
    )

    # -- write ----------------------------------------------------------------
    def lengths_of(truth: TreeTruth) -> np.ndarray:
        if config.true_branch_lengths:
            return F.grid_branch_lengths()[truth.bins]
        return np.full(truth.index.n_branches, 0.2)

    io.write_tree(
        out / "species.txt", species_edges, species, lengths_of(species_truth)
    )
    io.write_tree(
        out / "molecules.txt", molecule_edges, molecules, lengths_of(molecule_truth)
    )
    io.write_paper_counts(
        out / "species_papers.txt", "species", species.leaf_names(), species_papers
    )
    io.write_paper_counts(
        out / "molecules_papers.txt",
        "molecules",
        molecules.leaf_names(),
        molecule_papers,
    )

    # The truth is written where no inference run is ever pointed at it: a stray
    # flag would otherwise start the chain on the answer it is meant to find. The
    # one file that legitimately hands truth to the C++ is SIMULATE_PARAMETERS,
    # because the replicate comparison needs *both* implementations running the
    # same parameters. It is never referenced by an infer script.
    #
    # The `simulated_` prefix there is load-bearing, not decorative. stattools
    # picks a reader by *filename*: a name/value file is only matched up by
    # parameter name when the filename contains one of trace/simulated/meanVar/
    # statePosteriors/posteriorMode, and otherwise has to be a bare one-column
    # file of exactly the right length (TReadInitialValues.h:133).
    species_values = _tree_truth_values("species", species_truth, molecules, config)
    molecule_values = _tree_truth_values("molecules", molecule_truth, species, config)
    io.write_parameters(out / "truth_species.txt", species_values)
    io.write_parameters(out / "truth_molecules.txt", molecule_values)
    io.write_parameters(
        out / SIMULATE_PARAMETERS, {**species_values, **molecule_values}
    )

    io.write_field(out / "simulated_Y.txt", latent_field, species, molecules)
    io.write_node_states(
        out / "simulated_Z_species.txt", species_truth.states, species, molecules, 0
    )
    io.write_node_states(
        out / "simulated_Z_molecules.txt",
        molecule_truth.states.T,
        species,
        molecules,
        1,
    )
    io.write_observations(out / "simulated_lotus.tsv", lotus, species, molecules)
    io.write_observations(
        out / "simulated_simple_data.tsv", simple_data, species, molecules
    )

    counters = L.link_counters(species_field, molecule_field, latent_field)
    meta = {
        "seed": config.seed,
        "n_bins": F.N_BINS,
        "n_species_leaves": species.n_leaves,
        "n_molecule_leaves": molecules.n_leaves,
        "n_species_branches": species.n_branches,
        "n_molecule_branches": molecules.n_branches,
        "mean_log_nu": config.mean_log_nu,
        "var_log_nu": config.var_log_nu,
        "error_probability": omega,
        "epsilon": config.epsilon,
        "gamma": config.gamma,
        "error_rate": config.error_rate,
        "iterations": config.iterations,
        "burnin": config.burnin,
        "n_burnin_rounds": config.n_burnin_rounds,
        "true_branch_lengths": config.true_branch_lengths,
        "field_ones_fraction": float(latent_field.mean()),
        "species_tree_field_ones_fraction": float(species_field.mean()),
        "molecule_tree_field_ones_fraction": float(molecule_field.mean()),
        # The expected density is the product of the two adjusted rates, and
        # nothing else (ADR-0005, derivation 3). It is what the field's own
        # density can be checked against, and all it can be checked against.
        "expected_field_ones_fraction": float(
            L.adjusted_rate(species_truth.alphas, omega).mean()
            * L.adjusted_rate(molecule_truth.alphas, omega).mean()
        ),
        "link_counters": counters.tolist(),
        "lotus_records": int(lotus.sum()),
        "simple_data_disagreements": int((simple_data != latent_field).sum()),
        "species_bins": species_truth.bins.tolist(),
        "molecule_bins": molecule_truth.bins.tolist(),
    }
    (out / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")

    _write_scripts(out, config)
    return meta


def _tree_truth_values(
    tree_name: str,
    truth: TreeTruth,
    other: TreeIndex,
    config: ScenarioConfig,
) -> dict[str, object]:
    """One tree's answers, named the way the C++ names its parameters.

    A clique of this tree is named by a leaf of the *other* tree, which is what
    `other` is here for. Branches are named by the child node they hang below,
    and those belong to this tree.
    """
    values: dict[str, object] = {}
    for name, value in zip(truth.index.branch_names(), truth.bins):
        values[f"{tree_name}_branch_lengths_{name}"] = int(value)
    values[f"{tree_name}_mean_log_nu"] = config.mean_log_nu
    values[f"{tree_name}_var_log_nu"] = config.var_log_nu
    for name, value in zip(other.leaf_names(), truth.log_nus):
        values[f"{tree_name}_log_nu_{name}"] = float(value)
    for name, value in zip(other.leaf_names(), truth.alphas):
        values[f"{tree_name}_alpha_{name}"] = float(value)
    return values


_PREAMBLE = """\
#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
ROOT="$SCRIPT_DIR/../.."

MODE="${{ACOL_MODE:-release}}"
FLAGS="{flags}"

cd "$ROOT"
just build "$MODE" "$FLAGS"
ACOL="$ROOT/$(just bin "$MODE" "$FLAGS")"

cd "$SCRIPT_DIR"
"""


def _infer_script(
    name: str,
    flags: str,
    extra: str,
    config: ScenarioConfig,
    seed_offset: int,
    threads: str = "all",
) -> str:
    """A complete, standalone script: preamble plus one inference run."""
    return _PREAMBLE.format(flags=flags) + _infer_command(
        name, extra, config, seed_offset, threads
    )


def _infer_command(
    name: str,
    extra: str,
    config: ScenarioConfig,
    seed_offset: int,
    threads: str = "all",
) -> str:
    """Just the invocation, so several can share one preamble."""
    return (
        f'mkdir -p "{name}"\n'
        f'"$ACOL" infer \\\n'
        f"    --out ./{name}/acol \\\n"
        f"    --tree_species species.txt \\\n"
        f"    --tree_molecules molecules.txt \\\n"
        f"    --species_paper_counts species_papers.txt \\\n"
        f"    --molecules_paper_counts molecules_papers.txt \\\n"
        f"    --lotus simulated_lotus.tsv \\\n"
        f"    --simple_data simulated_simple_data.tsv \\\n"
        f"    --epsilon_simple_model {config.epsilon} \\\n"
        f"    --gamma {config.gamma} \\\n"
        f"    --epsilon {config.error_rate} \\\n"
        f"    --error_probability {config.error_probability} \\\n"
        f"    --iterations {config.iterations} \\\n"
        f"    --burnin {config.burnin} \\\n"
        f"    --numBurnin {config.n_burnin_rounds} \\\n"
        f"    --fixedSeed {config.seed + seed_offset} \\\n"
        f"    --numThreads {threads} \\\n"
        f"    --writeBurnin \\\n"
        f"    --write_joint_log_prob_density"
        + (f" \\\n    {extra}" if extra else "")
        + "\n"
    )


PIN_FIELD = "--set_Y simulated_Y.txt --Y.update false"
PIN_STATES = (
    "--set_species_Z simulated_Z_species.txt "
    "--set_molecules_Z simulated_Z_molecules.txt --Z.update false"
)

RUNGS = [
    ("rung1_pin_field_and_states", "s", f"{PIN_FIELD} {PIN_STATES}", 0),
    ("rung2_pin_field", "s", PIN_FIELD, 1),
    ("rung3_from_data_s", "s", "", 2),
    ("rung3_from_data_l", "l", "", 3),
    ("rung3_from_data_ls", "ls", "", 4),
]


def _write_scripts(out: pathlib.Path, config: ScenarioConfig) -> None:
    for name, flags, extra, offset in RUNGS:
        path = out / f"{name}.sh"
        path.write_text(_infer_script(name, flags, extra, config, offset))
        path.chmod(0o755)

    # The C++ simulator, run under the very parameters the reference drew its own
    # scenario from. `compare_fields.py` then holds the two side by side. Every
    # parameter of *both* trees is pinned, because both are active now.
    replicates = out / "replicates.sh"
    replicates.write_text(
        _PREAMBLE.format(flags="s") + f'N="${{1:-{DEFAULT_REPLICATES}}}"\n'
        'for i in $(seq 1 "$N"); do\n'
        '  mkdir -p "replicates/cpp_$i"\n'
        '  "$ACOL" simulate \\\n'
        '      --out "replicates/cpp_$i/acol" \\\n'
        "      --tree_species species.txt --tree_molecules molecules.txt \\\n"
        "      --species_paper_counts species_papers.txt \\\n"
        "      --molecules_paper_counts molecules_papers.txt \\\n"
        f"      --error_probability {config.error_probability} \\\n"
        "      --epsilon_simple_model 0.5 \\\n"
        "      --numThreads all --write_Y --write_Z \\\n"
        '      --fixedSeed "$((1000 + i))" \\\n'
        "      "
        + " \\\n      ".join(
            f"--{tree}_{p} {SIMULATE_PARAMETERS}"
            for tree in ("species", "molecules")
            for p in ("branch_lengths", "mean_log_nu", "var_log_nu", "log_nu", "alpha")
        )
        + "\n"
        "done\n"
        "\n"
        'echo "Now: uv run python compare_fields.py . --replicates $N"\n'
    )
    replicates.chmod(0o755)
