"""Assemble a complete validation scenario: truth, observations, and run scripts.

The three rungs are run in the order written and abandoned at the first failure.
Each pins strictly less than the one before, so a failure localises the fault:

1. field *and* internal states pinned — only the parameter updates move, so this
   is close to closed form and sets the empirical ceiling for the other two.
2. field pinned, internal states inferred — adds the Z Gibbs sweep.
3. nothing pinned, inferred from observations — the production path, run against
   the simple error model alone, then LOTUS alone, then both.

Throughout, the molecules dimension is pinned neutral (ADR-0001).
"""

from __future__ import annotations

import json
import pathlib
from dataclasses import dataclass

import numpy as np

from ..tree import Tree, TreeType
from . import field as F
from . import io
from .data import sample_paper_counts, simulate_lotus, simulate_simple_error
from .indexing import TreeIndex, build_tree_index

# Puts nu ~ 148, comfortably past the C++'s stationary threshold, so the
# molecules transition rows are exactly (0.5, 0.5).
NEUTRAL_LOG_NU = 5.0
NEUTRAL_ALPHA = 0.5

# Gibbs sweeps the C++ runs after its own one-pass initialisation. The reference
# does a single pass and stops; if the C++'s sweeps target the same distribution,
# the two agree, and if they do not, that is the bug being hunted.
SIMULATION_SWEEPS = 1000
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
PINNED_MOLECULES = "simulated_pinned_molecules.txt"
SIMULATE_PARAMETERS = "simulated_parameters.txt"


@dataclass
class ScenarioConfig:
    seed: int = 42
    n_species_nodes: int = 255
    n_molecule_nodes: int = 255
    mean_log_nu: float = -0.5
    var_log_nu: float = 0.25
    epsilon: float = 0.05
    gamma: float = 1.1
    error_rate: float = 0.001
    iterations: int = 10_000
    burnin: int = 1_000
    n_burnin_rounds: int = 10
    true_branch_lengths: bool = False
    """Write grid centres into the species tree file instead of flat lengths.

    Flat is the default: it makes the chain find the branch lengths rather than
    start on them, so 'never left the initialisation' cannot pass as recovery.
    """


def _edges(tree: Tree) -> list[tuple[str, str]]:
    frame = tree.to_dataframe()
    return list(zip(frame["child"].astype(str), frame["parent"].astype(str)))


def build_scenario(out: pathlib.Path, config: ScenarioConfig) -> dict:
    rng = np.random.default_rng(config.seed)
    out.mkdir(parents=True, exist_ok=True)

    species_tree = Tree(config.n_species_nodes, TreeType.balanced, "species")
    molecule_tree = Tree(config.n_molecule_nodes, TreeType.balanced, "molecules")
    species_edges, molecule_edges = _edges(species_tree), _edges(molecule_tree)
    species = build_tree_index(species_edges)
    molecules = build_tree_index(molecule_edges)

    # -- species truth: one clique per molecule leaf --------------------------
    n_cliques = molecules.n_leaves
    alphas = rng.beta(0.5, 0.5, size=n_cliques)
    log_nus = rng.normal(config.mean_log_nu, np.sqrt(config.var_log_nu), size=n_cliques)
    nus = np.exp(log_nus)
    if nus.max() >= F.STATIONARY_NU_THRESHOLD:
        raise ValueError(
            f"A species clique drew nu = {nus.max():.2f}, at or past the C++'s "
            f"stationary short-circuit ({F.STATIONARY_NU_THRESHOLD}). The reference "
            "does not replicate that approximation, so the two would diverge "
            "legitimately. Re-seed or lower var_log_nu."
        )

    species_bins = F.sample_binned_branch_lengths(rng, species.n_branches)

    states = F.sample_states(rng, species, species_bins, alphas, nus)
    latent_field = states[species.leaves]
    species_internal_states = states[species.internals]

    # -- molecules: neutral, so its internal states are arbitrary -------------
    neutral_nu = float(np.exp(NEUTRAL_LOG_NU))
    if neutral_nu <= F.STATIONARY_NU_THRESHOLD:
        raise ValueError(
            f"Neutral nu = {neutral_nu:.2f} does not exceed the stationary "
            f"threshold ({F.STATIONARY_NU_THRESHOLD}); the molecules dimension "
            "would not be neutral and the whole scenario would be invalid."
        )
    molecule_internal_states = (
        rng.random((species.n_leaves, molecules.n_internals)) < 0.5
    )
    # A second, unrelated draw. Under neutrality the species posteriors must not
    # be able to tell these two apart; `check_neutrality_invariant.sh` checks it.
    molecule_internal_states_alt = (
        rng.random((species.n_leaves, molecules.n_internals)) < 0.5
    )
    molecule_bins = np.full(molecules.n_branches, F.N_BINS // 2, dtype=np.int64)

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
    species_lengths = (
        F.grid_branch_lengths()[species_bins]
        if config.true_branch_lengths
        else np.full(species.n_branches, 0.2)
    )
    io.write_tree(out / "species.txt", species_edges, species, species_lengths)
    io.write_tree(
        out / "molecules.txt",
        molecule_edges,
        molecules,
        np.full(molecules.n_branches, 0.2),
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

    # Two separate files, deliberately: nothing that the C++ is pointed at may
    # contain species truth, or a stray flag would silently start the chain on
    # the answer it is meant to find.
    #
    # The `simulated_` prefix is load-bearing, not decorative. stattools picks a
    # reader by *filename*: a name/value file is only matched up by parameter
    # name when the filename contains one of trace/simulated/meanVar/
    # statePosteriors/posteriorMode, and otherwise has to be a bare one-column
    # file of exactly the right length (TReadInitialValues.h:133). Renaming these
    # files without keeping the marker makes the C++ reject them.
    io.write_parameters(out / PINNED_MOLECULES, _pinned_values(species, molecules))
    species_truth = _species_truth(
        species, molecules, alphas, log_nus, species_bins, config
    )
    io.write_parameters(out / "truth_species.txt", species_truth)

    # The one file that legitimately hands species truth to the C++: the
    # replicate comparison needs *both* implementations running the same
    # parameters. Never referenced by an infer script.
    io.write_parameters(
        out / SIMULATE_PARAMETERS,
        {**species_truth, **_pinned_values(species, molecules)},
    )

    io.write_field(out / "simulated_Y.txt", latent_field, species, molecules)
    io.write_internal_states(
        out / "simulated_Z_species.txt", species_internal_states, species, molecules, 0
    )
    io.write_internal_states(
        out / "simulated_Z_molecules.txt",
        molecule_internal_states,
        species,
        molecules,
        1,
    )
    io.write_internal_states(
        out / "simulated_Z_molecules_alt.txt",
        molecule_internal_states_alt,
        species,
        molecules,
        1,
    )
    io.write_observations(out / "simulated_lotus.tsv", lotus, species, molecules)
    io.write_observations(
        out / "simulated_simple_data.tsv", simple_data, species, molecules
    )

    meta = {
        "seed": config.seed,
        "n_bins": F.N_BINS,
        "n_species_leaves": species.n_leaves,
        "n_molecule_leaves": molecules.n_leaves,
        "n_species_branches": species.n_branches,
        "n_cliques": n_cliques,
        "mean_log_nu": config.mean_log_nu,
        "var_log_nu": config.var_log_nu,
        "epsilon": config.epsilon,
        "gamma": config.gamma,
        "error_rate": config.error_rate,
        "iterations": config.iterations,
        "burnin": config.burnin,
        "n_burnin_rounds": config.n_burnin_rounds,
        "true_branch_lengths": config.true_branch_lengths,
        "neutral_alpha": NEUTRAL_ALPHA,
        "neutral_log_nu": NEUTRAL_LOG_NU,
        "field_ones_fraction": float(latent_field.mean()),
        "lotus_records": int(lotus.sum()),
        "simple_data_disagreements": int((simple_data != latent_field).sum()),
        "species_bins": species_bins.tolist(),
        "molecule_bins": molecule_bins.tolist(),
    }
    (out / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")

    _write_scripts(out, config)
    return meta


def _pinned_values(species: TreeIndex, molecules: TreeIndex) -> dict[str, object]:
    """Everything the molecules dimension needs to be exactly uninformative."""
    values: dict[str, object] = {}
    for name in molecules.branch_names():
        values[f"molecules_branch_lengths_{name}"] = F.N_BINS // 2
    values["molecules_mean_log_nu"] = NEUTRAL_LOG_NU
    # Positive but negligible: every log_nu is pinned to the same value, so the
    # variance has no spread to estimate and only the prior keeps it off zero.
    values["molecules_var_log_nu"] = 1e-4
    for name in species.leaf_names():
        values[f"molecules_log_nu_{name}"] = NEUTRAL_LOG_NU
    for name in species.leaf_names():
        values[f"molecules_alpha_{name}"] = NEUTRAL_ALPHA
    return values


def _species_truth(
    species: TreeIndex,
    molecules: TreeIndex,
    alphas: np.ndarray,
    log_nus: np.ndarray,
    bins: np.ndarray,
    config: ScenarioConfig,
) -> dict[str, object]:
    """The answers, for the validator only. Never passed to the C++."""
    values: dict[str, object] = {}
    for name, value in zip(species.branch_names(), bins):
        values[f"species_branch_lengths_{name}"] = int(value)
    values["species_mean_log_nu"] = config.mean_log_nu
    values["species_var_log_nu"] = config.var_log_nu
    for name, value in zip(molecules.leaf_names(), log_nus):
        values[f"species_log_nu_{name}"] = float(value)
    for name, value in zip(molecules.leaf_names(), alphas):
        values[f"species_alpha_{name}"] = float(value)
    return values


_PIN_MOLECULES = " \\\n    ".join(
    f"--molecules_{p} {PINNED_MOLECULES} --molecules_{p}.update false"
    for p in ("branch_lengths", "mean_log_nu", "var_log_nu", "log_nu", "alpha")
)

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
        f"    --iterations {config.iterations} \\\n"
        f"    --burnin {config.burnin} \\\n"
        f"    --numBurnin {config.n_burnin_rounds} \\\n"
        f"    --fixedSeed {config.seed + seed_offset} \\\n"
        f"    --numThreads {threads} \\\n"
        f"    --writeBurnin \\\n"
        f"    --write_joint_log_prob_density \\\n"
        f"    {_PIN_MOLECULES}" + (f" \\\n    {extra}" if extra else "") + "\n"
    )


PIN_FIELD = "--set_Y simulated_Y.txt --Y.update false"
PIN_STATES = (
    "--set_species_Z simulated_Z_species.txt "
    "--set_molecules_Z simulated_Z_molecules.txt --Z.update false"
)
PIN_STATES_ALT = PIN_STATES.replace(
    "simulated_Z_molecules.txt", "simulated_Z_molecules_alt.txt"
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

    # The neutrality invariant: the same run twice, differing only in which
    # arbitrary draw of the molecules internal states is pinned. Self-contained
    # rather than diffed against rung 1, so it does not depend on rung 1 having
    # been run or on its thread count.
    #
    # Single-threaded on purpose. With --numThreads all the binary is not
    # reproducible even under --fixedSeed: two identical invocations produce
    # posterior means differing by ~0.2 posterior standard deviations. An exact
    # diff is only a valid test of neutrality when the run is deterministic.
    invariant = out / "check_neutrality_invariant.sh"
    invariant.write_text(
        _PREAMBLE.format(flags="s")
        + _infer_command(
            "neutrality_a", f"{PIN_FIELD} {PIN_STATES}", config, 0, threads="1"
        )
        + "\n"
        + _infer_command(
            "neutrality_b", f"{PIN_FIELD} {PIN_STATES_ALT}", config, 0, threads="1"
        )
        + "\n"
        "if diff -q neutrality_a/acol_meanVar.txt \\\n"
        "          neutrality_b/acol_meanVar.txt >/dev/null; then\n"
        "  echo 'PASS: the molecules dimension is neutral.'\n"
        "else\n"
        "  echo 'FAIL: species posteriors depend on arbitrary molecule states.'\n"
        "  echo '      Every result in this scenario is meaningless until fixed.'\n"
        "  diff neutrality_a/acol_meanVar.txt neutrality_b/acol_meanVar.txt | head -20\n"
        "  exit 1\n"
        "fi\n"
    )
    invariant.chmod(0o755)

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
        f"      --iterations {SIMULATION_SWEEPS} \\\n"
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
