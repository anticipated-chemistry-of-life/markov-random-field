"""Generate a validation scenario whose field owes nothing to the C++ binary.

Both trees are active. Each is sampled by one pass down its own nodes — root from
the stationary distribution, then every node given its parent, leaves included —
and the field is a noisy AND of the two leaf blocks. Nothing is neutralised: the
reference is exact with both trees running (see ADR-0005).

    uv run python simulate_independent.py --seed 42
    bash <scenario>/rung1_pin_field_and_states.sh
    uv run python validate_independent.py <scenario> rung1_pin_field_and_states
"""

import pathlib

import click

from src.independent.mass_spec import LAYOUTS, MassSpecOptions
from src.independent.scenario import (
    OUTSIDE_CONTRAST_NAME,
    OUTSIDE_SUFFIX,
    ScenarioConfig,
    build_scenario,
    rungs_of,
)


@click.command()
@click.option("--seed", type=int, default=42, show_default=True)
@click.option(
    "--number_of_species",
    type=int,
    default=255,
    show_default=True,
    help="Nodes in the species tree (2^k - 1).",
)
@click.option(
    "--number_of_molecules",
    type=int,
    default=255,
    show_default=True,
    help="Nodes in the molecules tree (2^k - 1).",
)
@click.option("--mean_log_nu", type=float, default=-0.5, show_default=True)
@click.option(
    "--var_log_nu",
    type=float,
    default=0.25,
    show_default=True,
    help="Variance, not standard deviation.",
)
@click.option(
    "--error_probability",
    type=float,
    default=0.05,
    show_default=True,
    help="Omega: the rate a tree field cell is corrupted before the AND, in (0, 0.5).",
)
@click.option(
    "--epsilon",
    type=float,
    default=0.05,
    show_default=True,
    help="Simple error model per-cell flip probability.",
)
@click.option("--gamma", type=float, default=1.1, show_default=True)
@click.option(
    "--error_rate",
    type=float,
    default=0.001,
    show_default=True,
    help="LOTUS probability of a record for an absent pair.",
)
@click.option("--iterations", type=int, default=10_000, show_default=True)
@click.option(
    "--burnin",
    type=int,
    default=1_000,
    show_default=True,
    help="Iterations per burn-in round.",
)
@click.option("--n_burnin_rounds", type=int, default=10, show_default=True)
@click.option(
    "--true_branch_lengths",
    is_flag=True,
    help="Start the chain on the true branch lengths instead of flat ones.",
)
@click.option(
    "--mass_spec",
    "layout",
    type=click.Choice(sorted(LAYOUTS)),
    default=None,
    help="Add simulated mass-spec data (LOTUS+MS) in this layout: dense is two runs "
    "per species leaf, sparse is eight species with one run each. The rungs then "
    "run on the 'lm' build.",
)
@click.option(
    "--outside_pool",
    type=click.IntRange(min=0),
    default=0,
    show_default=True,
    help="With --mass_spec: the number of outside molecules (ADR-0020), which have no "
    "leaf in the tree. Each is detected per run at the p_o that makes outside features "
    "twice the molecule features, and decoy rows draw from the leaves and the pool "
    "together. The reported scenario uses 4 times the molecule leaves. 0 is none, and "
    "leaves the scenario as it was.",
)
@click.option(
    "--out",
    type=click.Path(),
    default=None,
    help="Output directory (default: independent_y_s<n>_m<n>_seed<seed>).",
)
def main(out: str | None, layout: str | None, **kwargs) -> None:
    """Write a complete scenario: truth, observations, and per-rung run scripts."""
    if kwargs["outside_pool"] and not layout:
        raise click.UsageError("--outside_pool needs --mass_spec.")
    config = ScenarioConfig(
        seed=kwargs["seed"],
        n_species_nodes=kwargs["number_of_species"],
        n_molecule_nodes=kwargs["number_of_molecules"],
        mean_log_nu=kwargs["mean_log_nu"],
        var_log_nu=kwargs["var_log_nu"],
        error_probability=kwargs["error_probability"],
        epsilon=kwargs["epsilon"],
        gamma=kwargs["gamma"],
        error_rate=kwargs["error_rate"],
        iterations=kwargs["iterations"],
        burnin=kwargs["burnin"],
        n_burnin_rounds=kwargs["n_burnin_rounds"],
        true_branch_lengths=kwargs["true_branch_lengths"],
        mass_spec=LAYOUTS[layout] if layout else None,
        ms_options=MassSpecOptions(outside_pool=kwargs["outside_pool"]),
    )

    if out is None:
        out = (
            f"independent_y_s{config.n_species_nodes}"
            f"_m{config.n_molecule_nodes}_seed{config.seed}"
        )
        if layout:
            out += f"_ms_{layout}"
        if config.ms_options.outside_pool:
            out += OUTSIDE_SUFFIX
    path = pathlib.Path(out)
    meta = build_scenario(path, config)

    click.echo(f"Scenario written to {path.resolve()}")
    click.echo(
        f"  {meta['n_species_leaves']} x {meta['n_molecule_leaves']} field, "
        f"{meta['field_ones_fraction']:.1%} present "
        f"(the two adjusted rates predict {meta['expected_field_ones_fraction']:.1%})"
    )
    click.echo(
        f"  tree fields {meta['species_tree_field_ones_fraction']:.1%} and "
        f"{meta['molecule_tree_field_ones_fraction']:.1%} present, "
        f"omega {meta['error_probability']}"
    )
    click.echo(
        f"  {meta['n_species_branches']} species branches and "
        f"{meta['n_molecule_branches']} molecule branches, budgets "
        f"{sum(meta['species_bins'])} and {sum(meta['molecule_bins'])}"
    )
    click.echo(
        f"  {meta['lotus_records']} LOTUS records, "
        f"{meta['simple_data_disagreements']} simple-error disagreements"
    )
    if "mass_spec" in meta:
        ms = meta["mass_spec"]
        click.echo(
            f"  {ms['layout']} mass spec: {ms['n_runs']} runs of "
            f"{ms['n_species_with_runs']} species, {ms['n_features']} features "
            f"({ms['n_unknown_features']} unknown), lambda "
            f"{ms['expected_features_per_run']:.2f}, rho {ms['unknown_rate']:.4g}"
        )
        if "outside_pool" in ms:
            click.echo(
                f"  outside pool {ms['outside_pool']}: {ms['n_outside_features']} outside "
                f"features, p_o {ms['outside_detection']:.4g}, true tau {ms['outside_rate']:.4g}"
            )
    if not meta["true_branch_lengths"]:
        click.echo("  branch lengths start flat: the chain has to find them")
    click.echo("\nRun in order, stopping at the first failure:")
    for name, _, _, _ in rungs_of(config):
        if name == OUTSIDE_CONTRAST_NAME:
            click.echo("\nThe contrast, rung 3 at tau = 0, is a separate run:")
        click.echo(
            f"  bash {path}/{name}.sh && "
            f"uv run python validate_independent.py {path} {name}"
        )


if __name__ == "__main__":
    main()
