"""Generate a validation scenario whose field owes nothing to the C++ binary.

The field is sampled by one pass down the species tree — root from the stationary
distribution, then every node given its parent — with no molecules tree involved
at all. The molecules dimension is then pinned neutral so the C++ model reduces
exactly to that process (see ADR-0001).

    uv run python simulate_independent.py --seed 42
    bash <scenario>/rung1_pin_field_and_states.sh
    uv run python validate_independent.py <scenario> rung1_pin_field_and_states
"""

import pathlib

import click

from src.independent.scenario import RUNGS, ScenarioConfig, build_scenario


@click.command()
@click.option("--seed", type=int, default=42, show_default=True)
@click.option("--number_of_species", type=int, default=255, show_default=True,
              help="Nodes in the species tree (2^k - 1).")
@click.option("--number_of_molecules", type=int, default=255, show_default=True,
              help="Nodes in the molecules tree (2^k - 1).")
@click.option("--mean_log_nu", type=float, default=-0.5, show_default=True)
@click.option("--var_log_nu", type=float, default=0.25, show_default=True,
              help="Variance, not standard deviation.")
@click.option("--epsilon", type=float, default=0.05, show_default=True,
              help="Simple error model per-cell flip probability.")
@click.option("--gamma", type=float, default=1.1, show_default=True)
@click.option("--error_rate", type=float, default=0.001, show_default=True,
              help="LOTUS probability of a record for an absent pair.")
@click.option("--iterations", type=int, default=10_000, show_default=True)
@click.option("--burnin", type=int, default=1_000, show_default=True,
              help="Iterations per burn-in round.")
@click.option("--n_burnin_rounds", type=int, default=10, show_default=True)
@click.option("--true_branch_lengths", is_flag=True,
              help="Start the chain on the true branch lengths instead of flat ones.")
@click.option("--out", type=click.Path(), default=None,
              help="Output directory (default: independent_y_s<n>_m<n>_seed<seed>).")
def main(out: str | None, **kwargs) -> None:
    """Write a complete scenario: truth, observations, and per-rung run scripts."""
    config = ScenarioConfig(
        seed=kwargs["seed"],
        n_species_nodes=kwargs["number_of_species"],
        n_molecule_nodes=kwargs["number_of_molecules"],
        mean_log_nu=kwargs["mean_log_nu"],
        var_log_nu=kwargs["var_log_nu"],
        epsilon=kwargs["epsilon"],
        gamma=kwargs["gamma"],
        error_rate=kwargs["error_rate"],
        iterations=kwargs["iterations"],
        burnin=kwargs["burnin"],
        n_burnin_rounds=kwargs["n_burnin_rounds"],
        true_branch_lengths=kwargs["true_branch_lengths"],
    )

    if out is None:
        out = (
            f"independent_y_s{config.n_species_nodes}"
            f"_m{config.n_molecule_nodes}_seed{config.seed}"
        )
    path = pathlib.Path(out)
    meta = build_scenario(path, config)

    click.echo(f"Scenario written to {path.resolve()}")
    click.echo(
        f"  {meta['n_species_leaves']} x {meta['n_molecule_leaves']} field, "
        f"{meta['field_ones_fraction']:.1%} present"
    )
    click.echo(
        f"  {meta['n_cliques']} cliques, {meta['n_species_branches']} branches, "
        f"branch-length budget {sum(meta['species_bins'])}"
    )
    click.echo(
        f"  {meta['lotus_records']} LOTUS records, "
        f"{meta['simple_data_disagreements']} simple-error disagreements"
    )
    if not meta["true_branch_lengths"]:
        click.echo("  branch lengths start flat: the chain has to find them")
    click.echo("\nRun in order, stopping at the first failure:")
    click.echo(f"  bash {path}/check_neutrality_invariant.sh")
    for name, _, _, _ in RUNGS:
        click.echo(f"  bash {path}/{name}.sh && "
                   f"uv run python validate_independent.py {path} {name}")


if __name__ == "__main__":
    main()
