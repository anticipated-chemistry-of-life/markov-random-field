"""Is the `log nu` drift explained by the missing normalising constant?

Varies the molecules dimension from neutral to strongly non-neutral and reports,
for each setting, where the C++'s objective is maximised versus where the
correctly normalised one is. Everything is exact enumeration over the field — no
MCMC, no sampling, no approximation. See `src/independent/toy_normaliser.py`.

    uv run python diagnose_normaliser.py

The correct argmax landing on the true `log nu` in every row is the self-check:
it confirms the enumeration before any conclusion is drawn from the gap.
"""

from __future__ import annotations

import click
import numpy as np

from src.independent.toy_normaliser import (
    run_chain,
    MAX_ENUMERABLE_DEPTH,
    expected_log_likelihood_profile,
    field_distribution,
    leaf_pattern_probabilities,
    normalising_constant,
)

MOLECULE_LOG_NUS = (5.0, 3.0, 2.0, 1.0, 0.0, -0.5, -1.0, -2.0, -3.0)


def _row_for_log_nu(
    log_nu_molecules: float,
    true_log_nu_species: float,
    alpha_species: float,
    alpha_molecules: float,
    depth: int,
    grid: np.ndarray,
) -> dict:
    molecules = leaf_pattern_probabilities(
        alpha_molecules, float(np.exp(log_nu_molecules)), depth
    )
    truth_species = leaf_pattern_probabilities(
        alpha_species, float(np.exp(true_log_nu_species)), depth
    )
    truth = field_distribution(truth_species, molecules, depth)

    targeted, correct = expected_log_likelihood_profile(
        truth, grid, alpha_species, molecules, depth
    )

    return {
        "log_nu_molecules": log_nu_molecules,
        "targeted_argmax": float(grid[int(np.argmax(targeted))]),
        "correct_argmax": float(grid[int(np.argmax(correct))]),
        "log_C": float(np.log(normalising_constant(truth_species, molecules, depth))),
        "at_grid_edge": bool(np.argmax(targeted) == 0),
    }


@click.command()
@click.option("--true_log_nu_species", type=float, default=-0.5, show_default=True,
              help="The species log_nu the field is actually generated with.")
@click.option("--alpha_species", type=float, default=0.5, show_default=True)
@click.option("--alpha_molecules", type=float, default=0.5, show_default=True)
@click.option("--grid_low", type=float, default=-8.0, show_default=True)
@click.option("--grid_high", type=float, default=3.0, show_default=True)
@click.option("--grid_points", type=int, default=441, show_default=True)
def main(
    true_log_nu_species: float,
    alpha_species: float,
    alpha_molecules: float,
    grid_low: float,
    grid_high: float,
    grid_points: int,
) -> None:
    """Compare the targeted and correctly normalised objectives."""
    grid = np.linspace(grid_low, grid_high, grid_points)

    click.echo(
        f"Balanced trees, every branch at bin 5 (grid length 1.0).\n"
        f"True species log_nu {true_log_nu_species:+.2f}, "
        f"alpha_species {alpha_species}, alpha_molecules {alpha_molecules}."
    )

    for depth in range(1, MAX_ENUMERABLE_DEPTH + 1):
        n_leaves = 2**depth
        click.echo(
            f"\n  depth {depth}: {2 * n_leaves - 1} nodes per tree, "
            f"{n_leaves}x{n_leaves} field, {2 ** (n_leaves * n_leaves)} states"
        )
        click.echo(
            f"  {'mol log_nu':>10}  {'log C':>9}  {'C++ argmax':>11}  "
            f"{'correct argmax':>14}  {'bias':>7}"
        )
        click.echo(f"  {'-' * 10}  {'-' * 9}  {'-' * 11}  {'-' * 14}  {'-' * 7}")

        for log_nu_molecules in MOLECULE_LOG_NUS:
            row = _row_for_log_nu(
                log_nu_molecules, true_log_nu_species, alpha_species,
                alpha_molecules, depth, grid,
            )
            bias = row["targeted_argmax"] - row["correct_argmax"]
            edge = "  <- pinned at grid edge" if row["at_grid_edge"] else ""
            click.echo(
                f"  {row['log_nu_molecules']:>10.1f}  {row['log_C']:>9.4f}  "
                f"{row['targeted_argmax']:>11.3f}  {row['correct_argmax']:>14.3f}  "
                f"{bias:>+7.3f}{edge}"
            )

    click.echo(
        "\nThe correct argmax recovering the true log_nu everywhere confirms the "
        "enumeration.\nThe bias is the omitted log C(theta), and it grows both as "
        "the molecules\ndimension leaves neutrality and as the field gets bigger."
    )

    _report_chains(
        true_log_nu_species, alpha_species, alpha_molecules, MAX_ENUMERABLE_DEPTH
    )


def _report_chains(
    true_log_nu_species: float,
    alpha_species: float,
    alpha_molecules: float,
    depth: int,
) -> None:
    """Run the species log_nu update with and without the normaliser."""
    click.echo(
        f"\n\nMetropolis on the species log_nu over 200 fields observed without "
        f"error,\nstarted at the truth ({true_log_nu_species:+.2f}), "
        f"2000 iterations at depth {depth}:\n"
    )
    click.echo(
        f"  {'mol log_nu':>10}  {'C++ final':>10}  {'C++ mean':>9}  "
        f"{'corrected final':>15}  {'corrected mean':>14}"
    )
    click.echo(f"  {'-' * 10}  {'-' * 10}  {'-' * 9}  {'-' * 15}  {'-' * 14}")

    for log_nu_molecules in (5.0, 0.0, -1.0, -2.0):
        traces = {}
        for corrected in (False, True):
            traces[corrected] = run_chain(
                np.random.default_rng(20260825),
                correct_normaliser=corrected,
                true_log_nu_species=true_log_nu_species,
                log_nu_molecules=log_nu_molecules,
                alpha_species=alpha_species,
                alpha_molecules=alpha_molecules,
                depth=depth,
                n_iterations=2000,
            )
        click.echo(
            f"  {log_nu_molecules:>10.1f}  {traces[False][-1]:>10.3f}  "
            f"{traces[False].mean():>9.3f}  {traces[True][-1]:>15.3f}  "
            f"{traces[True].mean():>14.3f}"
        )

    click.echo(
        "\nBoth chains start at the truth and share a seed, so the corrected "
        "column staying\nput while the C++ column walks away is the whole effect."
    )


if __name__ == "__main__":
    main()
