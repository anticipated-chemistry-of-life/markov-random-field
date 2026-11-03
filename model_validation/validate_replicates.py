"""Score one rung over many independently simulated scenarios.

`validate_independent.py` gates one scenario's coverage against a binomial band,
which treats every clique as an independent trial. Once the node states are
inferred, log nu's cliques are not: one shared var log nu sets every clique's
interval at once. Across 32 seeds the coverage of an exact posterior with rung 2's
information scatters at an sd of 0.052-0.065, against the band's 0.027, so one
scenario says little about the sampler (#112). The mean over many scenarios does.
Its standard error comes from the scenarios' own scatter, as in `compare_fields.py`,
and the gate scores the mean's distance from the nominal coverage in those units.

Nominal coverage holds on average over the prior. A scenario fixes mean and var
log nu rather than drawing them, so an exact posterior's mean coverage need not
be exactly 0.9: on those 32 seeds it measured 0.900-0.917. That offset is small
against the standard error of a few dozen scenarios, but not against hundreds.

Alpha and log nu are scored alike, at every rung. Run from `model_validation/`:

    for s in $(seq 1 16); do uv run python simulate_independent.py --seed "$s"; done
    pixi run build release-s
    printf '%s\\n' independent_y_s255_m255_seed{1..16} |
        xargs -P 16 -I{} env ACOL_NO_BUILD=1 ACOL_THREADS=4 bash {}/rung2_pin_field.sh
    uv run python validate_replicates.py rung2_pin_field \\
        independent_y_s255_m255_seed{1..16} --gate 3
"""

from __future__ import annotations

import math
import pathlib

import click
import numpy as np

import validate_independent
from src.independent import recovery

TREES = validate_independent.TREES
COVERAGES = ("alpha", "log_nu")
MIN_SCENARIOS = 3
"""Two scenarios cannot estimate a spread."""


def coverage_scores(summaries: list[dict]) -> dict:
    """Each tree's mean coverage over the scenarios, in its own standard errors."""
    scores = {}
    for tree in TREES:
        for key in COVERAGES:
            blocks = [summary[tree][key] for summary in summaries]
            if not all(block["n"] for block in blocks):
                continue
            per_scenario = np.array([block["coverage"] for block in blocks])
            mean = float(per_scenario.mean())
            sd = float(per_scenario.std(ddof=1))
            standard_error = sd / math.sqrt(len(per_scenario))
            deviation = abs(mean - recovery.NOMINAL_COVERAGE)
            if standard_error > 0.0:
                sigmas = deviation / standard_error
            else:
                # Every scenario agrees exactly: right, or wrong by a distance no
                # scatter can explain.
                sigmas = 0.0 if deviation == 0.0 else math.inf
            scores[f"{tree} {key}"] = {
                "scenarios": len(per_scenario),
                "mean": mean,
                "sd": sd,
                "standard_error": standard_error,
                "sigmas": sigmas,
                "outside_band": sum(not block["passes"] for block in blocks),
            }
    return scores


def _report(rung: str, scores: dict) -> None:
    click.echo(f"\n=== {rung}: mean coverage over scenarios ===")
    for name, block in scores.items():
        click.echo(
            f"  {name:<18} mean {block['mean']:6.1%}  sd {block['sd']:.3f}  "
            f"se {block['standard_error']:.3f}  {block['sigmas']:5.2f} se from "
            f"{recovery.NOMINAL_COVERAGE:.0%}  ({block['outside_band']} of "
            f"{block['scenarios']} outside the one-scenario band)"
        )


def apply_gate(scores: dict, gate: float) -> bool:
    """Every mean coverage within `gate` standard errors of nominal."""
    failures = [
        f"{name} mean coverage {block['mean']:.1%} is {block['sigmas']:.2f} standard "
        f"errors from {recovery.NOMINAL_COVERAGE:.0%}"
        for name, block in scores.items()
        if block["sigmas"] > gate
    ]

    click.echo(f"\n  gate at {gate} standard errors:")
    for failure in failures:
        click.echo(f"    FAIL {failure}")
    if not failures:
        click.echo("    PASS")
    return not failures


@click.command()
@click.argument("rung")
@click.argument(
    "scenario_dirs", nargs=-1, required=True, type=click.Path(exists=True, file_okay=False)
)
@click.option(
    "--gate",
    type=float,
    default=None,
    help="Exit 1 when a mean coverage is more than this many standard errors from nominal.",
)
def main(rung: str, scenario_dirs: tuple[str, ...], gate: float | None) -> None:
    """Score RUNG in every SCENARIO_DIR and pool the coverages."""
    if len(scenario_dirs) < MIN_SCENARIOS:
        raise click.ClickException(
            f"Need at least {MIN_SCENARIOS} scenarios to estimate a spread, "
            f"got {len(scenario_dirs)}."
        )
    summaries = [validate_independent.score(pathlib.Path(d), rung) for d in scenario_dirs]
    for summary in summaries:
        for warning in summary["warnings"]:
            click.echo(f"  ! {summary['scenario']}: {warning}")

    scores = coverage_scores(summaries)
    _report(rung, scores)

    if gate is not None and not apply_gate(scores, gate):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
