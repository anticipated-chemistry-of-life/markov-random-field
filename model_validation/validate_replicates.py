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

The difference report is the second mode. A scenario that is misspecified on purpose has
no truth to gate against, so it gets a size of bias. The mode takes a scenario family
and a baseline family with the same seeds and layout. For each score it prints the mean
over the seeds, the baseline's mean and the difference. It has no gate. `--baseline` is
a glob for the baseline family's scenario directories. `--baseline-rung` is the
baseline's rung, and it defaults to RUNG. Rung 2 against rung 3 of one family:

    uv run python validate_replicates.py rung2_pin_field \\
        independent_y_s255_m255_seed{1..16}_ms_dense \\
        --baseline-rung rung3_from_data_lm
"""

from __future__ import annotations

import glob
import json
import math
import pathlib
from dataclasses import dataclass

import click
import numpy as np

import validate_independent
from src.independent import recovery

TREES = validate_independent.TREES
COVERAGES = ("alpha", "log_nu")
MIN_SCENARIOS = 3
"""Two scenarios cannot estimate a spread."""

SCORES = {
    "assignment: mean posterior on the true assignment": (
        "mass_spec",
        "assignment",
        "mean_on_true_assignment",
    ),
    "assignment: P(unknown), unknown features": (
        "mass_spec",
        "assignment",
        "mean_unknown_of_unknown_features",
    ),
    "field false presence": ("mass_spec", "field_false_presence", "mean"),
    "AUC, species with runs": ("mass_spec", "field", "with_runs", "auc"),
    "AUC, species without runs": ("mass_spec", "field", "without_runs", "auc"),
    **{
        f"{tree} {key.replace('_', ' ')} coverage": (tree, key, "coverage")
        for tree in TREES
        for key in COVERAGES
    },
    "kappa coverage": ("mass_spec", "kappa", "covered"),
}
"""Each score of the difference report by label, with its path in a validator summary.

The filter-probability coverages follow, one per stratum of runs present."""


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


def _entry_at(summary: dict, path: tuple[str, ...]) -> object:
    """The entry at `path` in a summary, or None where the summary has no such entry."""
    node: object = summary
    for key in path:
        if not isinstance(node, dict) or key not in node:
            return None
        node = node[key]
    return node


def _number_at(summary: dict, path: tuple[str, ...]) -> float:
    """The number at `path` in a summary, or NaN where the summary has no such entry."""
    entry = _entry_at(summary, path)
    if entry is None:
        return math.nan
    if not isinstance(entry, (int, float)):
        raise TypeError(f"The summary holds {entry!r} at {path}, not a number.")
    return float(entry)


def _stratum_coverage(summary: dict, stratum: str) -> float:
    rows = _entry_at(summary, ("mass_spec", "filter_probabilities", "strata"))
    if not isinstance(rows, list):
        return math.nan
    return next((row["coverage"] for row in rows if row["label"] == stratum), math.nan)


def scenario_scores(summary: dict) -> dict[str, float]:
    """One scenario's scores by label, in the order of the report. A score it lacks is NaN."""
    scores = {label: _number_at(summary, path) for label, path in SCORES.items()}
    for stratum, _, _ in recovery.RUN_STRATA:
        scores[f"filter probability coverage, runs {stratum}"] = _stratum_coverage(
            summary, stratum
        )
    return scores


def difference_scores(scenarios: list[dict], baselines: list[dict]) -> dict[str, dict]:
    """Each score's mean over the scenarios, its mean over the baselines, and the difference.

    `scenarios[i]` and `baselines[i]` are summaries of the same seed. A seed with no value
    for a score, in either family, leaves that seed out of both means, so the two means
    average the same seeds.
    """
    scenario_values = [scenario_scores(summary) for summary in scenarios]
    baseline_values = [scenario_scores(summary) for summary in baselines]
    rows = {}
    for label in scenario_values[0]:
        pairs = np.array(
            [(a[label], b[label]) for a, b in zip(scenario_values, baseline_values)]
        )
        pairs = pairs[~np.isnan(pairs).any(axis=1)]
        if not len(pairs):
            continue
        mean, baseline = pairs.mean(axis=0)
        rows[label] = {
            "seeds": len(pairs),
            "mean": float(mean),
            "baseline": float(baseline),
            "difference": float(mean - baseline),
        }
    return rows


def _report_difference(rung: str, baseline_rung: str, rows: dict[str, dict]) -> None:
    click.echo(f"\n=== {rung} against {baseline_rung}: mean over scenarios ===")
    width = max((len(label) for label in rows), default=0)
    click.echo(
        f"  {'score':<{width}} {'mean':>8} {'baseline':>9} {'difference':>11} {'seeds':>6}"
    )
    for label, row in rows.items():
        click.echo(
            f"  {label:<{width}} {row['mean']:8.4f} {row['baseline']:9.4f} "
            f"{row['difference']:+11.4f} {row['seeds']:6d}"
        )


def _echo_warnings(summaries: list[dict]) -> None:
    for summary in summaries:
        for warning in summary["warnings"]:
            click.echo(f"  ! {summary['scenario']} {summary['rung']}: {warning}")


@dataclass(frozen=True)
class Scenario:
    """A scenario directory with the two facts of its `meta.json` that pair it with a baseline."""

    path: pathlib.Path
    seed: int
    layout: str | None
    """The mass-spec layout, or None for a scenario without mass spec."""

    @classmethod
    def read(cls, path: pathlib.Path) -> Scenario:
        meta = json.loads((path / "meta.json").read_text())
        return cls(path, meta["seed"], meta.get("mass_spec", {}).get("layout"))


def _find_baselines(pattern: str, scenarios: list[Scenario]) -> list[pathlib.Path]:
    """The baseline scenario of each scenario's seed, found by glob.

    A match that holds no scenario, or whose seed is no scenario's seed, is left out.
    Two matches of one seed are an error: the pattern then reaches past one family.
    """
    found: dict[int, Scenario] = {}
    for path in sorted(map(pathlib.Path, glob.glob(pattern))):
        if not (path / "meta.json").is_file():
            continue
        match = Scenario.read(path)
        if match.seed in found:
            raise click.ClickException(
                f"Baseline scenarios {found[match.seed].path} and {path} both have seed "
                f"{match.seed}. The pattern matches more than one family."
            )
        found[match.seed] = match
    baselines = []
    for scenario in scenarios:
        baseline = found.get(scenario.seed)
        if baseline is None:
            raise click.ClickException(
                f"No baseline scenario has seed {scenario.seed}, "
                f"the seed of {scenario.path}."
            )
        if baseline.layout != scenario.layout:
            raise click.ClickException(
                f"Seed {scenario.seed} has layout {scenario.layout} in {scenario.path} "
                f"and {baseline.layout} in {baseline.path}. The baseline needs the "
                "layout of its scenario."
            )
        baselines.append(baseline.path)
    return baselines


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
@click.option(
    "--baseline",
    default=None,
    help="The baseline family of the difference report: a glob for its scenario "
    "directories. Quote it, so that the tool expands it. Each scenario's seed picks "
    "its match.",
)
@click.option(
    "--baseline-rung",
    default=None,
    help="The baseline's rung. It defaults to RUNG. Without --baseline, the baseline "
    "is this rung of the same scenarios.",
)
def main(
    rung: str,
    scenario_dirs: tuple[str, ...],
    gate: float | None,
    baseline: str | None,
    baseline_rung: str | None,
) -> None:
    """Score RUNG in every SCENARIO_DIR and pool the coverages.

    With --baseline or --baseline-rung, print the difference report instead: each score's
    mean over the seeds, the baseline's mean and the difference, with no gate.
    """
    paths = [pathlib.Path(d) for d in scenario_dirs]
    if baseline is not None or baseline_rung is not None:
        if gate is not None:
            raise click.UsageError("The difference report has no gate. Drop --gate.")
        baseline_rung = baseline_rung or rung
        baseline_paths = (
            _find_baselines(baseline, [Scenario.read(path) for path in paths])
            if baseline is not None
            else paths
        )
        summaries = [validate_independent.score(path, rung) for path in paths]
        baselines = [
            validate_independent.score(path, baseline_rung) for path in baseline_paths
        ]
        _echo_warnings([*summaries, *baselines])
        _report_difference(rung, baseline_rung, difference_scores(summaries, baselines))
        return

    if len(scenario_dirs) < MIN_SCENARIOS:
        raise click.ClickException(
            f"Need at least {MIN_SCENARIOS} scenarios to estimate a spread, "
            f"got {len(scenario_dirs)}."
        )
    summaries = [validate_independent.score(path, rung) for path in paths]
    _echo_warnings(summaries)

    scores = coverage_scores(summaries)
    _report(rung, scores)

    if gate is not None and not apply_gate(scores, gate):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
