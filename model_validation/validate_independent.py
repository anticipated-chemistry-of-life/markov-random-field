"""Score one rung's inference output against the independently simulated truth.

Reports, and does not by default gate on:

- **coverage** for `alpha` and `log_nu` — the fraction of the 128 cliques whose
  95% credible interval contains the truth. Coverage catches both bias and
  over-confident posteriors; RMSE alone catches only the first.
- **Spearman and mean absolute error** for branch lengths, which are discrete and
  budget-constrained, so a credible interval is a strange object for them.
- **movement**, because a chain that never left its initialisation can otherwise
  masquerade as recovery when the branch lengths start flat.

Thresholds are deliberately absent: rung 1 pins the field *and* the internal
states, so whatever it scores is the empirical ceiling, and the gates for the
looser rungs should be set relative to that rather than guessed in advance. Pass
`--gates <rung1-summary.json>` once you have it.
"""

from __future__ import annotations

import json
import pathlib

import click
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

CREDIBLE_MASS = 0.95


def _read_truth(path: pathlib.Path) -> dict[str, float]:
    frame = pd.read_csv(path, sep="\t")
    return dict(zip(frame["name"].astype(str), frame["value"].astype(float)))


def _read_trace(rung_dir: pathlib.Path, meta: dict) -> pd.DataFrame:
    """The main trace plus the species tree's own, joined column-wise.

    Branch lengths are updated by the tree rather than by a stattools updater, so
    they are written to `acol_species_trace.txt` and appear nowhere in
    `acol_trace.txt`. Both files share a thinning factor and row count.
    """
    frames = [pd.read_csv(rung_dir / "acol_trace.txt", sep="\t")]

    species_trace = rung_dir / "acol_species_trace.txt"
    if species_trace.exists():
        frames.append(pd.read_csv(species_trace, sep="\t"))

    if len({len(f) for f in frames}) != 1:
        raise click.ClickException(
            f"Trace files disagree on row count: {[len(f) for f in frames]}."
        )
    trace = pd.concat(frames, axis=1)

    # The trace includes the burn-in rounds, which would contaminate every
    # quantile. Trim by the configured proportion, then check the trimmed mean
    # against what the C++ itself reported rather than trusting the arithmetic.
    burnin_iterations = meta["n_burnin_rounds"] * meta["burnin"]
    total = burnin_iterations + meta["iterations"]
    keep_from = int(round(len(trace) * burnin_iterations / total))
    return trace.iloc[keep_from:].reset_index(drop=True)


def _check_trim(trace: pd.DataFrame, rung_dir: pathlib.Path) -> str | None:
    """Cross-check the burn-in trim against the C++'s own posterior means."""
    mean_var_path = rung_dir / "acol_meanVar.txt"
    if not mean_var_path.exists():
        return "acol_meanVar.txt missing; burn-in trim is unverified."

    reported = pd.read_csv(mean_var_path, sep="\t")
    reported = reported[reported["name"].isin(trace.columns)]
    if reported.empty:
        return "No shared parameters with acol_meanVar.txt; trim is unverified."

    ours = trace[reported["name"]].mean().to_numpy()
    theirs = reported["posterior_mean"].to_numpy()

    # Measured in posterior standard deviations, not as a relative difference: a
    # parameter whose posterior mean sits near zero produces an enormous relative
    # gap from a negligible absolute one, which says nothing about the trim.
    spread = np.maximum(np.sqrt(reported["posterior_variance"].to_numpy()), 1e-12)
    worst = float(np.max(np.abs(ours - theirs) / spread))
    if worst > 0.25:
        return (
            f"Burn-in trim looks wrong: our trimmed mean is up to {worst:.2f} "
            "posterior standard deviations from the C++'s own reported mean. "
            "Coverage below is not trustworthy."
        )
    return None


def _coverage(trace: pd.DataFrame, truth: dict[str, float], prefix: str) -> dict:
    names = [n for n in truth if n.startswith(prefix) and n in trace.columns]
    if not names:
        return {"n": 0}

    values = trace[names].to_numpy()
    tail = (1.0 - CREDIBLE_MASS) / 2.0
    lower, upper = np.quantile(values, [tail, 1.0 - tail], axis=0)
    actual = np.array([truth[n] for n in names])
    mean = values.mean(axis=0)

    inside = (actual >= lower) & (actual <= upper)
    return {
        "n": len(names),
        "coverage": float(inside.mean()),
        "rmse": float(np.sqrt(np.mean((mean - actual) ** 2))),
        "bias": float(np.mean(mean - actual)),
        "correlation": float(np.corrcoef(mean, actual)[0, 1]) if len(names) > 2 else None,
        "moved": float(np.mean([trace[n].nunique() > 1 for n in names])),
    }


def _branch_lengths(trace: pd.DataFrame, truth: dict[str, float], meta: dict) -> dict:
    prefix = "species_branch_lengths_"
    names = [n for n in truth if n.startswith(prefix) and n in trace.columns]
    if not names:
        return {"n": 0}

    values = trace[names].to_numpy()
    mean = values.mean(axis=0)
    actual = np.array([truth[n] for n in names])

    budget = sum(meta["species_bins"])
    sums = values.sum(axis=1)
    return {
        "n": len(names),
        "spearman": float(spearmanr(mean, actual).statistic),
        "mae": float(np.mean(np.abs(mean - actual))),
        "moved": float(np.mean([trace[n].nunique() > 1 for n in names])),
        "budget_expected": int(budget),
        "budget_observed": sorted({int(s) for s in sums}),
        "budget_conserved": bool(np.all(sums == budget)),
    }


def _scalar(trace: pd.DataFrame, truth: dict[str, float], name: str) -> dict:
    if name not in trace.columns or name not in truth:
        return {"present": False}
    values = trace[name].to_numpy()
    tail = (1.0 - CREDIBLE_MASS) / 2.0
    lower, upper = np.quantile(values, [tail, 1.0 - tail])
    return {
        "present": True,
        "truth": truth[name],
        "mean": float(values.mean()),
        "lower": float(lower),
        "upper": float(upper),
        "covered": bool(lower <= truth[name] <= upper),
    }


def _report(summary: dict) -> None:
    click.echo(f"\n=== {summary['rung']} ===")

    for label, key in (("alpha", "alpha"), ("log_nu", "log_nu")):
        block = summary[key]
        if not block["n"]:
            click.echo(f"  {label:<16} absent from trace")
            continue
        click.echo(
            f"  {label:<16} coverage {block['coverage']:6.1%}  "
            f"rmse {block['rmse']:7.4f}  bias {block['bias']:+7.4f}  "
            f"corr {block['correlation']:5.3f}  moved {block['moved']:5.1%}"
            f"   (n={block['n']})"
        )

    branches = summary["branch_lengths"]
    if branches["n"]:
        click.echo(
            f"  {'branch lengths':<16} spearman {branches['spearman']:5.3f}  "
            f"mae {branches['mae']:6.3f}  moved {branches['moved']:5.1%}"
            f"   (n={branches['n']})"
        )
        if not branches["budget_conserved"]:
            click.echo(
                f"    ! branch-length budget not conserved: expected "
                f"{branches['budget_expected']}, saw {branches['budget_observed']}"
            )

    for name in ("species_mean_log_nu", "species_var_log_nu"):
        block = summary[name]
        if not block["present"]:
            continue
        mark = "ok " if block["covered"] else "MISS"
        click.echo(
            f"  {name:<22} {mark} truth {block['truth']:+.4f}  "
            f"mean {block['mean']:+.4f}  "
            f"95% [{block['lower']:+.4f}, {block['upper']:+.4f}]"
        )

    for warning in summary["warnings"]:
        click.echo(f"  ! {warning}")


def _apply_gates(summary: dict, gates_path: pathlib.Path) -> bool:
    """Compare against a reference rung's scores, allowing a margin."""
    reference = json.loads(gates_path.read_text())
    failures = []

    for key in ("alpha", "log_nu"):
        if summary[key]["n"] and reference[key]["n"]:
            if summary[key]["coverage"] < reference[key]["coverage"] - 0.10:
                failures.append(
                    f"{key} coverage {summary[key]['coverage']:.1%} is more than "
                    f"10 points below the reference {reference[key]['coverage']:.1%}"
                )

    here, there = summary["branch_lengths"], reference["branch_lengths"]
    if here["n"] and there["n"] and here["spearman"] < there["spearman"] - 0.15:
        failures.append(
            f"branch-length spearman {here['spearman']:.3f} is well below the "
            f"reference {there['spearman']:.3f}"
        )

    click.echo(f"\n  gates vs {gates_path.name}:")
    for failure in failures:
        click.echo(f"    FAIL {failure}")
    if not failures:
        click.echo("    PASS")
    return not failures


@click.command()
@click.argument("scenario_dir", type=click.Path(exists=True, file_okay=False))
@click.argument("rung")
@click.option("--gates", type=click.Path(exists=True), default=None,
              help="A previous rung's summary JSON to score this one against.")
def main(scenario_dir: str, rung: str, gates: str | None) -> None:
    """Score RUNG's output inside SCENARIO_DIR against the simulated truth."""
    base = pathlib.Path(scenario_dir)
    rung_dir = base / rung
    if not (rung_dir / "acol_trace.txt").exists():
        raise click.ClickException(
            f"No trace at {rung_dir / 'acol_trace.txt'}. Run bash {base}/{rung}.sh first."
        )

    meta = json.loads((base / "meta.json").read_text())
    truth = _read_truth(base / "truth_species.txt")
    trace = _read_trace(rung_dir, meta)

    warnings = []
    if len(trace) < 50:
        warnings.append(
            f"Only {len(trace)} samples survive the burn-in trim. Credible "
            "intervals from this few draws are meaningless; raise --iterations "
            "or check that meta.json matches the run scripts."
        )
    trim_warning = _check_trim(trace, rung_dir)
    if trim_warning:
        warnings.append(trim_warning)

    summary = {
        "rung": rung,
        "scenario": str(base),
        "retained_samples": len(trace),
        "alpha": _coverage(trace, truth, "species_alpha_"),
        "log_nu": _coverage(trace, truth, "species_log_nu_"),
        "branch_lengths": _branch_lengths(trace, truth, meta),
        "species_mean_log_nu": _scalar(trace, truth, "species_mean_log_nu"),
        "species_var_log_nu": _scalar(trace, truth, "species_var_log_nu"),
        "warnings": warnings,
    }

    _report(summary)

    out_path = rung_dir / "validation_summary.json"
    out_path.write_text(json.dumps(summary, indent=2) + "\n")
    click.echo(f"\n  summary written to {out_path}")

    if gates and not _apply_gates(summary, pathlib.Path(gates)):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
