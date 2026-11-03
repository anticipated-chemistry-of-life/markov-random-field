"""Score one rung's inference output against the independently simulated truth.

Both trees are scored. Each draws its own node state and each is inferred, so a
validator that read one of them would be blind to half the model (ADR-0005).

Reports, and does not by default gate on:

- **coverage** for `alpha` and `log_nu` — the fraction of a tree's cliques whose
  95% credible interval contains the truth. Coverage catches both bias and
  over-confident posteriors; RMSE alone catches only the first.
- **Spearman and mean absolute error** for branch lengths, which are discrete and
  budget-constrained, so a credible interval is a strange object for them.
- **movement**, because a chain that never left its initialisation can otherwise
  masquerade as recovery when the branch lengths start flat.

On the LOTUS+MS scenarios it also scores the mass-spec half of the model:

- **field**: AUC and a reliability table of P(Y=1) against the true Y, for species
  with runs and species without. Reported, no threshold.
- **assignment**: accuracy of the assignment posterior mode, the mean posterior on the
  true assignment, and the mean unknown posterior for unknown features and for
  features of a detected molecule. The posterior's format is checked first.
- **filter probabilities**: 90% credible-interval coverage, stratified by the number
  of runs the molecule is present in.
- **kappa**: coverage.

A coverage passes inside the 95% binomial band around 0.9 (`recovery.coverage_band`).

Thresholds are deliberately absent: rung 1 pins the field *and* the internal
states, so whatever it scores is the empirical ceiling, and the gates for the
looser rungs should be set relative to that rather than guessed in advance. Pass
`--gates <rung1-summary.json>` once you have it.
"""

from __future__ import annotations

import json
import pathlib
from dataclasses import asdict

import click
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from src.independent import mass_spec, recovery

CREDIBLE_MASS = 0.95


TREES = ("species", "molecules")

# Each tree's own bins, under the name meta.json gives them.
BINS_KEY = {"species": "species_bins", "molecules": "molecule_bins"}


def _read_truth(base: pathlib.Path) -> dict[str, float]:
    """Every tree's truth, from the one file per tree the scenario writes."""
    truth: dict[str, float] = {}
    for tree in TREES:
        path = base / f"truth_{tree}.txt"
        if not path.exists():
            continue
        frame = pd.read_csv(path, sep="\t")
        truth.update(zip(frame["name"].astype(str), frame["value"].astype(float)))
    if not truth:
        raise click.ClickException(f"No truth_<tree>.txt in {base}.")
    return truth


def _read_trace(rung_dir: pathlib.Path, meta: dict) -> pd.DataFrame:
    """The main trace plus each tree's own, joined column-wise.

    Branch lengths are updated by the tree rather than by a stattools updater, so
    they are written to `acol_<tree>_trace.txt` and appear nowhere in
    `acol_trace.txt`. Every file shares a thinning factor and row count.
    """
    frames = [pd.read_csv(rung_dir / "acol_trace.txt", sep="\t")]

    for tree in TREES:
        tree_trace = rung_dir / f"acol_{tree}_trace.txt"
        if tree_trace.exists():
            frames.append(pd.read_csv(tree_trace, sep="\t"))

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
        "correlation": float(np.corrcoef(mean, actual)[0, 1])
        if len(names) > 2
        else None,
        "moved": float(np.mean([trace[n].nunique() > 1 for n in names])),
    }


def _branch_lengths(
    trace: pd.DataFrame, truth: dict[str, float], meta: dict, tree: str
) -> dict:
    prefix = f"{tree}_branch_lengths_"
    names = [n for n in truth if n.startswith(prefix) and n in trace.columns]
    if not names:
        return {"n": 0}

    values = trace[names].to_numpy()
    mean = values.mean(axis=0)
    actual = np.array([truth[n] for n in names])

    budget = sum(meta[BINS_KEY[tree]])
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


def _assignment_scores(
    base: pathlib.Path, rung_dir: pathlib.Path
) -> recovery.AssignmentScore:
    """Check the assignment posterior against the data file, then score it."""
    path = rung_dir / "acol_assignment_posterior.txt"
    if not path.exists():
        raise click.ClickException(
            f"No assignment posterior at {path}. The rung must pass "
            "--write_assignment_posterior."
        )
    names = {column: str for column in (*recovery.FEATURE_KEY, "molecule")}
    # `fraction` stays unparsed here: a bad cell is reported by feature, not by pandas.
    posterior = pd.read_csv(path, sep="\t", dtype=names)
    data = pd.read_csv(
        base / mass_spec.DATA_FILE, sep="\t", usecols=list(names), dtype=names
    )
    truth = pd.read_csv(base / mass_spec.ASSIGNMENT_FILE, sep="\t", dtype=names)
    try:
        recovery.check_assignment_format(data, posterior)
        return recovery.assignment_scores(truth, posterior)
    except recovery.AssignmentFormatError as error:
        raise click.ClickException(f"{path}: {error}") from error


def _mass_spec_scores(
    base: pathlib.Path, rung_dir: pathlib.Path, meta: dict, trace: pd.DataFrame
) -> dict:
    """Field, filter-probability and kappa recovery on a LOTUS+MS scenario."""
    truth_y = pd.read_csv(base / "simulated_Y.txt", sep="\t")
    posterior_path = rung_dir / "acol_Y_posterior.tsv"
    if not posterior_path.exists():
        raise click.ClickException(f"No Y posterior at {posterior_path}.")
    posterior = pd.read_csv(posterior_path, sep="\t")

    data = pd.read_csv(base / mass_spec.DATA_FILE, sep="\t", usecols=["species"])
    species_with_runs = set(data["species"].astype(str))
    runs_per_species = meta["mass_spec"]["n_runs"] // len(species_with_runs)

    parameters = pd.read_csv(base / "simulated_parameters.txt", sep="\t")
    parameters = dict(zip(parameters["name"].astype(str), parameters["value"]))
    prefix = f"{mass_spec.FILTER_PARAMETER}_"
    filters = [
        n for n in parameters if n.startswith(prefix) and n in trace.columns
    ]
    molecules = [n[len(prefix) :] for n in filters]

    field = recovery.field_scores(truth_y, posterior, species_with_runs)
    filter_rows = recovery.coverage_by_runs(
        trace[filters].to_numpy(),
        np.array([float(parameters[n]) for n in filters]),
        recovery.runs_present(truth_y, species_with_runs, runs_per_species, molecules),
    )

    kappa_name = mass_spec.CONTAMINATION_PARAMETER
    kappa: dict = {"present": kappa_name in trace.columns}
    if kappa["present"]:
        draws = trace[[kappa_name]].to_numpy()
        truth = float(parameters[kappa_name])
        tail = (1.0 - recovery.CREDIBLE_MASS) / 2.0
        lower, upper = np.quantile(draws[:, 0], [tail, 1.0 - tail])
        covered = bool(recovery.is_covered(draws, np.array([truth]))[0])
        # One interval is one Bernoulli trial, and every band at n=1 is [0, 1]: a
        # miss would pass. So kappa passes when its interval holds the truth.
        kappa.update(
            truth=truth,
            mean=float(draws.mean()),
            lower=float(lower),
            upper=float(upper),
            covered=covered,
            passes=covered,
        )

    assignment = _assignment_scores(base, rung_dir)

    return {
        "field": {key: asdict(score) for key, score in field.items()},
        "assignment": asdict(assignment),
        "filter_probabilities": {
            "n": len(filters),
            "strata": [asdict(row) for row in filter_rows],
        },
        "kappa": kappa,
    }


def _report_mass_spec(scores: dict) -> None:
    click.echo("\n  -- mass spec --")
    for key, label in (("with_runs", "with runs"), ("without_runs", "without runs")):
        block = scores["field"][key]
        click.echo(
            f"  field, species {label:<12} AUC {block['auc']:5.3f}   (n={block['n']} cells)"
        )
        click.echo(f"    {'P(Y=1) bin':<12} {'n':>7} {'mean pred':>10} {'observed':>9}")
        for row in block["reliability"]:
            click.echo(
                f"    [{row['low']:.1f}, {row['high']:.1f}) {row['n']:>7} "
                f"{row['predicted']:>10.3f} {row['observed']:>9.3f}"
            )

    assignment = scores["assignment"]
    click.echo(f"  assignment posterior (n={assignment['n']} features):")
    click.echo(f"    accuracy of the posterior mode        {assignment['accuracy']:6.1%}")
    click.echo(f"    mean posterior on the true assignment {assignment['mean_on_true_assignment']:6.3f}")
    click.echo(
        f"    mean unknown posterior, unknown features   {assignment['mean_unknown_of_unknown_features']:6.3f}"
        f"   (n={assignment['n_unknown_features']})"
    )
    click.echo(
        f"    mean unknown posterior, molecule features  {assignment['mean_unknown_of_molecule_features']:6.3f}"
        f"   (n={assignment['n_molecule_features']})"
    )

    click.echo("  filter probabilities, 90% interval coverage by runs present:")
    for row in scores["filter_probabilities"]["strata"]:
        low, high = recovery.coverage_band(row["n"])
        click.echo(
            f"    runs {row['label']:<5} {'PASS' if row['passes'] else 'FAIL'}  "
            f"coverage {row['coverage']:6.1%}  band [{low:.1%}, {high:.1%}]  "
            f"(n={row['n']})"
        )

    kappa = scores["kappa"]
    if kappa["present"]:
        click.echo(
            f"  kappa  {'PASS' if kappa['passes'] else 'FAIL'}  truth {kappa['truth']:.4f}  "
            f"mean {kappa['mean']:.4f}  90% [{kappa['lower']:.4f}, {kappa['upper']:.4f}]"
        )


def _report(summary: dict) -> None:
    click.echo(f"\n=== {summary['rung']} ===")

    for tree in TREES:
        click.echo(f"\n  -- {tree} --")
        for label in ("alpha", "log_nu"):
            block = summary[tree][label]
            if not block["n"]:
                click.echo(f"  {label:<16} absent from trace")
                continue
            click.echo(
                f"  {label:<16} coverage {block['coverage']:6.1%}  "
                f"rmse {block['rmse']:7.4f}  bias {block['bias']:+7.4f}  "
                f"corr {block['correlation']:5.3f}  moved {block['moved']:5.1%}"
                f"   (n={block['n']})"
            )

        branches = summary[tree]["branch_lengths"]
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

        for scalar in ("mean_log_nu", "var_log_nu"):
            block = summary[tree][scalar]
            if not block["present"]:
                continue
            mark = "ok " if block["covered"] else "MISS"
            click.echo(
                f"  {tree + '_' + scalar:<22} {mark} truth {block['truth']:+.4f}  "
                f"mean {block['mean']:+.4f}  "
                f"95% [{block['lower']:+.4f}, {block['upper']:+.4f}]"
            )

    if "mass_spec" in summary:
        _report_mass_spec(summary["mass_spec"])

    for warning in summary["warnings"]:
        click.echo(f"  ! {warning}")


def _apply_gates(summary: dict, gates_path: pathlib.Path) -> bool:
    """Compare against a reference rung's scores, allowing a margin."""
    reference = json.loads(gates_path.read_text())
    failures = []

    for tree in TREES:
        for key in ("alpha", "log_nu"):
            here, there = summary[tree][key], reference[tree][key]
            if here["n"] and there["n"] and here["coverage"] < there["coverage"] - 0.10:
                failures.append(
                    f"{tree} {key} coverage {here['coverage']:.1%} is more than "
                    f"10 points below the reference {there['coverage']:.1%}"
                )

        here, there = summary[tree]["branch_lengths"], reference[tree]["branch_lengths"]
        if here["n"] and there["n"] and here["spearman"] < there["spearman"] - 0.15:
            failures.append(
                f"{tree} branch-length spearman {here['spearman']:.3f} is well below "
                f"the reference {there['spearman']:.3f}"
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
@click.option(
    "--gates",
    type=click.Path(exists=True),
    default=None,
    help="A previous rung's summary JSON to score this one against.",
)
def main(scenario_dir: str, rung: str, gates: str | None) -> None:
    """Score RUNG's output inside SCENARIO_DIR against the simulated truth."""
    base = pathlib.Path(scenario_dir)
    rung_dir = base / rung
    if not (rung_dir / "acol_trace.txt").exists():
        raise click.ClickException(
            f"No trace at {rung_dir / 'acol_trace.txt'}. Run bash {base}/{rung}.sh first."
        )

    meta = json.loads((base / "meta.json").read_text())
    truth = _read_truth(base)
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
        "warnings": warnings,
    }
    for tree in TREES:
        summary[tree] = {
            "alpha": _coverage(trace, truth, f"{tree}_alpha_"),
            "log_nu": _coverage(trace, truth, f"{tree}_log_nu_"),
            "branch_lengths": _branch_lengths(trace, truth, meta, tree),
            "mean_log_nu": _scalar(trace, truth, f"{tree}_mean_log_nu"),
            "var_log_nu": _scalar(trace, truth, f"{tree}_var_log_nu"),
        }

    if "mass_spec" in meta:
        summary["mass_spec"] = _mass_spec_scores(base, rung_dir, meta, trace)

    _report(summary)

    out_path = rung_dir / "validation_summary.json"
    out_path.write_text(json.dumps(summary, indent=2) + "\n")
    click.echo(f"\n  summary written to {out_path}")

    if gates and not _apply_gates(summary, pathlib.Path(gates)):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
