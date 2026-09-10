"""Read the mixing cost of dropping the block update out of the runs `mixing_cost/run.sh` made.

ADR-0005 argued that the leaf layer needs an eight-state block update: at a small
error probability a field cell at one pins both tree fields, and the triple is
metastable under single-variable draws. It also said the failure would present as
slow mixing rather than as a bug. This is the measurement that turns that from an
assumption into a number.

Two things are read, and they answer different questions.

**The joint density trace** says how fast a chain forgets where it was. Every
factor of the ADR-0005 factorisation is a proper density, so the total moves only
with the configuration and the parameters. Its integrated autocorrelation time is
what one independent draw costs.

It is read factor by factor and not only as a total, which is what the columns
are for (`src/field/joint_density.h`). The factors mix at very different speeds:
the two node-state factors carry the phylogenetic parameters, which are the
slowest thing in the model, while the **data** factor is a function of the field
and the data-source parameters alone. So the total is the whole chain's number
and the data factor is the field's own, and the second is the one that answers
the question this measurement asks. Each is reported next to a **drift** -- how
far the second half of the trace sits from the first, in trace standard
deviations -- because an autocorrelation time taken over a chain still climbing
out of its initialisation is a lower bound rather than an answer.

**Both tree field posteriors** say where the chain went. The field's own
posterior cannot: derivation 3 of ADR-0005 puts the field's density at the
product `a~_s * a~_m`, so the two trees can trade against each other at constant
product and leave the field exactly where it was. Two runs that differ along that
ridge and nowhere else are the signature the block update existed to prevent, so
the two densities are reported as a point on it -- a coordinate along the ridge
and one across it.

The data set is the control variable and the seed is not. `mixing_cost/run.sh`
says why.

    uv run python compare_mixing.py <workdir> [--json summary.json]
"""

from __future__ import annotations

import dataclasses
import json
import pathlib
import statistics

import click
import numpy as np

from src.mixing.posteriors import (
    RidgePoint,
    TreeFieldPosterior,
    compare_posteriors,
    read_tree_field_posterior,
    ridge_point,
    ridge_shift,
)
from src.mixing.traces import (
    TraceSummary,
    autocorrelation,
    read_joint_density_factors,
    summarise_trace,
)

REFERENCE, CURRENT = "reference", "current"
BINARIES = (REFERENCE, CURRENT)

# The lags the autocorrelation function is printed at. Far enough apart to show the
# shape of the decay, and short enough that the estimate at the last of them is
# still made from most of the trace.
ACF_LAGS = (1, 2, 5, 10, 25, 50, 100, 200)

# The column that is the whole chain's number, and the one that is the field's own.
# `data` is `log p(L, D | Y)`: a function of the field and the data-source
# parameters, and of nothing the two phylogenies carry.
TOTAL = "joint_density"
FIELD_FACTOR = "data"


def _chains(workdir: pathlib.Path, binary: str) -> list[tuple[str, pathlib.Path]]:
    """Every chain one binary ran, as `(seed, directory)` in seed order."""
    directories = sorted((workdir / binary).glob("seed_*"), key=lambda p: int(p.name.split("_")[1]))
    if not directories:
        raise click.ClickException(f"{workdir / binary} holds no seed_* chain directory.")
    return [(directory.name.split("_")[1], directory) for directory in directories]


def _pooled(posteriors) -> np.ndarray:
    """The replicates' posteriors averaged cell by cell.

    Each replicate counts the same iterations over the same leaf-pair space, so
    the mean over replicates is the posterior a chain of all of them together
    would have reported.
    """
    return np.mean([p.fractions for p in posteriors], axis=0)


def _spread(values: list[float]) -> float:
    """The standard deviation over replicates.

    One replicate has no spread to report, and nan says that in the column rather
    than a zero, which would read as replicates that agreed exactly.
    """
    return statistics.stdev(values) if len(values) > 1 else float("nan")


def _report_traces(chains, meta, report) -> dict[str, dict[str, list[np.ndarray]]]:
    """Every chain's joint density trace: the total per chain, then every factor."""
    prefix, burn_in_rows = meta["prefix"], int(meta["burn_in_rows"])

    click.echo("== the joint density trace, as a total ==")
    click.echo(f"{'binary':<10} {'seed':>5} {'samples':>8} {'mean':>14} {'acf(1)':>8} "
               f"{'tau':>9} {'ESS':>8} {'drift':>7} {'s/chain':>8} {'ESS/s':>8}")
    factors: dict[str, dict[str, list[np.ndarray]]] = {}
    summaries: dict[str, list[TraceSummary]] = {}
    for binary in BINARIES:
        summaries[binary] = []
        factors[binary] = {}
        for seed, directory in chains[binary]:
            columns = read_joint_density_factors(
                directory / f"{prefix}_joint_density.txt", burn_in_rows=burn_in_rows
            )
            for name, values in columns.items():
                factors[binary].setdefault(name, []).append(values)

            summary = summarise_trace(columns[TOTAL])
            summaries[binary].append(summary)
            seconds = float(meta["seconds"][binary][seed])
            click.echo(f"{binary:<10} {seed:>5} {summary.n_samples:>8} {summary.mean:>14.1f} "
                       f"{summary.acf_lag_one:>8.3f} {summary.autocorrelation_time:>9.1f} "
                       f"{summary.effective_sample_size:>8.1f} {summary.drift:>7.2f} "
                       f"{seconds:>8.1f} {summary.effective_sample_size / seconds:>8.3f}")
            report["traces"].setdefault(binary, {})[seed] = {
                "n_samples": summary.n_samples,
                "mean": summary.mean,
                "sd": summary.sd,
                "acf_lag_one": summary.acf_lag_one,
                "autocorrelation_time": summary.autocorrelation_time,
                "effective_sample_size": summary.effective_sample_size,
                "drift": summary.drift,
                "seconds": seconds,
                "effective_samples_per_second": summary.effective_sample_size / seconds,
            }

    _report_factors(factors, report)
    _report_totals(chains, summaries, report)
    return factors


def _report_factors(factors, report) -> None:
    """The same traces, one row per factor of the ADR-0005 factorisation."""
    click.echo("\n== the same trace, factor by factor ==")
    click.echo(f"{'factor':<22} {'binary':<10} {'acf(1)':>8} {'median tau':>11} {'ESS':>8} "
               f"{'drift':>7} {'ratio':>7}")
    report["factors"] = {}
    for name in factors[REFERENCE]:
        medians = {}
        for binary in BINARIES:
            per_chain = [summarise_trace(values) for values in factors[binary][name]]
            medians[binary] = statistics.median([s.autocorrelation_time for s in per_chain])
            report["factors"].setdefault(name, {})[binary] = {
                "acf_lag_one": [s.acf_lag_one for s in per_chain],
                "autocorrelation_time": [s.autocorrelation_time for s in per_chain],
                "effective_sample_size": [s.effective_sample_size for s in per_chain],
                "drift": [s.drift for s in per_chain],
                "median_autocorrelation_time": medians[binary],
            }
            click.echo(
                f"{name if binary == REFERENCE else '':<22} {binary:<10} "
                f"{statistics.mean([s.acf_lag_one for s in per_chain]):>8.3f} "
                f"{medians[binary]:>11.1f} "
                f"{statistics.median([s.effective_sample_size for s in per_chain]):>8.1f} "
                f"{statistics.mean([s.drift for s in per_chain]):>7.2f} "
                f"{medians[binary] / medians[REFERENCE]:>7.2f}"
            )
        report["factors"][name]["ratio_current_against_reference"] = (
            medians[CURRENT] / medians[REFERENCE]
        )

    # The autocorrelation function itself, and not only the time it collapses to.
    # A single tau hides the shape: two chains can carry the same tau with a
    # short sharp decay on one side and a long shallow one on the other.
    for name in (TOTAL, FIELD_FACTOR):
        click.echo(f"\n  the {name} autocorrelation function, averaged over the replicates")
        click.echo("  " + f"{'binary':<10}"
                   + "".join(f"{'lag ' + str(lag):>10}" for lag in ACF_LAGS))
        for binary in BINARIES:
            curves = [autocorrelation(v, max_lag=max(ACF_LAGS)) for v in factors[binary][name]]
            mean_curve = np.mean(curves, axis=0)
            click.echo("  " + f"{binary:<10}"
                       + "".join(f"{mean_curve[lag]:>10.3f}" for lag in ACF_LAGS))
            report.setdefault("acf", {}).setdefault(name, {})[binary] = {
                str(lag): float(mean_curve[lag]) for lag in ACF_LAGS
            }


def _report_totals(chains, summaries, report) -> None:
    """What one independent draw of the whole model costs, in iterations and in seconds."""
    taus = {b: [s.autocorrelation_time for s in summaries[b]] for b in BINARIES}
    ess_per_second = {
        b: [report["traces"][b][seed]["effective_samples_per_second"] for seed, _ in chains[b]]
        for b in BINARIES
    }
    seconds = {
        b: [report["traces"][b][seed]["seconds"] for seed, _ in chains[b]] for b in BINARIES
    }
    click.echo(f"\n  median total tau   reference {statistics.median(taus[REFERENCE]):9.1f}   "
               f"current {statistics.median(taus[CURRENT]):9.1f}   "
               f"ratio {statistics.median(taus[CURRENT]) / statistics.median(taus[REFERENCE]):.2f}")
    click.echo(f"  median ESS/s       reference "
               f"{statistics.median(ess_per_second[REFERENCE]):9.3f}   "
               f"current {statistics.median(ess_per_second[CURRENT]):9.3f}   "
               f"ratio {statistics.median(ess_per_second[CURRENT]) / statistics.median(ess_per_second[REFERENCE]):.2f}")
    click.echo(f"  seconds per chain  reference {statistics.median(seconds[REFERENCE]):9.1f}   "
               f"current {statistics.median(seconds[CURRENT]):9.1f}")
    report["summary"] = {
        "median_autocorrelation_time": {b: statistics.median(taus[b]) for b in BINARIES},
        "median_effective_samples_per_second": {
            b: statistics.median(ess_per_second[b]) for b in BINARIES
        },
    }


def _report_posteriors(chains, meta, trees, report) -> dict[str, dict[str, list]]:
    """Each tree field's posterior density, and the two binaries cell by cell."""
    prefix, n_cells = meta["prefix"], int(meta["n_cells"])

    click.echo("\n== both tree field posteriors ==")
    posteriors: dict[str, dict[str, list]] = {}
    for binary in BINARIES:
        posteriors[binary] = {}
        for tree in trees:
            posteriors[binary][tree] = [
                read_tree_field_posterior(
                    directory / f"{prefix}_{tree}_tree_field_posterior.txt", n_cells
                )
                for _, directory in chains[binary]
            ]

    click.echo(f"{'tree':<12} {'binary':<10} {'density':>9} {'spread':>9}")
    for tree in trees:
        for binary in BINARIES:
            densities = [p.mean for p in posteriors[binary][tree]]
            click.echo(f"{tree:<12} {binary:<10} {statistics.mean(densities):>9.5f} "
                       f"{_spread(densities):>9.5f}")
            report["posteriors"].setdefault(tree, {})[binary] = {
                "densities": densities,
                "mean_density": statistics.mean(densities),
                "spread_over_replicates": _spread(densities),
            }

    click.echo(f"\n{'tree':<12} {'cells':>7} {'mean |diff|':>12} {'max |diff|':>11} "
               f"{'rms':>8} {'corr':>7}")
    for tree in trees:
        # The replicates pooled, so what is compared is each binary's answer and
        # not one chain's noise against another's.
        left = TreeFieldPosterior(_pooled(posteriors[REFERENCE][tree]))
        right = TreeFieldPosterior(_pooled(posteriors[CURRENT][tree]))
        agreement = compare_posteriors(left, right)
        click.echo(f"{tree:<12} {agreement.n_cells:>7} {agreement.mean_abs_diff:>12.5f} "
                   f"{agreement.max_abs_diff:>11.5f} {agreement.rms_diff:>8.5f} "
                   f"{agreement.correlation:>7.4f}")
        report["posteriors"][tree]["agreement"] = dataclasses.asdict(agreement)

    return posteriors


def _report_ridge(chains, posteriors, trees, report) -> None:
    """Where each run sits on the ADR-0005 ridge, and how far the two binaries sit apart on it."""
    click.echo("\n== the ADR-0005 ridge ==")
    click.echo(f"{'binary':<10} {'seed':>5} {'a~_s':>9} {'a~_m':>9} {'product':>10} "
               f"{'along':>9} {'across':>9}")
    points: dict[str, list[RidgePoint]] = {}
    for binary in BINARIES:
        points[binary] = []
        for index, (seed, _) in enumerate(chains[binary]):
            point = ridge_point(posteriors[binary][trees[0]][index],
                                posteriors[binary][trees[1]][index])
            points[binary].append(point)
            click.echo(f"{binary:<10} {seed:>5} {point.species:>9.5f} "
                       f"{point.molecules:>9.5f} {point.product:>10.6f} "
                       f"{point.log_ratio:>9.4f} {point.log_product:>9.4f}")
        report["ridge"][binary] = [
            {"species": p.species, "molecules": p.molecules, "product": p.product,
             "log_ratio": p.log_ratio, "log_product": p.log_product}
            for p in points[binary]
        ]

    # A chain that rides the ridge scatters along it between seeds while the
    # product stays put, so the two spreads are reported apart.
    click.echo(f"\n{'binary':<10} {'spread along':>13} {'spread across':>14}")
    for binary in BINARIES:
        along = _spread([p.log_ratio for p in points[binary]])
        across = _spread([p.log_product for p in points[binary]])
        click.echo(f"{binary:<10} {along:>13.4f} {across:>14.4f}")
        report["ridge"].setdefault("spread", {})[binary] = {"along": along, "across": across}

    centre = {
        b: RidgePoint(
            species=statistics.mean([p.species for p in points[b]]),
            molecules=statistics.mean([p.molecules for p in points[b]]),
        )
        for b in BINARIES
    }
    shift = ridge_shift(centre[CURRENT], centre[REFERENCE])
    click.echo("\n  current against reference, averaged over seeds:")
    click.echo(f"    along the ridge  {shift.along:+.4f}   (the two trees trading)")
    click.echo(f"    across it        {shift.across:+.4f}   (the field's own density moving)")
    report["ridge"]["shift_current_against_reference"] = {
        "along": shift.along, "across": shift.across
    }


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("workdir", type=click.Path(exists=True, file_okay=False, path_type=pathlib.Path))
@click.option("--json", "json_path", type=click.Path(path_type=pathlib.Path),
              help="Write the report as JSON as well as printing it.")
def main(workdir: pathlib.Path, json_path: pathlib.Path | None) -> None:
    meta = json.loads((workdir / "meta.json").read_text())
    trees = list(meta["tree_names"])
    # RidgePoint names its two coordinates `species` and `molecules`, so the order meta.json gives
    # is checked rather than assumed: the ridge would still be a ridge with the two swapped, and the
    # `along` coordinate would silently change sign.
    if trees != ["species", "molecules"]:
        raise click.ClickException(
            f"The ridge is the species tree against the molecule tree, in that order, but "
            f"meta.json names {trees}."
        )

    report: dict = {"meta": meta, "traces": {}, "posteriors": {}, "ridge": {}}
    chains = {binary: _chains(workdir, binary) for binary in BINARIES}

    click.echo(f"The data set is fixed and the seed is not: {meta['n_replicates']} chains per "
               f"binary, {meta['iterations']} iterations each, on one data set simulated by the "
               f"reference binary at seed {meta['simulate_seed']}.")
    click.echo(f"  reference {meta['reference_revision']}   current {meta['current_revision']}")
    click.echo(f"  {meta['n_species_leaves']} x {meta['n_molecules_leaves']} = "
               f"{meta['n_cells']} leaf pairs, omega = {meta['error_probability']}\n")

    _report_traces(chains, meta, report)
    posteriors = _report_posteriors(chains, meta, trees, report)
    _report_ridge(chains, posteriors, trees, report)

    if json_path is not None:
        json_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        click.echo(f"\nwrote {json_path}")


if __name__ == "__main__":
    main()
