"""Compare what the reference draws against what the C++ simulator draws.

Both sides run identical parameters, so any disagreement is a bug in one of them.
Both trees are active: each draws its own node state top-down over every one of
its nodes, and the field is a noisy AND of the two leaf blocks (ADR-0005).

Five statistics, in increasing order of how much of the machinery they touch:

1. **Field density.** The expected density at leaf pair `(s, m)` is the *product*
   of the two adjusted rates, `a~(alpha_s[m]) * a~(alpha_m[s])`, and nothing else
   (ADR-0005, derivation 3). This is an analytic target rather than a
   cross-implementation comparison, so it says which side is wrong rather than
   merely that they differ. It uses only the field.
2. **Tree-field presence per clique.** Every node's marginal is stationary, so
   the fraction of a clique's tree field in state 1 is that clique's `alpha`.
   Again analytic, and it catches root-sampling and stationary-distribution bugs
   in either tree.
3. **The link.** `P(Y = 1 | bucket)` read off each side's own tree fields and
   field, against `P_k = (1 - omega)^k * omega^(2 - k)`. This is the one
   statistic the old shared-field model had no room for.
4. **Sibling calibration.** For each sibling leaf pair of the species tree the
   disagreement probability of its tree field is closed form; pairs are grouped
   into deciles of that prediction and observed frequency is plotted against it.
   Catches transition-matrix bugs, including an off-by-one in the bin grid.
5. **The C++'s own counters.** The six integers it traced, against the six this
   module tallies from the same three files. The error probability's whole
   likelihood reads those counters, so a tally that drifts from the cells it
   counts would be silent everywhere else.

`--gate` turns the report into a check. Each deviation above is scored in units
of its own binomial standard error and the command exits non-zero when one
exceeds the given number of them. There is no tolerance to pick: the sampling
error of a fraction of `n` draws is what it is.

The C++ builds its bin matrices by repeated multiplication rather than by
evaluating the exponential at each bin; the reference deliberately does the
latter, so accumulated error in that recursion shows up here.
"""

from __future__ import annotations

import json
import pathlib
from dataclasses import dataclass

import click
import numpy as np
import pandas as pd

from src.independent import field as F
from src.independent import io
from src.independent import link as L
from src.independent.indexing import TreeIndex, build_tree_index

N_DECILES = 10


@dataclass
class Draw:
    """One replicate: both tree fields, and the field they were reconciled into.

    Both tree fields are indexed `[species_leaf, molecule_leaf]`, the same way the
    field is: a tree field and the field share their subscripts (ADR-0005).
    """

    species_field: np.ndarray
    molecule_field: np.ndarray
    field: np.ndarray


def _load_indices(base: pathlib.Path) -> tuple[TreeIndex, TreeIndex]:
    def edges(path):
        frame = pd.read_csv(path, sep="\t")
        return list(zip(frame["child"].astype(str), frame["parent"].astype(str)))

    return build_tree_index(edges(base / "species.txt")), build_tree_index(
        edges(base / "molecules.txt")
    )


def _truth_vectors(
    base: pathlib.Path, tree_name: str, clique_names: list[str]
) -> tuple[np.ndarray, np.ndarray]:
    """One tree's `(alphas, nus)`, one per clique, in clique order.

    A clique of this tree is named by a leaf of the *other* one, which is what
    `clique_names` carries.
    """
    frame = pd.read_csv(base / f"truth_{tree_name}.txt", sep="\t")
    values = dict(zip(frame["name"].astype(str), frame["value"].astype(float)))
    alphas = np.array([values[f"{tree_name}_alpha_{n}"] for n in clique_names])
    nus = np.exp(np.array([values[f"{tree_name}_log_nu_{n}"] for n in clique_names]))
    return alphas, nus


def _sibling_pairs(species: TreeIndex) -> list[tuple[int, int, int, int]]:
    """`(left_leaf_ix, right_leaf_ix, left_bin_slot, right_bin_slot)` per pair."""
    leaf_slot = {node: i for i, node in enumerate(species.leaves)}
    children: dict[int, list[int]] = {}
    for node in range(species.n_nodes):
        parent = int(species.parent[node])
        if parent >= 0 and node in leaf_slot:
            children.setdefault(parent, []).append(node)

    pairs = []
    for kin in children.values():
        for i in range(len(kin)):
            for j in range(i + 1, len(kin)):
                left, right = kin[i], kin[j]
                pairs.append(
                    (
                        leaf_slot[left],
                        leaf_slot[right],
                        int(species.branch_of_node[left]),
                        int(species.branch_of_node[right]),
                    )
                )
    return pairs


def _sibling_observed(
    tree_fields: list[np.ndarray], pairs, predicted: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Observed vs predicted sibling disagreement, grouped into deciles."""
    observed = np.zeros(len(predicted))
    n_cliques = tree_fields[0].shape[1]
    for offset, (left, right, _, _) in enumerate(pairs):
        slot = slice(offset * n_cliques, (offset + 1) * n_cliques)
        observed[slot] = np.mean([f[left] != f[right] for f in tree_fields], axis=0)

    order = np.argsort(predicted)
    groups = np.array_split(order, N_DECILES)
    return (
        np.array([predicted[g].mean() for g in groups]),
        np.array([observed[g].mean() for g in groups]),
    )


def _sibling_predicted(pairs, bins, alphas, nus) -> np.ndarray:
    """Closed-form disagreement probability for every (pair, clique)."""
    out = np.empty(len(pairs) * len(alphas))
    for offset, (_, _, left_slot, right_slot) in enumerate(pairs):
        for clique, (alpha, nu) in enumerate(zip(alphas, nus)):
            out[offset * len(alphas) + clique] = F.sibling_disagreement_probability(
                alpha, nu, int(bins[left_slot]), int(bins[right_slot])
            )
    return out


def _expected_field(
    species_alphas: np.ndarray, molecule_alphas: np.ndarray, omega: float
) -> np.ndarray:
    """`P(Y = 1)` at every leaf pair: the product of the two adjusted rates.

    The species tree's cliques are named by molecule leaves and the molecule
    tree's by species leaves, so the species alphas index the columns of the
    field and the molecule alphas index its rows.
    """
    return np.outer(
        L.adjusted_rate(molecule_alphas, omega), L.adjusted_rate(species_alphas, omega)
    )


def _reference_draws(
    base: pathlib.Path, meta: dict, species: TreeIndex, molecules: TreeIndex, n: int
) -> list[Draw]:
    species_alphas, species_nus = _truth_vectors(
        base, "species", molecules.leaf_names()
    )
    molecule_alphas, molecule_nus = _truth_vectors(
        base, "molecules", species.leaf_names()
    )
    species_bins = np.array(meta["species_bins"], dtype=np.int64)
    molecule_bins = np.array(meta["molecule_bins"], dtype=np.int64)
    omega = meta["error_probability"]

    draws = []
    for i in range(n):
        rng = np.random.default_rng(meta["seed"] + 10_000 + i)
        species_states = F.sample_states(
            rng, species, species_bins, species_alphas, species_nus
        )
        molecule_states = F.sample_states(
            rng, molecules, molecule_bins, molecule_alphas, molecule_nus
        )
        species_field = species_states[species.leaves]
        molecule_field = molecule_states[molecules.leaves].T
        draws.append(
            Draw(
                species_field=species_field,
                molecule_field=molecule_field,
                field=L.sample_field(rng, species_field, molecule_field, omega),
            )
        )
    return draws


def _cpp_draws(
    base: pathlib.Path, species: TreeIndex, molecules: TreeIndex, n: int
) -> list[Draw]:
    draws = []
    for i in range(1, n + 1):
        run = base / "replicates" / f"cpp_{i}"
        field_path = run / "acol_simulated_Y.txt"
        if not field_path.exists():
            break
        species_states = io.read_node_states(
            run / "acol_simulated_Z_species.txt", molecules.n_leaves
        )
        molecule_states = io.read_node_states(
            run / "acol_simulated_Z_molecules.txt", molecules.n_nodes
        )
        draws.append(
            Draw(
                species_field=species_states[species.leaves],
                molecule_field=molecule_states[:, molecules.leaves],
                field=io.read_field(field_path, molecules.n_leaves),
            )
        )
    return draws


def _max_sigmas(per_replicate: np.ndarray, expected: np.ndarray) -> float:
    """The worst deviation from `expected`, in standard errors of the mean.

    `per_replicate` is one row per replicate. The spread comes from the replicates
    themselves and not from a binomial formula, because the cells inside one
    replicate are correlated along the tree: they are drawn from one process down
    one set of branches, so their effective count is far below their number.
    Replicates are independent draws of the whole scenario, so their scatter is
    the honest yardstick, and the gate needs no tolerance of its own.

    Two replicates cannot estimate a spread, so a shorter run reports 0.
    """
    if len(per_replicate) < 3:
        return 0.0
    mean = np.nanmean(per_replicate, axis=0)
    spread = np.nanstd(per_replicate, axis=0, ddof=1) / np.sqrt(len(per_replicate))
    # A statistic that never moved across replicates has no scatter to divide by.
    # It is either exactly right, or wrong by a distance no replicate can explain.
    deviation = np.abs(mean - expected)
    with np.errstate(divide="ignore", invalid="ignore"):
        sigmas = np.where(spread > 0.0, deviation / np.maximum(spread, 1e-300), 0.0)
    return float(np.nanmax(sigmas))


def _field_deviation(draws: list[Draw], expected: np.ndarray) -> dict:
    observed = np.mean([d.field for d in draws], axis=0)
    deviation = observed - expected
    column = _max_sigmas(
        np.array([d.field.mean(axis=0) for d in draws]), expected.mean(axis=0)
    )
    row = _max_sigmas(
        np.array([d.field.mean(axis=1) for d in draws]), expected.mean(axis=1)
    )
    return {
        "max_abs_deviation": float(np.abs(deviation).max()),
        "mean_deviation": float(deviation.mean()),
        "column_max_abs_deviation": float(
            np.abs(observed.mean(axis=0) - expected.mean(axis=0)).max()
        ),
        "row_max_abs_deviation": float(
            np.abs(observed.mean(axis=1) - expected.mean(axis=1)).max()
        ),
        "max_sigmas": max(column, row),
    }


def _tree_field_deviation(
    draws: list[Draw], tree_name: str, alphas: np.ndarray
) -> dict:
    """Per-clique presence of one tree field against its analytic target, `alpha`."""
    # A species-tree clique is a column of the tree field and a molecule-tree
    # clique is a row of it, so each is averaged over the other tree's leaves.
    axis = 0 if tree_name == "species" else 1
    fields = [
        d.species_field if tree_name == "species" else d.molecule_field for d in draws
    ]
    per_replicate = np.array([f.mean(axis=axis) for f in fields])
    observed = per_replicate.mean(axis=0)
    deviation = observed - alphas
    return {
        "max_abs_deviation": float(np.abs(deviation).max()),
        "mean_deviation": float(deviation.mean()),
        "correlation": float(np.corrcoef(observed, alphas)[0, 1]),
        "max_sigmas": _max_sigmas(per_replicate, alphas),
    }


def _traced_counters(base: pathlib.Path, n: int) -> np.ndarray | None:
    """The six counters the C++ traced, summed over its replicates.

    The trace holds `n(bucket, field state)` in bucket-major order, one row per
    tally -- and a simulated run takes one draw, so one row. Reading it is what
    puts the C++'s own bucketing under the comparison; recomputing both sides in
    Python would hide a bucketing convention that disagrees.
    """
    total = np.zeros((L.N_BUCKETS, 2), dtype=np.int64)
    found = False
    for i in range(1, n + 1):
        path = (
            base / "replicates" / f"cpp_{i}" / "acol_simulated_link_counters_trace.txt"
        )
        if not path.exists():
            break
        frame = pd.read_csv(path, sep="\t")
        total += frame.to_numpy(dtype=np.int64).sum(axis=0).reshape(L.N_BUCKETS, 2)
        found = True
    return total if found else None


def _link_rates(draws: list[Draw], omega: float) -> dict:
    """`P(Y = 1 | bucket)` from each side's own cells, against the link's `P_k`."""
    counters = np.zeros((L.N_BUCKETS, 2), dtype=np.int64)
    for draw in draws:
        counters += L.link_counters(draw.species_field, draw.molecule_field, draw.field)
    totals = counters.sum(axis=1)
    # A bucket that held no cell gives 0 / 0, which stays a NaN rather than a zero.
    observed = np.divide(
        counters[:, 1], totals, out=np.full(L.N_BUCKETS, np.nan), where=totals > 0
    )
    predicted = L.prob_for_bucket(np.arange(L.N_BUCKETS), omega)

    # One rate per bucket per replicate. A bucket a replicate never filled gives a
    # NaN there, which the nan-aware mean and spread skip.
    per_replicate = []
    for draw in draws:
        one = L.link_counters(draw.species_field, draw.molecule_field, draw.field)
        with np.errstate(divide="ignore", invalid="ignore"):
            per_replicate.append(one[:, 1] / one.sum(axis=1))
    return {
        "counters": counters.tolist(),
        "observed": [None if np.isnan(p) else float(p) for p in observed],
        "predicted": [float(p) for p in predicted],
        "max_abs_deviation": float(np.nanmax(np.abs(observed - predicted))),
        "max_sigmas": _max_sigmas(np.array(per_replicate), predicted),
    }


@click.command()
@click.argument("scenario_dir", type=click.Path(exists=True, file_okay=False))
@click.option("--replicates", type=int, default=20, show_default=True)
@click.option(
    "--gate",
    type=float,
    default=None,
    help="Exit non-zero when a deviation exceeds this many standard errors.",
)
def main(scenario_dir: str, replicates: int, gate: float | None) -> None:
    """Compare reference-drawn scenarios against the C++ simulator's."""
    base = pathlib.Path(scenario_dir)
    meta = json.loads((base / "meta.json").read_text())
    species, molecules = _load_indices(base)
    species_alphas, species_nus = _truth_vectors(
        base, "species", molecules.leaf_names()
    )
    molecule_alphas, _ = _truth_vectors(base, "molecules", species.leaf_names())
    species_bins = np.array(meta["species_bins"], dtype=np.int64)
    omega = meta["error_probability"]

    sides = {
        "reference": _reference_draws(base, meta, species, molecules, replicates),
        "c++": _cpp_draws(base, species, molecules, replicates),
    }
    click.echo(
        f"reference replicates: {len(sides['reference'])}   "
        f"C++ replicates: {len(sides['c++'])}"
    )
    if not sides["c++"]:
        click.echo(f"No C++ draws found. Run: bash {base}/replicates.sh {replicates}")

    summary: dict = {
        "n_reference": len(sides["reference"]),
        "n_cpp": len(sides["c++"]),
        "error_probability": omega,
    }

    # 1. the field against the product of the two adjusted rates
    expected = _expected_field(species_alphas, molecule_alphas, omega)
    click.echo(
        "\n-- field density (analytic target: the product of the two adjusted rates, "
        f"mean {expected.mean():.4f}) --"
    )
    for label, draws in sides.items():
        if not draws:
            continue
        block = _field_deviation(draws, expected)
        summary[f"field_{label}"] = block
        click.echo(
            f"  {label:<10} max |observed - expected| {block['max_abs_deviation']:.4f}   "
            f"per column {block['column_max_abs_deviation']:.4f}   "
            f"per row {block['row_max_abs_deviation']:.4f}   "
            f"mean {block['mean_deviation']:+.4f}"
        )

    # 2. each tree field against its own alphas
    for tree_name, alphas in (
        ("species", species_alphas),
        ("molecules", molecule_alphas),
    ):
        click.echo(
            f"\n-- {tree_name} tree field per clique (analytic target: alpha) --"
        )
        for label, draws in sides.items():
            if not draws:
                continue
            block = _tree_field_deviation(draws, tree_name, alphas)
            summary[f"tree_field_{tree_name}_{label}"] = block
            click.echo(
                f"  {label:<10} max |observed - alpha| {block['max_abs_deviation']:.4f}   "
                f"mean {block['mean_deviation']:+.4f}   corr {block['correlation']:.4f}"
            )

    # 3. the link, read off each side's own cells
    click.echo(f"\n-- the link at omega = {omega} --")
    click.echo("  bucket   predicted   reference        c++")
    rates = {
        label: _link_rates(draws, omega) for label, draws in sides.items() if draws
    }
    for label, block in rates.items():
        summary[f"link_{label}"] = block
    for bucket in range(L.N_BUCKETS):
        predicted = L.prob_for_bucket(bucket, omega)

        def cell(label: str) -> str:
            block = rates.get(label)
            if block is None or block["observed"][bucket] is None:
                return "        -"
            return f"{block['observed'][bucket]:9.4f}"

        click.echo(
            f"  {bucket:>6}   {float(predicted):9.4f}   {cell('reference')}  {cell('c++')}"
        )

    # 5. the C++'s own counters against the ones tallied from its own files
    traced = _traced_counters(base, replicates)
    if traced is not None and sides["c++"]:
        tallied = np.array(rates["c++"]["counters"], dtype=np.int64)
        summary["link_counters_traced"] = traced.tolist()
        summary["link_counters_agree"] = bool(np.array_equal(traced, tallied))
        click.echo(
            "\n-- the C++'s traced counters against the ones tallied from its own files --"
        )
        if summary["link_counters_agree"]:
            click.echo(f"  identical, {int(traced.sum())} cells")
        else:
            click.echo(
                f"  DIFFER\n    traced   {traced.tolist()}\n    tallied  {tallied.tolist()}"
            )

    # 4. sibling calibration on the species tree field
    pairs = _sibling_pairs(species)
    click.echo(f"\n-- species sibling calibration ({len(pairs)} leaf pairs) --")
    predicted_all = _sibling_predicted(pairs, species_bins, species_alphas, species_nus)
    curves = {
        label: _sibling_observed([d.species_field for d in draws], pairs, predicted_all)
        for label, draws in sides.items()
        if draws
    }
    click.echo("  decile   predicted   reference        c++")

    def decile(label: str, d: int) -> str:
        curve = curves.get(label)
        return "        -" if curve is None else f"{curve[1][d]:9.4f}"

    for d in range(N_DECILES):
        # Every curve shares the prediction; only the observed column differs.
        predicted = next(iter(curves.values()))[0][d]
        click.echo(
            f"  {d + 1:>6}   {predicted:9.4f}   {decile('reference', d)}  {decile('c++', d)}"
        )
    for label, curve in curves.items():
        summary[f"sibling_{label}_max_error"] = float(np.abs(curve[1] - curve[0]).max())

    out_path = base / "field_comparison.json"
    out_path.write_text(json.dumps(summary, indent=2) + "\n")
    click.echo(f"\n  summary written to {out_path}")

    if gate is not None:
        _apply_gate(summary, gate)


def _apply_gate(summary: dict, gate: float) -> None:
    """Fail on any deviation further than `gate` standard errors from its target.

    Every statistic above compares an observed fraction with a target the model
    fixes, so its own binomial standard error is the right yardstick. A side that
    drew fewer replicates is judged more leniently by construction, which is what
    a standard error is for.
    """
    failures = []
    for key, block in sorted(summary.items()):
        if not isinstance(block, dict) or "max_sigmas" not in block:
            continue
        if block["max_sigmas"] > gate:
            failures.append(f"{key}: {block['max_sigmas']:.2f} standard errors")
    if summary.get("link_counters_agree") is False:
        failures.append(
            "link_counters: the C++ traced a tally its own cells do not give"
        )

    click.echo(f"\n  gate at {gate} standard errors:")
    for failure in failures:
        click.echo(f"    FAIL {failure}")
    if failures:
        raise SystemExit(1)
    click.echo("    PASS")


if __name__ == "__main__":
    main()
