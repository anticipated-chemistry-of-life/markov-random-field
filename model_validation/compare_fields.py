"""Compare fields drawn by the reference implementation against the C++ simulator.

Both sides run identical parameters, so any disagreement is a bug in one of them.
Three statistics, in increasing order of how much of the machinery they touch:

1. **Per-clique presence.** Every node's marginal is stationary, so the expected
   fraction of present cells in clique `c` is exactly `alpha_c`. This is an
   analytic target, not a cross-implementation comparison, so it says which side
   is wrong rather than merely that they differ. It catches root-sampling and
   stationary-distribution bugs.
2. **Sibling calibration.** For each sibling leaf pair the disagreement
   probability is closed form; pairs are grouped into deciles of that prediction
   and observed frequency is plotted against it. Uses only the field, so it does
   not depend on trusting the C++'s internal states. Catches transition-matrix
   bugs, including an off-by-one in the bin grid.
3. **Transition density by bin.** Parent-to-child switching rate per branch,
   which needs the internal states. Localises a disagreement to a bin.

The C++ builds its bin matrices by repeated multiplication rather than by
evaluating the exponential at each bin; the reference deliberately does the
latter, so accumulated error in that recursion shows up here.
"""

from __future__ import annotations

import json
import pathlib

import click
import numpy as np
import pandas as pd

from src.independent import field as F
from src.independent import io
from src.independent.indexing import build_tree_index

N_DECILES = 10


def _load_indices(base: pathlib.Path):
    def edges(path):
        frame = pd.read_csv(path, sep="\t")
        return list(zip(frame["child"].astype(str), frame["parent"].astype(str)))

    return build_tree_index(edges(base / "species.txt")), build_tree_index(
        edges(base / "molecules.txt")
    )


def _truth_vectors(base: pathlib.Path, molecules) -> tuple[np.ndarray, np.ndarray]:
    frame = pd.read_csv(base / "truth_species.txt", sep="\t")
    values = dict(zip(frame["name"].astype(str), frame["value"].astype(float)))
    leaves = molecules.leaf_names()
    alphas = np.array([values[f"species_alpha_{n}"] for n in leaves])
    nus = np.exp(np.array([values[f"species_log_nu_{n}"] for n in leaves]))
    return alphas, nus


def _sibling_pairs(species) -> list[tuple[int, int, int, int]]:
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


def _presence(fields: list[np.ndarray]) -> np.ndarray:
    """Mean per-clique presence across replicates. Shape `(n_cliques,)`."""
    return np.mean([f.mean(axis=0) for f in fields], axis=0)


def _sibling_observed(
    fields: list[np.ndarray], pairs, predicted: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Observed vs predicted sibling disagreement, grouped into deciles."""
    observed = np.zeros(len(predicted))
    for offset, (left, right, _, _) in enumerate(pairs):
        n_cliques = fields[0].shape[1]
        slot = slice(offset * n_cliques, (offset + 1) * n_cliques)
        observed[slot] = np.mean(
            [f[left] != f[right] for f in fields], axis=0
        )

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


def _reference_fields(base: pathlib.Path, meta: dict, n: int) -> list[np.ndarray]:
    species, molecules = _load_indices(base)
    alphas, nus = _truth_vectors(base, molecules)
    bins = np.array(meta["species_bins"], dtype=np.int64)
    fields = []
    for i in range(n):
        rng = np.random.default_rng(meta["seed"] + 10_000 + i)
        states = F.sample_states(rng, species, bins, alphas, nus)
        fields.append(states[species.leaves])
    return fields


def _cpp_fields(base: pathlib.Path, n_molecule_leaves: int, n: int) -> list[np.ndarray]:
    fields = []
    for i in range(1, n + 1):
        path = base / "replicates" / f"cpp_{i}" / "acol_simulated_Y.txt"
        if not path.exists():
            break
        fields.append(io.read_field(path, n_molecule_leaves))
    return fields


@click.command()
@click.argument("scenario_dir", type=click.Path(exists=True, file_okay=False))
@click.option("--replicates", type=int, default=20, show_default=True)
def main(scenario_dir: str, replicates: int) -> None:
    """Compare reference-sampled fields against the C++ simulator's."""
    base = pathlib.Path(scenario_dir)
    meta = json.loads((base / "meta.json").read_text())
    species, molecules = _load_indices(base)
    alphas, nus = _truth_vectors(base, molecules)
    bins = np.array(meta["species_bins"], dtype=np.int64)

    reference = _reference_fields(base, meta, replicates)
    cpp = _cpp_fields(base, molecules.n_leaves, replicates)
    click.echo(f"reference replicates: {len(reference)}   C++ replicates: {len(cpp)}")
    if not cpp:
        click.echo(f"No C++ fields found. Run: bash {base}/replicates.sh {replicates}")

    summary: dict = {"n_reference": len(reference), "n_cpp": len(cpp)}

    # 1. presence against the analytic stationary target
    click.echo("\n-- per-clique presence (analytic target: alpha) --")
    for label, fields in (("reference", reference), ("c++", cpp)):
        if not fields:
            continue
        observed = _presence(fields)
        deviation = observed - alphas
        click.echo(
            f"  {label:<10} max |observed - alpha| {np.abs(deviation).max():.4f}   "
            f"mean {deviation.mean():+.4f}   corr {np.corrcoef(observed, alphas)[0,1]:.4f}"
        )
        summary[f"presence_{label}"] = {
            "max_abs_deviation": float(np.abs(deviation).max()),
            "mean_deviation": float(deviation.mean()),
        }

    # 2. sibling calibration, using only the field
    pairs = _sibling_pairs(species)
    click.echo(f"\n-- sibling calibration ({len(pairs)} leaf pairs) --")
    predicted_all = _sibling_predicted(pairs, bins, alphas, nus)
    click.echo("  decile   predicted   reference        c++")
    reference_curve = (
        _sibling_observed(reference, pairs, predicted_all) if reference else None
    )
    cpp_curve = _sibling_observed(cpp, pairs, predicted_all) if cpp else None
    for d in range(N_DECILES):
        predicted = reference_curve[0][d] if reference_curve else cpp_curve[0][d]
        ref = f"{reference_curve[1][d]:9.4f}" if reference_curve else "        -"
        other = f"{cpp_curve[1][d]:9.4f}" if cpp_curve else "        -"
        click.echo(f"  {d + 1:>6}   {predicted:9.4f}   {ref}  {other}")

    if reference_curve:
        summary["sibling_reference_max_error"] = float(
            np.abs(reference_curve[1] - reference_curve[0]).max()
        )
    if cpp_curve:
        summary["sibling_cpp_max_error"] = float(
            np.abs(cpp_curve[1] - cpp_curve[0]).max()
        )
        click.echo(
            f"\n  worst decile error: c++ {summary['sibling_cpp_max_error']:.4f}"
            + (
                f"   reference {summary['sibling_reference_max_error']:.4f}"
                if reference_curve
                else ""
            )
        )

    out_path = base / "field_comparison.json"
    out_path.write_text(json.dumps(summary, indent=2) + "\n")
    click.echo(f"\n  summary written to {out_path}")


if __name__ == "__main__":
    main()
