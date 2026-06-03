"""Plot inference results against true (simulated) parameter values."""

import pathlib
import re

import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from dictances import cosine as cosine_dist
from dictances import mae as dict_mae
from dictances import mse as dict_mse
from dictances import pearson as dict_pearson


def _load_tsv(path: pathlib.Path) -> pd.DataFrame | None:
    if path.exists():
        return pd.read_csv(path, sep="\t")
    return None


DIM_COLORS = {"species": "#2196F3", "molecules": "#FF9800", "other": "#9E9E9E"}


def _param_type(name: str) -> str:
    if re.match(r"^gamma_", name):
        return "gamma"
    if name == "epsilon":
        return "epsilon"
    if re.search(r"_alpha_", name):
        return "alpha"
    if re.search(r"_log_nu_", name) and not re.search(r"mean_log_nu|var_log_nu", name):
        return "log_nu"
    if re.search(r"_mean_log_nu$", name):
        return "mean_log_nu"
    if re.search(r"_var_log_nu$", name):
        return "var_log_nu"
    if re.search(r"_branch_lengths_", name):
        return "branch_lengths"
    return "other"


def _tree_dim(name: str) -> str:
    """Return which tree dimension a parameter belongs to."""
    if name.startswith("species_") or name == "gamma_species":
        return "species"
    if name.startswith("molecules_") or name == "gamma_molecules":
        return "molecules"
    return "other"


def _scatter_true_vs_inferred(
    ax: plt.Axes,
    names: list[str],
    true_vals: np.ndarray,
    inferred_means: np.ndarray,
    inferred_sds: np.ndarray,
    title: str,
    n_max_points: int = 500,
) -> None:
    """Scatter: true value on x, inferred posterior mean ± 2 SD on y, colored by dimension."""
    dims = np.array([_tree_dim(n) for n in names])

    if len(true_vals) > n_max_points:
        rng = np.random.default_rng(0)
        idx = rng.choice(len(true_vals), n_max_points, replace=False)
        true_vals, inferred_means, inferred_sds, dims = (
            true_vals[idx],
            inferred_means[idx],
            inferred_sds[idx],
            dims[idx],
        )

    plotted_dims: set[str] = set()
    for dim, color in DIM_COLORS.items():
        mask = dims == dim
        if not mask.any():
            continue
        label = dim if dim not in plotted_dims else None
        ax.errorbar(
            true_vals[mask],
            inferred_means[mask],
            yerr=2 * inferred_sds[mask],
            fmt="o",
            color=color,
            markersize=4,
            alpha=0.6,
            linewidth=0.5,
            capsize=2,
            label=label,
        )
        plotted_dims.add(dim)

    lo = min(true_vals.min(), inferred_means.min())
    hi = max(true_vals.max(), inferred_means.max())
    margin = (hi - lo) * 0.05 or 0.1
    diag = np.array([lo - margin, hi + margin])
    ax.plot(diag, diag, "k--", linewidth=1, label="y = x")
    ax.set_xlabel("True value")
    ax.set_ylabel("Inferred posterior mean")
    ax.set_title(f"{title} (n={len(true_vals)})")
    ax.legend(fontsize=8)


def _load_y_file(path: pathlib.Path) -> dict[int, float] | None:
    if not path.exists():
        return None
    df = pd.read_csv(path, sep="\t")
    return dict(zip(df["position"], df["fraction_of_one"].astype(float)))


def _compute_y_metrics(true_dict: dict, pred_dict: dict) -> dict[str, float]:
    return {
        "MAE": dict_mae(true_dict, pred_dict),
        "MSE": dict_mse(true_dict, pred_dict),
        "Pearson r": dict_pearson(true_dict, pred_dict),
        "Cosine similarity": 1.0 - cosine_dist(true_dict, pred_dict),
    }


def _plot_y_metrics(ax: plt.Axes, metrics: dict[str, float]) -> None:
    ax.axis("off")
    lines = [f"{'Metric':<30}{'Value':>12}", "-" * 43]
    for k, v in metrics.items():
        lines.append(f"{k:<30}{v:>12.4f}")
    ax.text(
        0.05,
        0.95,
        "\n".join(lines),
        transform=ax.transAxes,
        verticalalignment="top",
        fontfamily="monospace",
        fontsize=10,
    )
    ax.set_title("Y distribution distances (dictances)")


def _scatter_y(
    ax: plt.Axes,
    true_dict: dict[int, float],
    pred_dict: dict[int, float],
    n_max_points: int = 30000,
) -> None:
    all_keys = list(set(true_dict))
    true_vec = np.array([true_dict.get(k, 0.0) for k in all_keys])
    pred_vec = np.array([pred_dict.get(k, 0.0) for k in all_keys])

    if len(true_vec) > n_max_points:
        rng = np.random.default_rng(0)
        idx = rng.choice(len(true_vec), n_max_points, replace=False)
        true_vec, pred_vec = true_vec[idx], pred_vec[idx]

    ax.scatter(true_vec, pred_vec, alpha=0.3, s=6, color="#4CAF50")
    lo = min(true_vec.min(), pred_vec.min())
    hi = max(true_vec.max(), pred_vec.max())
    margin = (hi - lo) * 0.05 or 0.1
    diag = np.array([lo - margin, hi + margin])
    ax.plot(diag, diag, "k--", linewidth=1)
    ax.set_xlabel("True fraction_of_one")
    ax.set_ylabel("Predicted fraction_of_one")
    ax.set_title(f"Y probabilities (n={len(true_vec)})")


@click.command()
@click.argument("scenario_dir", type=click.Path(exists=True, file_okay=False))
@click.option(
    "--out",
    type=click.Path(),
    default=None,
    help="Directory to save plots (default: <scenario_dir>/plots).",
)
@click.option(
    "--show/--no-show",
    default=False,
    help="Display plots interactively after saving.",
)
def main(scenario_dir: str, out: str | None, show: bool) -> None:
    base = pathlib.Path(scenario_dir)
    out_dir = pathlib.Path(out) if out else base / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- load inferred values ---
    inferred_df = _load_tsv(base / "test_out" / "acol_meanVar.txt")
    if inferred_df is None:
        raise click.ClickException(
            f"No inference output at {base / 'test_out' / 'acol_meanVar.txt'}. "
            "Run inference first."
        )
    inferred_df.columns = pd.Index(["name", "post_mean", "post_var"])
    inferred_df["post_sd"] = np.sqrt(inferred_df["post_var"].clip(lower=0))

    # --- load true values: acol_simulated.txt is primary, input file is fallback ---
    true_df = _load_tsv(base / "acol_simulated.txt")
    fallback_df = _load_tsv(base / "acol_input_simulated.txt")
    if true_df is None and fallback_df is None:
        raise click.ClickException("No true-value file found. Run simulation first.")
    if true_df is None:
        true_df = fallback_df
    elif fallback_df is not None:
        # merge: simulated takes priority, input fills in anything missing
        true_df = pd.concat([fallback_df, true_df]).drop_duplicates("name", keep="last")

    true_df.columns = pd.Index(["name", "true_value"])
    true_df["true_value"] = pd.to_numeric(true_df["true_value"], errors="coerce")

    # --- merge and classify ---
    merged = inferred_df.merge(true_df, on="name", how="inner")
    merged["ptype"] = merged["name"].apply(_param_type)

    plot_order = [
        "gamma",
        "epsilon",
        "alpha",
        "log_nu",
        "mean_log_nu",
        "var_log_nu",
        "branch_lengths",
    ]
    types_present = [pt for pt in plot_order if (merged["ptype"] == pt).any()]

    if not types_present:
        click.echo(
            "Nothing to plot: no parameters matched between inference and simulation output."
        )
        return

    ncols = min(3, len(types_present))
    nrows = (len(types_present) + ncols - 1) // ncols
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(6 * ncols, 5 * nrows), squeeze=False
    )
    axes_flat = axes.flatten()

    for ax_idx, ptype in enumerate(types_present):
        subset = merged[merged["ptype"] == ptype]
        _scatter_true_vs_inferred(
            axes_flat[ax_idx],
            subset["name"].tolist(),
            subset["true_value"].to_numpy(dtype=float),
            subset["post_mean"].to_numpy(dtype=float),
            subset["post_sd"].to_numpy(dtype=float),
            title=ptype.replace("_", " ").title(),
        )

    for ax in axes_flat[len(types_present) :]:
        ax.set_visible(False)

    fig.suptitle(f"Inference results — {base.name}", fontsize=13, y=1.01)
    fig.tight_layout()

    out_file = out_dir / "inference_summary.pdf"
    fig.savefig(out_file, bbox_inches="tight")
    click.echo(f"Saved: {out_file}")

    if show:
        plt.show()
    plt.close(fig)

    # --- Y distribution comparison ---
    sim_y = _load_y_file(base / "acol_simulated_Y.txt")
    post_y = _load_y_file(base / "test_out" / "acol_Y_posterior.txt")
    if sim_y is not None and post_y is not None:
        y_metrics = _compute_y_metrics(sim_y, post_y)
        fig_y, (ax_scatter, ax_metrics) = plt.subplots(1, 2, figsize=(12, 5))
        _scatter_y(ax_scatter, sim_y, post_y)
        _plot_y_metrics(ax_metrics, y_metrics)
        fig_y.suptitle(f"Y comparison — {base.name}", fontsize=13)
        fig_y.tight_layout()
        out_file_y = out_dir / "y_distribution_distances.pdf"
        fig_y.savefig(out_file_y, bbox_inches="tight")
        click.echo(f"Saved: {out_file_y}")
        if show:
            plt.show()
        plt.close(fig_y)
    else:
        click.echo(
            "Skipping Y comparison: missing acol_simulated_Y.txt or acol_Y_posterior.txt."
        )


if __name__ == "__main__":
    main()
