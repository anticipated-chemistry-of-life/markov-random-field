"""Plot inference results against true (simulated) parameter values."""

import pathlib
import re

import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _load_tsv(path: pathlib.Path) -> pd.DataFrame | None:
    if path.exists():
        return pd.read_csv(path, sep="\t")
    return None


def _param_type(name: str) -> str:
    """Classify a parameter name into a broad category."""
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


def _scatter_true_vs_inferred(
    ax: plt.Axes,
    true_vals: np.ndarray,
    inferred_means: np.ndarray,
    inferred_sds: np.ndarray,
    title: str,
    n_max_points: int = 500,
) -> None:
    """Scatter plot: true value on x, inferred posterior mean on y, ±2 SD error bars."""
    if len(true_vals) > n_max_points:
        rng = np.random.default_rng(0)
        idx = rng.choice(len(true_vals), n_max_points, replace=False)
        true_vals, inferred_means, inferred_sds = (
            true_vals[idx],
            inferred_means[idx],
            inferred_sds[idx],
        )

    ax.errorbar(
        true_vals,
        inferred_means,
        yerr=2 * inferred_sds,
        fmt="o",
        markersize=3,
        alpha=0.5,
        linewidth=0.5,
        capsize=2,
    )
    lo = min(true_vals.min(), inferred_means.min())
    hi = max(true_vals.max(), inferred_means.max())
    margin = (hi - lo) * 0.05 or 0.1
    diag = np.array([lo - margin, hi + margin])
    ax.plot(diag, diag, "k--", linewidth=1, label="y = x")
    ax.set_xlabel("True value")
    ax.set_ylabel("Inferred posterior mean")
    ax.set_title(title)
    ax.legend(fontsize=8)


def _bar_inferred_only(
    ax: plt.Axes,
    names: list[str],
    inferred_means: np.ndarray,
    inferred_sds: np.ndarray,
    title: str,
) -> None:
    """Bar chart for scalar parameters that have no matching true value."""
    x = np.arange(len(names))
    ax.bar(x, inferred_means, color="steelblue", alpha=0.7)
    ax.errorbar(x, inferred_means, yerr=2 * inferred_sds, fmt="none", color="black", capsize=4)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=15, ha="right")
    ax.set_ylabel("Posterior mean ± 2 SD")
    ax.set_title(title)


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

    # --- load data ---
    inferred_df = _load_tsv(base / "test_out" / "acol_meanVar.txt")
    if inferred_df is None:
        raise click.ClickException(
            f"No inference output found at {base / 'test_out' / 'acol_meanVar.txt'}.\n"
            "Run the inference first."
        )
    inferred_df.columns = pd.Index(["name", "post_mean", "post_var"])
    inferred_df["post_sd"] = np.sqrt(inferred_df["post_var"].clip(lower=0))

    # Collect true values from multiple sources (simulated output takes priority)
    true_df = _load_tsv(base / "acol_input_simulated.txt")
    for extra in ["acol_simulated.txt", "acol_molecules_simulated.txt", "acol_species_simulated.txt"]:
        extra_df = _load_tsv(base / extra)
        if extra_df is not None and true_df is not None:
            true_df = pd.concat([true_df, extra_df]).drop_duplicates("name", keep="last")
        elif extra_df is not None:
            true_df = extra_df

    if true_df is not None:
        true_df.columns = pd.Index(["name", "true_value"])
        true_df["true_value"] = pd.to_numeric(true_df["true_value"], errors="coerce")

    # Tag each inferred parameter with its type
    inferred_df["ptype"] = inferred_df["name"].apply(_param_type)

    # --- build matched and unmatched sets ---
    if true_df is not None:
        merged = inferred_df.merge(true_df, on="name", how="left")
    else:
        merged = inferred_df.copy()
        merged["true_value"] = float("nan")

    matched = merged.dropna(subset=["true_value"])
    unmatched_scalars = merged[merged["true_value"].isna() & merged["ptype"].isin(["gamma", "epsilon"])]

    # --- plot matched parameters by type ---
    comparable_types = [pt for pt in ["alpha", "log_nu", "branch_lengths", "mean_log_nu", "var_log_nu"]
                        if (matched["ptype"] == pt).any()]

    n_plots = len(comparable_types) + (1 if not unmatched_scalars.empty else 0)
    if n_plots == 0:
        click.echo("Nothing to plot: no parameters with known true values found.")
        return

    ncols = min(3, n_plots)
    nrows = (n_plots + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5 * nrows), squeeze=False)
    axes_flat = axes.flatten()

    ax_idx = 0
    for ptype in comparable_types:
        subset = matched[matched["ptype"] == ptype]
        _scatter_true_vs_inferred(
            axes_flat[ax_idx],
            subset["true_value"].to_numpy(dtype=float),
            subset["post_mean"].to_numpy(dtype=float),
            subset["post_sd"].to_numpy(dtype=float),
            title=ptype.replace("_", " ").title(),
        )
        ax_idx += 1

    if not unmatched_scalars.empty:
        _bar_inferred_only(
            axes_flat[ax_idx],
            unmatched_scalars["name"].tolist(),
            unmatched_scalars["post_mean"].to_numpy(dtype=float),
            unmatched_scalars["post_sd"].to_numpy(dtype=float),
            title="Inferred scalars (gamma / epsilon)",
        )
        ax_idx += 1

    for ax in axes_flat[ax_idx:]:
        ax.set_visible(False)

    fig.suptitle(f"Inference results — {base.name}", fontsize=13, y=1.01)
    fig.tight_layout()

    out_file = out_dir / "inference_summary.pdf"
    fig.savefig(out_file, bbox_inches="tight")
    click.echo(f"Saved: {out_file}")

    if show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
