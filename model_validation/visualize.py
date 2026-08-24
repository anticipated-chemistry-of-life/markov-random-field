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
from matplotlib.colors import PowerNorm
from scipy.stats import gaussian_kde
from sklearn.metrics import confusion_matrix, matthews_corrcoef


def _load_tsv(path: pathlib.Path) -> pd.DataFrame | None:
    if path.exists():
        return pd.read_csv(path, sep="\t")
    return None


DIM_COLORS = {"species": "#2196F3", "molecules": "#FF9800", "other": "#9E9E9E"}

# scalar parameters: few enough per run that the posterior itself is worth
# showing, so these are drawn as a trace KDE instead of a true-vs-inferred scatter
KDE_PTYPES = {"gamma", "epsilon", "epsilon_simple_model", "mean_log_nu", "var_log_nu"}


def _param_type(name: str) -> str:
    if re.match(r"^gamma_", name):
        return "gamma"
    # checked before "epsilon" so the exact-match test below cannot shadow it
    if name == "epsilon_simple_model":
        return "epsilon_simple_model"
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
    """Scatter: inferred posterior mean ± 2 SD on x, true value on y, colored by dimension."""
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
            inferred_means[mask],
            true_vals[mask],
            xerr=2 * inferred_sds[mask],
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
    ax.set_xlabel("Inferred posterior mean")
    ax.set_ylabel("True value")
    ax.set_title(f"{title} (n={len(true_vals)})")
    ax.legend(fontsize=8)


def _kde_posterior(
    ax: plt.Axes,
    names: list[str],
    true_by_name: dict[str, float],
    trace: pd.DataFrame,
    title: str,
    n_grid: int = 512,
) -> None:
    """Posterior KDE of the MCMC trace with a vertical line at the true value."""
    n_samples = 0
    for name in names:
        color = DIM_COLORS[_tree_dim(name)]
        samples = trace[name].to_numpy(dtype=float)
        samples = samples[np.isfinite(samples)]
        true_val = true_by_name[name]
        n_samples = max(n_samples, len(samples))

        if len(samples) < 2 or np.ptp(samples) == 0.0:
            # parameter was held fixed (update = 0): no density to estimate
            ax.axvline(
                samples[0] if len(samples) else np.nan,
                color=color,
                linewidth=2,
                label=f"{name} (fixed)",
            )
        else:
            kde = gaussian_kde(samples)
            lo = min(samples.min(), true_val)
            hi = max(samples.max(), true_val)
            margin = (hi - lo) * 0.1 or 0.1
            grid = np.linspace(lo - margin, hi + margin, n_grid)
            density = kde(grid)
            ax.plot(
                grid, density, color=color, linewidth=1.8, label=f"{name} posterior"
            )
            ax.fill_between(grid, density, color=color, alpha=0.2)

        ax.axvline(
            true_val,
            color=color,
            linestyle="--",
            linewidth=1.5,
            label=f"{name} true = {true_val:.4g}",
        )

    ax.set_xlabel("Value")
    ax.set_ylabel("Posterior density")
    ax.set_ylim(bottom=0)
    ax.set_title(f"{title} (n={n_samples} samples)")
    ax.legend(fontsize=8)


def _hpd_inclusion_level(probs: np.ndarray, true_bins: np.ndarray) -> np.ndarray:
    """Smallest HPD credible mass that already contains the true bin, per branch.

    Bins are added to the credible set in decreasing posterior order; the returned value
    is the cumulative mass at which the true bin joins. A calibrated posterior spreads
    these uniformly over [0, 1], so the share below any level q is the coverage of the
    q credible set.
    """
    order = np.argsort(-probs, axis=1, kind="stable")
    cumulative = np.cumsum(np.take_along_axis(probs, order, axis=1), axis=1)
    rank = np.argmax(order == true_bins[:, None], axis=1)
    return cumulative[np.arange(len(true_bins)), rank]


def _load_branch_lengths(
    base: pathlib.Path,
) -> tuple[pd.DataFrame, np.ndarray] | None:
    """Per-branch true bin and posterior over bins, pooled across all trees.

    Branch lengths are discrete, so inference reports a distribution over bins in
    acol_<tree>_statePosteriors.txt rather than a mean/variance in acol_meanVar.txt;
    the true bins likewise live in per-tree acol_<tree>_simulated.txt files.
    """
    names: list[str] = []
    true_bins: list[int] = []
    posteriors: list[np.ndarray] = []

    for state_file in sorted(base.glob("acol_*_statePosteriors.txt")):
        tree = state_file.name[len("acol_") : -len("_statePosteriors.txt")]
        true_df = _load_tsv(base.parent / f"acol_{tree}_simulated.txt")
        if true_df is None:
            continue
        true_df.columns = pd.Index(["name", "true_bin"])
        state_df = pd.read_csv(state_file, sep="\t")
        state_cols = [c for c in state_df.columns if c != "name"]

        merged = state_df.merge(true_df, on="name", how="inner")
        merged = merged[merged["name"].apply(_param_type) == "branch_lengths"]
        if merged.empty:
            continue
        names += merged["name"].tolist()
        true_bins += merged["true_bin"].astype(int).tolist()
        posteriors.append(merged[state_cols].to_numpy(dtype=float))

    if not names:
        return None
    bin_counts = {p.shape[1] for p in posteriors}
    if len(bin_counts) > 1:
        raise click.ClickException(
            f"Trees disagree on the number of branch-length bins: {sorted(bin_counts)}."
        )

    probs = np.vstack(posteriors)
    df = pd.DataFrame({"name": names, "true_bin": true_bins})

    n_bins = probs.shape[1]
    valid = ((df["true_bin"] >= 0) & (df["true_bin"] < n_bins)).to_numpy()
    if not valid.all():
        click.echo(
            f"Dropping {(~valid).sum()} branch(es) with a true bin outside [0, {n_bins})."
        )
        probs, df = probs[valid], df[valid].reset_index(drop=True)
        if df.empty:
            return None

    bins = np.arange(n_bins, dtype=float)
    post_mean = probs @ bins
    df["dim"] = [_tree_dim(n) for n in df["name"]]
    df["post_mean"] = post_mean
    df["post_sd"] = np.sqrt(np.clip(probs @ bins**2 - post_mean**2, 0.0, None))
    df["post_mode"] = probs.argmax(axis=1)
    df["prob_of_true"] = probs[np.arange(len(df)), df["true_bin"].to_numpy()]
    df["hpd_level"] = _hpd_inclusion_level(probs, df["true_bin"].to_numpy())
    return df, probs


def _compute_branch_metrics(df: pd.DataFrame) -> dict[str, float]:
    err = df["post_mean"] - df["true_bin"]
    return {
        "Branches": float(len(df)),
        "RMSE (posterior mean)": float(np.sqrt((err**2).mean())),
        "MAE (posterior mean)": float(err.abs().mean()),
        "Bias (inferred - true)": float(err.mean()),
        "Pearson r": float(df["post_mean"].corr(df["true_bin"].astype(float))),
        "Mode exact-match rate": float((df["post_mode"] == df["true_bin"]).mean()),
        "Mean P(true bin)": float(df["prob_of_true"].mean()),
        "Coverage of 50% HPD": float((df["hpd_level"] <= 0.5).mean()),
        "Coverage of 95% HPD": float((df["hpd_level"] <= 0.95).mean()),
    }


def _branch_metric_sections(df: pd.DataFrame) -> dict[str, dict[str, float]]:
    """Pooled metrics, plus a per-dimension breakdown once more than one tree is present."""
    sections = {"All branches": _compute_branch_metrics(df)}
    dims = [d for d in DIM_COLORS if (df["dim"] == d).any()]
    if len(dims) > 1:
        for dim in dims:
            sections[dim.title()] = _compute_branch_metrics(df[df["dim"] == dim])
    return sections


def _branch_posterior_heatmap(
    ax: plt.Axes,
    probs: np.ndarray,
    true_bins: np.ndarray,
    n_max_rows: int = 400,
) -> None:
    """Posterior mass per bin, one row per branch, sorted by the true bin.

    A well-inferred run puts the bright band on top of the red true-bin line.
    """
    if len(true_bins) > n_max_rows:
        rng = np.random.default_rng(0)
        idx = np.sort(rng.choice(len(true_bins), n_max_rows, replace=False))
        probs, true_bins = probs[idx], true_bins[idx]

    order = np.argsort(true_bins, kind="stable")
    probs, true_bins = probs[order], true_bins[order]
    n_rows, n_bins = probs.shape

    im = ax.imshow(
        probs,
        aspect="auto",
        origin="lower",
        cmap="viridis",
        interpolation="nearest",
        # square-root scale: posteriors are spiky, and a linear scale leaves every
        # branch but the most concentrated one indistinguishable from zero
        norm=PowerNorm(gamma=0.5, vmin=0.0, vmax=probs.max()),
        extent=(-0.5, n_bins - 0.5, -0.5, n_rows - 0.5),
    )
    ax.plot(
        true_bins, np.arange(n_rows), color="#F44336", linewidth=1.2, label="true bin"
    )
    ax.set_xlabel("Branch length bin")
    ax.set_ylabel("Branch (sorted by true bin)")
    ax.set_title(f"Posterior over bins (n={n_rows})")
    ax.legend(fontsize=8, loc="lower right")
    ax.figure.colorbar(im, ax=ax, label="Posterior probability")


def _branch_coverage_curve(
    ax: plt.Axes, hpd_levels: np.ndarray, n_grid: int = 101
) -> None:
    """Share of branches whose true bin falls inside the nominal HPD credible set."""
    nominal = np.linspace(0.0, 1.0, n_grid)
    empirical = (hpd_levels[None, :] <= nominal[:, None]).mean(axis=1)

    ax.plot(
        nominal, empirical, color="#1565C0", linewidth=2, label="empirical coverage"
    )
    ax.plot([0, 1], [0, 1], "k--", linewidth=1, label="perfect calibration")
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("Nominal HPD credible level")
    ax.set_ylabel("Fraction of branches covered")
    ax.set_title(f"Branch length calibration (n={len(hpd_levels)})")
    ax.legend(fontsize=8, loc="lower right")


def _load_y_file(path: pathlib.Path) -> dict[int, float] | None:
    if not path.exists():
        return None
    df = pd.read_csv(path, sep="\t")
    return dict(zip(df["position"], df["fraction_of_one"].astype(float)))


def _load_y_state(path: pathlib.Path) -> dict[int, int] | None:
    """True binary Y per position (the simulated latent state)."""
    if not path.exists():
        return None
    df = pd.read_csv(path, sep="\t")
    return dict(zip(df["position"], (df["fraction_of_one"].values > 0.5).astype(int)))


def _compute_y_metrics(true_dict: dict, pred_dict: dict) -> dict[str, float]:
    return {
        "MAE": dict_mae(true_dict, pred_dict),
        "MSE": dict_mse(true_dict, pred_dict),
        "Cosine similarity": 1.0 - cosine_dist(true_dict, pred_dict),
    }


def _binarize_y(
    true_state: dict[int, int],
    pred_prob: dict[int, float],
    threshold: float,
) -> tuple[np.ndarray, np.ndarray]:
    """True binary Y and the thresholded prediction, aligned on the true positions."""
    keys = list(true_state)
    y_true = np.array([true_state[k] for k in keys], dtype=int)
    probs = np.array([pred_prob.get(k, 0.0) for k in keys], dtype=float)
    return y_true, (probs >= threshold).astype(int)


def _compute_y_classification_metrics(
    y_true: np.ndarray, y_pred: np.ndarray
) -> dict[str, float]:
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    precision = tp / (tp + fp) if tp + fp else float("nan")
    recall = tp / (tp + fn) if tp + fn else float("nan")
    denom = precision + recall
    return {
        "MCC": matthews_corrcoef(y_true, y_pred),
        "Accuracy": (tp + tn) / len(y_true),
        "Precision": precision,
        "Recall (sensitivity)": recall,
        "Specificity": tn / (tn + fp) if tn + fp else float("nan"),
        "F1": 2 * precision * recall / denom if denom else float("nan"),
    }


def _plot_metrics_table(
    ax: plt.Axes, sections: dict[str, dict[str, float]], title: str
) -> None:
    ax.axis("off")
    lines: list[str] = []
    for section, metrics in sections.items():
        if lines:
            lines.append("")
        lines += [section, "-" * 43]
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
    ax.set_title(title)


def _plot_confusion_matrix(ax: plt.Axes, cm: np.ndarray, threshold: float) -> None:
    """2x2 confusion matrix heatmap annotated with counts and share of all positions."""
    total = cm.sum()
    ax.imshow(cm, cmap="Blues")
    for i in range(2):
        for j in range(2):
            count = int(cm[i, j])
            ax.text(
                j,
                i,
                f"{count:,}\n{100 * count / total:.1f}%",
                ha="center",
                va="center",
                color="white" if count > cm.max() / 2 else "black",
                fontsize=11,
            )
    ax.set_xticks([0, 1], ["Y=0", "Y=1"])
    ax.set_yticks([0, 1], ["Y=0", "Y=1"])
    ax.set_xlabel(f"Predicted (P(Y=1) >= {threshold:g})")
    ax.set_ylabel("True Y")
    ax.set_title(f"Confusion matrix (n={total:,})")


def _logistic_y(
    ax: plt.Axes,
    true_state: dict[int, int],
    pred_prob: dict[int, float],
    n_max_points: int = 30000,
    n_bins: int = 10,
) -> None:
    """Logistic-regression-style plot: predicted P(Y=1) vs true binary Y.

    Raw points sit at y=0/1 (jittered for visibility); the binned empirical
    P(Y=1) overlays a calibration curve that should track the diagonal.
    """
    keys = list(true_state)
    y_true = np.array([true_state[k] for k in keys], dtype=float)
    x_pred = np.array([pred_prob.get(k, 0.0) for k in keys], dtype=float)

    rng = np.random.default_rng(0)
    if len(keys) > n_max_points:
        idx = rng.choice(len(keys), n_max_points, replace=False)
        x_pred, y_true = x_pred[idx], y_true[idx]

    jitter = rng.uniform(-0.03, 0.03, size=len(y_true))
    ax.scatter(x_pred, y_true + jitter, alpha=0.2, s=6, color="#4CAF50", label="data")

    # binned empirical P(Y=1): mean true Y within each predicted-probability bin
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    which = np.clip(np.digitize(x_pred, edges) - 1, 0, n_bins - 1)
    centers, means = [], []
    for b in range(n_bins):
        m = which == b
        if m.any():
            centers.append(x_pred[m].mean())
            means.append(y_true[m].mean())
    ax.plot(
        centers,
        means,
        "o-",
        color="#1565C0",
        lw=2,
        ms=5,
        label="binned empirical P(Y=1)",
    )

    ax.plot([0, 1], [0, 1], "k--", lw=1, label="perfect calibration")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.12, 1.12)
    ax.set_xlabel("Predicted P(Y=1)")
    ax.set_ylabel("True Y (acol_simulated_Y.txt)")
    ax.set_title(f"Y calibration (n={len(y_true)})")
    ax.legend(fontsize=8, loc="center right")


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
@click.option(
    "--y-threshold",
    type=click.FloatRange(0.0, 1.0),
    default=0.5,
    show_default=True,
    help="P(Y=1) cutoff used to binarize the posterior for the confusion matrix / MCC.",
)
def main(scenario_dir: str, out: str | None, show: bool, y_threshold: float) -> None:
    base = pathlib.Path(scenario_dir)
    out_dir = pathlib.Path(out) if out else base / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- load inferred values ---
    inferred_df = _load_tsv(base / "acol_meanVar.txt")
    if inferred_df is None:
        raise click.ClickException(
            f"No inference output at {base / 'acol_meanVar.txt'}. Run inference first."
        )
    inferred_df.columns = pd.Index(["name", "post_mean", "post_var"])
    inferred_df["post_sd"] = np.sqrt(inferred_df["post_var"].clip(lower=0))

    # --- load the MCMC trace (posterior samples of the scalar parameters) ---
    trace_df = _load_tsv(base / "acol_trace.txt")
    if trace_df is None:
        click.echo(
            f"No trace at {base / 'acol_trace.txt'}: "
            "falling back to scatter plots for the scalar parameters."
        )

    # --- load true values: acol_simulated.txt is primary, input file is fallback ---
    true_df = _load_tsv(base.parent / "acol_simulated.txt")
    fallback_df = _load_tsv(base.parent / "acol_input_simulated.txt")
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
        "epsilon_simple_model",
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
        title = ptype.replace("_", " ").title()
        traced = (
            [n for n in subset["name"] if n in trace_df.columns]
            if trace_df is not None and ptype in KDE_PTYPES
            else []
        )
        if traced:
            _kde_posterior(
                axes_flat[ax_idx],
                traced,
                dict(zip(subset["name"], subset["true_value"])),
                trace_df,
                title=title,
            )
            continue
        _scatter_true_vs_inferred(
            axes_flat[ax_idx],
            subset["name"].tolist(),
            subset["true_value"].to_numpy(dtype=float),
            subset["post_mean"].to_numpy(dtype=float),
            subset["post_sd"].to_numpy(dtype=float),
            title=title,
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
    sim_y = _load_y_file(base.parent / "acol_simulated_Y.txt")
    sim_state = _load_y_state(base.parent / "acol_simulated_Y.txt")
    post_y = _load_y_file(base / "acol_Y_posterior.txt")
    if sim_y is not None and sim_state is not None and post_y is not None:
        y_true, y_pred = _binarize_y(sim_state, post_y, y_threshold)
        cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
        sections = {
            "Y distribution distances (dictances)": _compute_y_metrics(sim_y, post_y),
            f"Y classification (cutoff {y_threshold:g})": (
                _compute_y_classification_metrics(y_true, y_pred)
            ),
        }
        fig_y, (ax_scatter, ax_cm, ax_metrics) = plt.subplots(
            1, 3, figsize=(18, 5), gridspec_kw={"width_ratios": [1.2, 0.85, 1.1]}
        )
        _logistic_y(ax_scatter, sim_state, post_y)
        _plot_confusion_matrix(ax_cm, cm, y_threshold)
        _plot_metrics_table(ax_metrics, sections, title="Y metrics")
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

    # --- branch length posteriors vs true bins ---
    branch_data = _load_branch_lengths(base)
    if branch_data is None:
        click.echo(
            "Skipping branch lengths: no acol_<tree>_statePosteriors.txt in "
            f"{base} matched an acol_<tree>_simulated.txt in {base.parent}."
        )
    else:
        branches, branch_probs = branch_data
        fig_b, axes_b = plt.subplots(2, 2, figsize=(14, 11))
        _scatter_true_vs_inferred(
            axes_b[0, 0],
            branches["name"].tolist(),
            branches["true_bin"].to_numpy(dtype=float),
            branches["post_mean"].to_numpy(dtype=float),
            branches["post_sd"].to_numpy(dtype=float),
            title="Branch lengths",
        )
        axes_b[0, 0].set_xlabel("Inferred posterior mean (bin)")
        axes_b[0, 0].set_ylabel("True bin")
        _branch_posterior_heatmap(
            axes_b[0, 1], branch_probs, branches["true_bin"].to_numpy()
        )
        _branch_coverage_curve(axes_b[1, 0], branches["hpd_level"].to_numpy())
        _plot_metrics_table(
            axes_b[1, 1],
            _branch_metric_sections(branches),
            title="Branch length metrics",
        )
        fig_b.suptitle(f"Branch lengths — {base.name}", fontsize=13)
        fig_b.tight_layout()
        out_file_b = out_dir / "branch_lengths.pdf"
        fig_b.savefig(out_file_b, bbox_inches="tight")
        click.echo(f"Saved: {out_file_b}")
        if show:
            plt.show()
        plt.close(fig_b)


if __name__ == "__main__":
    main()
