"""Plot inference results against true (simulated) parameter values."""

import pathlib
import re
from dataclasses import dataclass

import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import PowerNorm
from dictances import cosine as cosine_dist
from dictances import mae as dict_mae
from dictances import mse as dict_mse
from dictances import pearson as dict_pearson
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

# the parameter types that get a panel, in the order the panels are laid out;
# anything classifying outside this list is never drawn
PLOT_ORDER = (
    "gamma",
    "epsilon",
    "epsilon_simple_model",
    "alpha",
    "log_nu",
    "mean_log_nu",
    "var_log_nu",
    "branch_lengths",
)

# The names below are what the C++ binary writes when it is run the way
# s_balanced_255_m_grass_200 was; other scenarios name their files differently,
# so each one is a CLI option and these are only its default.
DEFAULT_RUN_PREFIX = "acol"
DEFAULT_TRUE_VALUES = ("acol_input_simulated.txt", "acol_simulated.txt")
DEFAULT_TRUE_BRANCH_LENGTHS = "acol_{tree}_simulated.txt"
DEFAULT_TRUE_Y = "acol_simulated_Y.txt"

STATE_POSTERIORS_SUFFIX = "_statePosteriors.txt"


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


@dataclass(frozen=True)
class ResolvedInputs:
    """Every file the run reads, resolved and checked before anything renders."""

    mean_var: pathlib.Path
    trace: pathlib.Path
    y_posterior: pathlib.Path
    true_values: tuple[pathlib.Path, ...]
    true_y: pathlib.Path
    true_scalars: dict[str, float]
    # (tree, its statePosteriors file, its true-bin file), trees without a
    # true-bin file already dropped
    branch_truth: tuple[tuple[str, pathlib.Path, pathlib.Path], ...]


def _is_explicit(ctx: click.Context, name: str) -> bool:
    """Whether the option was typed on the command line rather than defaulted."""
    return ctx.get_parameter_source(name) is click.core.ParameterSource.COMMANDLINE


def _require_existing(paths: list[pathlib.Path], option: str) -> None:
    """A path named on the command line must exist: a typo is a failure, not a skip."""
    missing = [p for p in paths if not p.exists()]
    if missing:
        raise click.ClickException(
            f"{option}: no such file: " + ", ".join(str(p) for p in missing)
        )


def _parse_true_scalars(
    ctx: click.Context, param: click.Parameter, value: tuple[str, ...]
) -> dict[str, float]:
    scalars: dict[str, float] = {}
    for item in value:
        name, sep, raw = item.partition("=")
        if not sep or not name:
            raise click.BadParameter(
                f"expected name=value, got {item!r}", ctx=ctx, param=param
            )
        try:
            scalars[name] = float(raw)
        except ValueError:
            raise click.BadParameter(
                f"{name}: {raw!r} is not a number", ctx=ctx, param=param
            ) from None
    return scalars


def _resolve_inputs(
    ctx: click.Context,
    base: pathlib.Path,
    run_prefix: str,
    truth_dir: str | None,
    true_values: tuple[str, ...],
    true_branch_lengths: str,
    true_y: str,
    true_scalars: dict[str, float],
) -> ResolvedInputs:
    """Resolve the filename options against the run and truth directories.

    Everything is checked here rather than at the point of use, so a mistyped
    path fails before the first figure is drawn. A path named on the command line
    must exist; a defaulted one keeps the behaviour it had when the name was
    hardcoded, which is mandatory for the inference output and skipped with a
    message for the rest.
    """
    truth_root = pathlib.Path(truth_dir) if truth_dir is not None else base.parent
    if _is_explicit(ctx, "truth_dir") and not truth_root.is_dir():
        raise click.ClickException(f"--truth-dir: no such directory: {truth_root}")

    def under_truth(name: str) -> pathlib.Path:
        path = pathlib.Path(name)
        return path if path.is_absolute() else truth_root / path

    mean_var = base / f"{run_prefix}_meanVar.txt"
    if not mean_var.exists():
        raise click.ClickException(
            f"No inference output at {mean_var}. Run inference first."
        )

    value_paths = [under_truth(name) for name in true_values]
    if _is_explicit(ctx, "true_values"):
        _require_existing(value_paths, "--true-values")

    y_path = under_truth(true_y)
    if _is_explicit(ctx, "true_y"):
        _require_existing([y_path], "--true-y")

    # the prefix may carry a directory part, but the tree name is read off the
    # file's own name, so strip it back to the stem before slicing
    stem = pathlib.PurePath(run_prefix).name
    found: list[tuple[str, pathlib.Path, pathlib.Path]] = []
    missing: list[pathlib.Path] = []
    for state_file in sorted(base.glob(f"{run_prefix}_*{STATE_POSTERIORS_SUFFIX}")):
        tree = state_file.name[len(stem) + 1 : -len(STATE_POSTERIORS_SUFFIX)]
        # a template without {tree} names one pooled file covering every tree
        truth = under_truth(true_branch_lengths.replace("{tree}", tree))
        if truth.exists():
            found.append((tree, state_file, truth))
        else:
            missing.append(truth)
    if missing and _is_explicit(ctx, "true_branch_lengths"):
        _require_existing(missing, "--true-branch-lengths")

    return ResolvedInputs(
        mean_var=mean_var,
        trace=base / f"{run_prefix}_trace.txt",
        y_posterior=base / f"{run_prefix}_Y_posterior.txt",
        true_values=tuple(value_paths),
        true_y=y_path,
        true_scalars=true_scalars,
        branch_truth=tuple(found),
    )


def _load_true_values(
    paths: tuple[pathlib.Path, ...],
    scalars: dict[str, float],
    known_names: set[str],
) -> pd.DataFrame:
    """The true value of every parameter, layered across files then overridden.

    Files are applied in the order given and the last one to carry a name wins,
    which is how the default pair keeps `acol_simulated.txt` ahead of the
    `acol_input_simulated.txt` fallback. A --true-scalar outranks every file: it
    exists to supply a truth that no file holds, or to correct one that does.
    """
    frames: list[pd.DataFrame] = []
    for path in paths:
        frame = _load_tsv(path)
        if frame is None:
            continue
        frame.columns = pd.Index(["name", "true_value"])
        frames.append(frame)

    if not frames:
        tried = ", ".join(str(p) for p in paths)
        raise click.ClickException(
            f"No true-value file found (tried: {tried}). Run simulation first."
        )

    true_df = pd.concat(frames).drop_duplicates("name", keep="last")
    true_df["true_value"] = pd.to_numeric(true_df["true_value"], errors="coerce")

    for name in scalars:
        # the flag's whole point is to put a value on a plot, so a name that
        # would silently vanish before reaching one is an error, not a no-op
        if name not in known_names:
            raise click.ClickException(
                f"--true-scalar {name}: not a parameter in the inference output."
            )
        if _param_type(name) not in PLOT_ORDER:
            raise click.ClickException(
                f"--true-scalar {name}: classifies as "
                f"{_param_type(name)!r}, which is never plotted."
            )

    if scalars:
        injected = pd.DataFrame(
            {"name": list(scalars), "true_value": list(scalars.values())}
        )
        true_df = pd.concat([true_df, injected]).drop_duplicates("name", keep="last")

    return true_df.reset_index(drop=True)


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
    branch_truth: tuple[tuple[str, pathlib.Path, pathlib.Path], ...],
) -> tuple[pd.DataFrame, np.ndarray] | None:
    """Per-branch true bin and posterior over bins, pooled across all trees.

    Branch lengths are discrete, so inference reports a distribution over bins in
    the per-tree statePosteriors file rather than a mean/variance in the meanVar
    file; the true bins come from whichever file --true-branch-lengths resolved to
    for that tree, which may be one pooled file shared by all of them.
    """
    names: list[str] = []
    true_bins: list[int] = []
    posteriors: list[np.ndarray] = []

    for _tree, state_file, truth_file in branch_truth:
        true_df = _load_tsv(truth_file)
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
    return dict(zip(df["position"], df["Y_state"].astype(int)))


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
    ax.set_ylabel("True Y")
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
@click.option(
    "--run-prefix",
    default=DEFAULT_RUN_PREFIX,
    show_default=True,
    help="Filename stem of the inference output in <scenario_dir>, i.e. whatever "
    "was passed to the binary's --out.",
)
@click.option(
    "--truth-dir",
    type=click.Path(),
    default=None,
    help="Directory holding the simulation output that the truth options name "
    "(default: <scenario_dir>/..).",
)
@click.option(
    "--true-values",
    multiple=True,
    default=DEFAULT_TRUE_VALUES,
    show_default=True,
    help="True parameter values, under --truth-dir. Repeatable: files are layered "
    "in the order given and the last one to carry a name wins.",
)
@click.option(
    "--true-branch-lengths",
    default=DEFAULT_TRUE_BRANCH_LENGTHS,
    show_default=True,
    help="True branch-length bins, under --truth-dir. {tree} expands to each tree "
    "found; a name without it is one pooled file shared by every tree.",
)
@click.option(
    "--true-y",
    default=DEFAULT_TRUE_Y,
    show_default=True,
    help="True simulated Y, under --truth-dir.",
)
@click.option(
    "--true-scalar",
    "true_scalars",
    multiple=True,
    callback=_parse_true_scalars,
    metavar="NAME=VALUE",
    help="True value for a parameter no truth file holds. Repeatable; outranks "
    "the files. Errors if NAME is not an inferred parameter that gets a panel.",
)
@click.pass_context
def main(
    ctx: click.Context,
    scenario_dir: str,
    out: str | None,
    show: bool,
    y_threshold: float,
    run_prefix: str,
    truth_dir: str | None,
    true_values: tuple[str, ...],
    true_branch_lengths: str,
    true_y: str,
    true_scalars: dict[str, float],
) -> None:
    base = pathlib.Path(scenario_dir)
    inputs = _resolve_inputs(
        ctx,
        base,
        run_prefix,
        truth_dir,
        true_values,
        true_branch_lengths,
        true_y,
        true_scalars,
    )
    out_dir = pathlib.Path(out) if out else base / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- load inferred values (_resolve_inputs already checked the file is there) ---
    inferred_df = _load_tsv(inputs.mean_var)
    assert inferred_df is not None
    inferred_df.columns = pd.Index(["name", "post_mean", "post_var"])
    inferred_df["post_sd"] = np.sqrt(inferred_df["post_var"].clip(lower=0))

    # --- load the MCMC trace (posterior samples of the scalar parameters) ---
    trace_df = _load_tsv(inputs.trace)
    if trace_df is None:
        click.echo(
            f"No trace at {inputs.trace}: "
            "falling back to scatter plots for the scalar parameters."
        )

    true_df = _load_true_values(
        inputs.true_values, inputs.true_scalars, set(inferred_df["name"])
    )

    # --- merge and classify ---
    merged = inferred_df.merge(true_df, on="name", how="inner")
    merged["ptype"] = merged["name"].apply(_param_type)

    types_present = [pt for pt in PLOT_ORDER if (merged["ptype"] == pt).any()]

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
    sim_y = _load_y_file(inputs.true_y)
    sim_state = _load_y_state(inputs.true_y)
    post_y = _load_y_file(inputs.y_posterior)
    if sim_y is not None and sim_state is not None and post_y is not None:
        if sim_y and set(sim_y.values()) <= {0.0, 1.0}:
            click.echo(
                f"Note: fraction_of_one in {inputs.true_y} is binary, so the "
                "distance metrics restate the classification metrics beside them."
            )
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
            f"Skipping Y comparison: missing {inputs.true_y} or {inputs.y_posterior}."
        )

    # --- branch length posteriors vs true bins ---
    branch_data = _load_branch_lengths(inputs.branch_truth)
    if branch_data is None:
        click.echo(
            f"Skipping branch lengths: no {run_prefix}_<tree>"
            f"{STATE_POSTERIORS_SUFFIX} in {base} matched a true-bin file."
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
