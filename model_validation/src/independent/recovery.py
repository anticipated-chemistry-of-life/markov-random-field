"""Scores for how well inference recovers the field and the detection parameters.

Pure functions over arrays and frames, so the validator script only reads files.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import binom, rankdata

from .mass_spec import UNKNOWN

NOMINAL_COVERAGE = 0.9
BAND_MASS = 0.95
CREDIBLE_MASS = 0.9


def coverage_band(n: int) -> tuple[float, float]:
    """The central 95% of Binomial(n, 0.9), as fractions of n."""
    tail = (1.0 - BAND_MASS) / 2.0
    low = binom.ppf(tail, n, NOMINAL_COVERAGE)
    high = binom.ppf(1.0 - tail, n, NOMINAL_COVERAGE)
    return float(low / n), float(high / n)


def coverage_passes(coverage: float, n: int) -> bool:
    low, high = coverage_band(n)
    return low <= coverage <= high


def is_covered(draws: np.ndarray, truth: np.ndarray) -> np.ndarray:
    """Whether each column's central 90% credible interval holds its truth.

    `draws` is `[sample, parameter]`.
    """
    tail = (1.0 - CREDIBLE_MASS) / 2.0
    lower, upper = np.quantile(draws, [tail, 1.0 - tail], axis=0)
    return (truth >= lower) & (truth <= upper)


@dataclass(frozen=True)
class CoverageRow:
    label: str
    n: int
    coverage: float
    passes: bool


# Inclusive run-count ranges; the last is open-ended.
RUN_STRATA = (("0", 0, 0), ("1", 1, 1), ("2", 2, 2), ("3-9", 3, 9), ("10+", 10, None))


def coverage_by_runs(
    draws: np.ndarray, truth: np.ndarray, runs_present: np.ndarray
) -> list[CoverageRow]:
    """Coverage per stratum of `runs_present`, the number of runs a molecule is present in."""
    covered = is_covered(draws, truth)
    rows = []
    for label, low, high in RUN_STRATA:
        mask = runs_present >= low
        if high is not None:
            mask &= runs_present <= high
        n = int(mask.sum())
        if n:
            coverage = float(covered[mask].mean())
            rows.append(CoverageRow(label, n, coverage, coverage_passes(coverage, n)))
    return rows


def auc(truth: np.ndarray, score: np.ndarray) -> float:
    """Area under the ROC curve; ties count half. NaN when only one class is present."""
    positive = np.asarray(truth) == 1
    n_positive = int(positive.sum())
    n_negative = positive.size - n_positive
    if n_positive == 0 or n_negative == 0:
        return float("nan")
    ranks = rankdata(score)
    return float(
        (ranks[positive].sum() - n_positive * (n_positive + 1) / 2)
        / (n_positive * n_negative)
    )


@dataclass(frozen=True)
class ReliabilityRow:
    low: float
    high: float
    n: int
    predicted: float
    """Mean predicted P(Y=1) in the bin."""
    observed: float
    """Fraction of the bin's cells that are truly 1."""


def reliability_table(
    truth: np.ndarray, predicted: np.ndarray, n_bins: int = 10
) -> list[ReliabilityRow]:
    """Predicted P(Y=1) against the true Y, in equal-width bins; empty bins are skipped."""
    # The top bin includes 1.0.
    bin_of = np.minimum((predicted * n_bins).astype(int), n_bins - 1)
    rows = []
    for b in range(n_bins):
        mask = bin_of == b
        if mask.any():
            rows.append(
                ReliabilityRow(
                    b / n_bins,
                    (b + 1) / n_bins,
                    int(mask.sum()),
                    float(predicted[mask].mean()),
                    float(truth[mask].mean()),
                )
            )
    return rows


@dataclass(frozen=True)
class FieldScore:
    n: int
    auc: float
    reliability: list[ReliabilityRow]


def field_scores(
    truth: pd.DataFrame, posterior: pd.DataFrame, species_with_runs: set[str]
) -> dict[str, FieldScore]:
    """AUC and a reliability table for P(Y=1), for species with and without runs.

    `truth` lists every cell with its `Y_state`. `posterior` is the C++ Y posterior,
    which lists only the cells that carry one: any other cell was never a one.
    """
    merged = truth[["species", "molecules", "Y_state"]].merge(
        posterior[["species", "molecules", "fraction_of_one"]],
        on=["species", "molecules"],
        how="left",
    )
    merged["fraction_of_one"] = merged["fraction_of_one"].fillna(0.0)
    has_runs = merged["species"].isin(species_with_runs).to_numpy()

    scores = {}
    for key, mask in (("with_runs", has_runs), ("without_runs", ~has_runs)):
        y = merged["Y_state"].to_numpy()[mask]
        p = merged["fraction_of_one"].to_numpy()[mask]
        scores[key] = FieldScore(int(mask.sum()), auc(y, p), reliability_table(y, p))
    return scores


def runs_present(
    truth: pd.DataFrame,
    species_with_runs: set[str],
    runs_per_species: int,
    molecules: list[str],
) -> np.ndarray:
    """For each of `molecules`, the number of runs whose species truly carries it."""
    carried = truth[(truth["Y_state"] == 1) & truth["species"].isin(species_with_runs)]
    species_count = carried.groupby("molecules").size()
    return (
        species_count.reindex(molecules, fill_value=0).to_numpy() * runs_per_species
    )


# --------------------------------------------------------------------------
# The assignment posterior
# --------------------------------------------------------------------------

FEATURE_KEY = ["species", "extract", "feature"]
SUM_TOLERANCE = 1e-9
"""The C++ writes each fraction as the shortest text that reads back exactly, so a
feature's fractions sum to 1 to rounding of the sum alone."""


class AssignmentFormatError(ValueError):
    """The assignment posterior does not match the MS data file."""


def _feature_label(key: tuple) -> str:
    species, extract, feature = key
    return f"feature {feature!r} (species {species!r}, extract {extract!r})"


def check_assignment_format(data: pd.DataFrame, posterior: pd.DataFrame) -> None:
    """Every (run, feature) of `data` is present with all its candidates and the unknown
    molecule, nothing else is, and each feature's fractions are in [0, 1] and sum to 1.

    `data` is the MS data file (one row per candidate); `posterior` the C++ file.
    Raises `AssignmentFormatError`, naming the first offending feature.
    """
    expected: dict[tuple, set[str]] = {}
    for *key, molecule in data[FEATURE_KEY + ["molecule"]].itertuples(index=False):
        expected.setdefault(tuple(key), set()).add(str(molecule))
    for candidates in expected.values():
        candidates.add(UNKNOWN)

    seen: dict[tuple, list[str]] = {}
    for *key, molecule in posterior[FEATURE_KEY + ["molecule"]].itertuples(index=False):
        seen.setdefault(tuple(key), []).append(str(molecule))

    for key in expected:
        if key not in seen:
            raise AssignmentFormatError(
                f"Assignment posterior has no rows for {_feature_label(key)}."
            )
    for key, molecules in seen.items():
        if key not in expected:
            raise AssignmentFormatError(
                f"Assignment posterior has rows for {_feature_label(key)}, "
                "which the MS data file does not have."
            )
        if len(molecules) != len(set(molecules)):
            raise AssignmentFormatError(
                f"Assignment posterior repeats a molecule for {_feature_label(key)}."
            )
        missing = expected[key] - set(molecules)
        if missing:
            raise AssignmentFormatError(
                f"Assignment posterior for {_feature_label(key)} lacks "
                f"{sorted(missing)}."
            )
        extra = set(molecules) - expected[key]
        if extra:
            raise AssignmentFormatError(
                f"Assignment posterior for {_feature_label(key)} has {sorted(extra)}, "
                "which are not its candidates."
            )

    fractions = posterior["fraction"].astype(float)
    for key, group in posterior.assign(fraction=fractions).groupby(FEATURE_KEY, sort=False):
        values = group["fraction"].to_numpy()
        if not ((values >= 0) & (values <= 1)).all():
            raise AssignmentFormatError(
                f"Assignment posterior for {_feature_label(key)} has a fraction "
                "outside [0, 1]."
            )
        total = float(values.sum())
        if abs(total - 1.0) > SUM_TOLERANCE:
            raise AssignmentFormatError(
                f"Assignment posterior for {_feature_label(key)} sums to {total!r}, "
                "not 1."
            )


@dataclass(frozen=True)
class AssignmentScore:
    n: int
    """Features scored."""
    accuracy: float
    """Fraction of features whose posterior mode is the true explanation."""
    mean_true: float
    """Mean posterior mass on the true explanation."""
    n_unknown_features: int
    mean_unknown_of_unknown_features: float
    """Mean posterior on the unknown molecule, over features that are truly unknown."""
    n_molecule_features: int
    mean_unknown_of_molecule_features: float
    """Mean posterior on the unknown molecule, over features of a detected molecule."""


def assignment_scores(truth: pd.DataFrame, posterior: pd.DataFrame) -> AssignmentScore:
    """How well the assignment posterior recovers the true assignment.

    Both frames are in the C++ layout (`species extract feature molecule fraction`),
    the truth with fraction 1 on the true explanation. The mode of a tie is the
    first row the file lists.
    """
    merged = truth.merge(
        posterior, on=FEATURE_KEY + ["molecule"], suffixes=("_true", "_posterior")
    )
    if len(merged) != len(truth):
        raise AssignmentFormatError("The posterior does not hold every true row.")

    grouped = merged.groupby(FEATURE_KEY, sort=False)
    is_true = merged["fraction_true"] == 1
    mode_is_true = grouped["fraction_posterior"].transform("idxmax") == merged.index
    mode_hit = (mode_is_true & is_true).groupby([merged[k] for k in FEATURE_KEY]).any()

    mass_on_truth = merged["fraction_posterior"].where(is_true, 0.0)
    mean_true = mass_on_truth.groupby([merged[k] for k in FEATURE_KEY]).sum()

    unknown = merged[merged["molecule"] == UNKNOWN].set_index(FEATURE_KEY)
    unknown_truth = unknown["fraction_true"] == 1
    unknown_mass = unknown["fraction_posterior"]

    def mean(values: pd.Series) -> float:
        return float(values.mean()) if len(values) else float("nan")

    return AssignmentScore(
        n=len(mean_true),
        accuracy=float(mode_hit.mean()),
        mean_true=float(mean_true.mean()),
        n_unknown_features=int(unknown_truth.sum()),
        mean_unknown_of_unknown_features=mean(unknown_mass[unknown_truth]),
        n_molecule_features=int((~unknown_truth).sum()),
        mean_unknown_of_molecule_features=mean(unknown_mass[~unknown_truth]),
    )
