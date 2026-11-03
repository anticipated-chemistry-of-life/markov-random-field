"""Scores for how well inference recovers the field and the detection parameters.

Pure functions over arrays and frames, so the validator script only reads files.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import binom, rankdata

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
