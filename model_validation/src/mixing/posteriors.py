"""Where each tree field's posterior landed, and where the two sit on the ridge.

The field's own posterior cannot answer the question this measurement asks: `omega` and the two
alphas trade at constant `a~_s * a~_m`, so a chain riding that ridge moves the two tree fields in
opposite directions and leaves the field where it was (ADR-0005, derivation 3). So both tree field
posteriors are read, and read as a point on the ridge. `mixing_cost/findings.md` reports what they
said.
"""

from __future__ import annotations

import math
import pathlib
from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class TreeFieldPosterior:
    """How often each leaf pair held a one, over the whole leaf-pair space."""

    fractions: np.ndarray

    @property
    def n_cells(self) -> int:
        return int(self.fractions.size)

    @property
    def mean(self) -> float:
        """The tree field's posterior density: the mean over every leaf pair."""
        return float(self.fractions.mean())


def read_tree_field_posterior(path: str | pathlib.Path, n_cells: int) -> TreeFieldPosterior:
    """Read `*_tree_field_posterior.txt` into the whole leaf-pair space.

    The file holds one row per leaf pair that is a one now or was counted a one
    at least once; a cell that is neither is left out (`src/tree/io/write_tree_field.h`).
    An absent row is therefore a posterior of zero and not a missing measurement,
    which is why the rows are scattered by `position` into a space of `n_cells`
    zeros rather than read in file order.
    """
    frame = pd.read_csv(path, sep="\t", usecols=["position", "fraction_of_one"])
    positions = frame["position"].to_numpy(dtype=np.int64)

    if positions.size and (positions.min() < 0 or positions.max() >= n_cells):
        raise ValueError(
            f"{path} names a cell outside the leaf-pair space of {n_cells} cells."
        )
    if np.unique(positions).size != positions.size:
        raise ValueError(f"{path} names the same cell twice.")

    fractions = np.zeros(n_cells, dtype=float)
    fractions[positions] = frame["fraction_of_one"].to_numpy(dtype=float)
    return TreeFieldPosterior(fractions=fractions)


@dataclass(frozen=True)
class PosteriorAgreement:
    """Two tree field posteriors, cell by cell."""

    n_cells: int
    mean_left: float
    mean_right: float
    mean_abs_diff: float
    max_abs_diff: float
    rms_diff: float
    correlation: float


def compare_posteriors(left: TreeFieldPosterior, right: TreeFieldPosterior) -> PosteriorAgreement:
    if left.n_cells != right.n_cells:
        raise ValueError(
            f"The two posteriors cover a different leaf-pair space: {left.n_cells} "
            f"cells against {right.n_cells}."
        )

    difference = left.fractions - right.fractions
    # Two posteriors that are both flat have no correlation to report; nan says
    # so rather than a zero that would read as disagreement.
    if left.fractions.std() == 0.0 or right.fractions.std() == 0.0:
        correlation = float("nan")
    else:
        correlation = float(np.corrcoef(left.fractions, right.fractions)[0, 1])

    return PosteriorAgreement(
        n_cells=left.n_cells,
        mean_left=left.mean,
        mean_right=right.mean,
        mean_abs_diff=float(np.abs(difference).mean()),
        max_abs_diff=float(np.abs(difference).max()),
        rms_diff=float(np.sqrt(np.mean(difference**2))),
        correlation=correlation,
    )


@dataclass(frozen=True)
class RidgePoint:
    """One run's position on the ADR-0005 ridge: the two tree field densities."""

    species: float
    molecules: float

    @property
    def product(self) -> float:
        """What the field sees. Constant along the ridge."""
        return self.species * self.molecules

    @property
    def log_product(self) -> float:
        """The coordinate across the ridge: it moves only when the field's density does."""
        return math.log(self.species) + math.log(self.molecules)

    @property
    def log_ratio(self) -> float:
        """The coordinate along the ridge: it moves when the two trees trade."""
        return math.log(self.species) - math.log(self.molecules)


def ridge_point(species: TreeFieldPosterior, molecules: TreeFieldPosterior) -> RidgePoint:
    return RidgePoint(species=species.mean, molecules=molecules.mean)


@dataclass(frozen=True)
class RidgeShift:
    """How far two runs sit apart, split into the two directions that mean different things."""

    along: float
    across: float


def ridge_shift(left: RidgePoint, right: RidgePoint) -> RidgeShift:
    """Signed, so that the direction reads: a positive `along` puts `left` further
    towards the species tree than `right`.

    `across` is the one that says the two runs disagree about the field. `along`
    on its own is the ridge, and two runs that differ there and nowhere else are
    the signature the block update existed to prevent.
    """
    return RidgeShift(
        along=left.log_ratio - right.log_ratio,
        across=left.log_product - right.log_product,
    )
