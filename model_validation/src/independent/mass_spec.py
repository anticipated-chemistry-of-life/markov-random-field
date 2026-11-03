"""The LOTUS+MS scenarios' mass-spec half: runs, features, candidates and scores.

The generative story follows the C++ model's detection table and ADR-0018, and
shares no code with the C++:

* a present molecule is detected in a run with its filter probability, an absent
  one with the contamination probability kappa;
* each detected molecule produces exactly one feature;
* each run also gets Poisson(lambda) *unknown* features, lambda being the expected
  number of molecule features per run under the true field, filter probabilities
  and kappa. So about half of all features are unknown, as on real data;
* every feature proposes n_f ~ Uniform{1..20} structures, the same distribution
  for both kinds, so the list length carries no evidence and the unknown weight
  rho * n_f is exact. A molecule feature lists its true molecule plus other
  molecule leaves, an unknown feature lists molecule leaves only; both draw
  uniformly without replacement. Every proposed structure is a molecule leaf;
* a wrong row scores x ~ Uniform(0,1) and a true row x ~ Beta(5,1), so the log
  correct density is ln 5 + 4 ln x and the log incorrect density is 0.

rho = lambda / (number of molecule leaves): a consequence of the story, not a
choice. The inference run is given rho and nothing else of this module's truth.
"""

from __future__ import annotations

import pathlib
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .indexing import TreeIndex

DATA_FILE = "simulated_ms_data.tsv"
ASSIGNMENT_FILE = "simulated_ms_assignment.tsv"

UNKNOWN = "<unknown>"
"""The reserved name the C++ writes the unknown molecule under."""

MAX_CANDIDATES = 20
CONTAMINATION = 0.01
N_FILTER_LEVELS = 255
"""Filter probabilities live on the byte grid `i / 255`, `i` in `0..254`: the C++
caps the index at 254 (`MAX_FILTER_PROBABILITY_INDEX`), so exactly 1 is not a level."""

CORRECT_SHAPE = 5.0
"""A true row scores Beta(CORRECT_SHAPE, 1): its density is `5 x^4`."""

FILTER_PARAMETER = "filter_proba"
CONTAMINATION_PARAMETER = "contamination_proba"
LAMBDA_PARAMETER = "ms_lambda"
UNKNOWN_RATE_PARAMETER = "ms_unknown_rate"


@dataclass(frozen=True)
class Layout:
    """How many runs there are, and for which species."""

    name: str
    runs_per_species: int
    n_species_with_runs: int | None
    """`None`: every species leaf has runs."""


LAYOUTS = {
    "dense": Layout("dense", runs_per_species=2, n_species_with_runs=None),
    "sparse": Layout("sparse", runs_per_species=1, n_species_with_runs=8),
}


@dataclass(frozen=True)
class Feature:
    species: int
    """Species leaf number."""
    extract: str
    feature: str
    candidates: np.ndarray
    """Molecule leaf numbers of the proposed structures."""
    log_correct: np.ndarray
    true_molecule: int | None
    """`None`: an unknown feature, which has no true row."""


@dataclass
class MassSpecTruth:
    layout: Layout
    filter_indices: np.ndarray
    contamination: float
    expected_features_per_run: float
    unknown_rate: float
    features: list[Feature]
    species_with_runs: np.ndarray

    @property
    def n_runs(self) -> int:
        return len(self.species_with_runs) * self.layout.runs_per_species

    @property
    def filter_probabilities(self) -> np.ndarray:
        return self.filter_indices / 255.0


def simulate(
    rng: np.random.Generator,
    field: np.ndarray,
    layout: Layout,
    contamination: float = CONTAMINATION,
) -> MassSpecTruth:
    """Draw the filter probabilities, then every run's features, under `field`.

    `field` is `[species_leaf, molecule_leaf]`.
    """
    n_species, n_molecules = field.shape
    if n_molecules < MAX_CANDIDATES:
        raise ValueError(
            f"A feature can propose {MAX_CANDIDATES} structures, so the molecule "
            f"tree needs that many leaves; it has {n_molecules}."
        )
    if layout.n_species_with_runs is None:
        species_with_runs = np.arange(n_species)
    else:
        if layout.n_species_with_runs > n_species:
            raise ValueError(
                f"The {layout.name} layout has runs for {layout.n_species_with_runs} "
                f"species, but the tree has {n_species} species leaves."
            )
        species_with_runs = np.sort(
            rng.choice(n_species, layout.n_species_with_runs, replace=False)
        )

    filter_indices = rng.integers(0, N_FILTER_LEVELS, size=n_molecules)
    p_detect = np.where(
        field[species_with_runs], filter_indices / 255.0, contamination
    )
    # Every species has the same number of runs, so the mean over species is the
    # mean over runs.
    expected = float(p_detect.sum(axis=1).mean())
    if expected <= 0.0:
        raise ValueError("No molecule can be detected, so there is nothing to simulate.")

    features: list[Feature] = []
    for row, species in enumerate(species_with_runs):
        for run in range(layout.runs_per_species):
            detected = np.flatnonzero(rng.random(n_molecules) < p_detect[row])
            sources = np.concatenate(
                [detected, np.full(rng.poisson(expected), -1)]
            ).astype(int)
            for number, source in enumerate(rng.permutation(sources)):
                features.append(
                    _draw_feature(
                        rng, int(species), f"extract_{run + 1}", f"feature_{number}",
                        int(source), n_molecules,
                    )
                )
    return MassSpecTruth(
        layout=layout,
        filter_indices=filter_indices,
        contamination=contamination,
        expected_features_per_run=expected,
        unknown_rate=expected / n_molecules,
        features=features,
        species_with_runs=species_with_runs,
    )


def _draw_feature(
    rng: np.random.Generator,
    species: int,
    extract: str,
    feature: str,
    source: int,
    n_molecules: int,
) -> Feature:
    """One feature of molecule `source`, or an unknown one when `source` is -1."""
    n_candidates = int(rng.integers(1, MAX_CANDIDATES + 1))
    if source < 0:
        candidates = rng.choice(n_molecules, n_candidates, replace=False)
    else:
        others = np.delete(np.arange(n_molecules), source)
        candidates = np.concatenate(
            [[source], rng.choice(others, n_candidates - 1, replace=False)]
        )
    candidates = rng.permutation(candidates)

    # (0, 1], so no score is exactly 0 and no log density is infinite.
    u = 1.0 - rng.random(n_candidates)
    # Beta(5, 1) by its inverse CDF; a wrong row is the uniform itself.
    x = np.where(candidates == source, u ** (1.0 / CORRECT_SHAPE), u)
    return Feature(
        species=species,
        extract=extract,
        feature=feature,
        candidates=candidates,
        log_correct=np.log(CORRECT_SHAPE) + (CORRECT_SHAPE - 1.0) * np.log(x),
        true_molecule=None if source < 0 else source,
    )


def write_data(
    path: pathlib.Path,
    truth: MassSpecTruth,
    species: TreeIndex,
    molecules: TreeIndex,
) -> None:
    """The loader's six columns. The incorrect density is flat, so its log is 0."""
    species_names = species.leaf_names()
    molecule_names = np.array(molecules.leaf_names())
    frames = [
        pd.DataFrame(
            {
                "species": species_names[f.species],
                "extract": f.extract,
                "feature": f.feature,
                "molecule": molecule_names[f.candidates],
                "log_correct_density": f.log_correct,
                "log_incorrect_density": 0.0,
            }
        )
        for f in truth.features
    ]
    pd.concat(frames, ignore_index=True).to_csv(path, sep="\t", index=False)


def write_assignment(
    path: pathlib.Path,
    truth: MassSpecTruth,
    species: TreeIndex,
    molecules: TreeIndex,
) -> None:
    """The true assignment, in the layout of the C++ assignment posterior.

    One row per species, extract, feature and explanation -- every candidate and
    the unknown molecule -- with the fraction 1 on the true explanation and 0
    elsewhere, so the two files join row by row.
    """
    species_names = species.leaf_names()
    molecule_names = molecules.leaf_names()
    rows = []
    for f in truth.features:
        for molecule in f.candidates:
            rows.append(
                (
                    species_names[f.species],
                    f.extract,
                    f.feature,
                    molecule_names[molecule],
                    int(molecule == f.true_molecule),
                )
            )
        rows.append(
            (
                species_names[f.species],
                f.extract,
                f.feature,
                UNKNOWN,
                int(f.true_molecule is None),
            )
        )
    pd.DataFrame(
        rows, columns=["species", "extract", "feature", "molecule", "fraction"]
    ).to_csv(path, sep="\t", index=False)


def truth_values(truth: MassSpecTruth, molecules: TreeIndex) -> dict[str, object]:
    """The true parameters under the names the C++ trace gives them.

    A filter probability is the byte index the C++ stores; the probability is the
    index over 255.
    """
    values: dict[str, object] = {
        f"{FILTER_PARAMETER}_{name}": int(index)
        for name, index in zip(molecules.leaf_names(), truth.filter_indices)
    }
    values[CONTAMINATION_PARAMETER] = truth.contamination
    values[LAMBDA_PARAMETER] = truth.expected_features_per_run
    values[UNKNOWN_RATE_PARAMETER] = truth.unknown_rate
    return values
