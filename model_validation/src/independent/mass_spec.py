"""The LOTUS+MS scenarios' mass-spec half: runs, features, candidates and scores.

The generative story follows the C++ model's detection table and ADR-0018, and
shares no code with the C++:

* a present molecule is detected in a run with its filter probability, an absent
  one with the contamination probability kappa;
* each detected molecule produces exactly one feature;
* each run also gets Poisson(lambda) *unknown* features, lambda being the expected
  number of molecule features per run under the true field, filter probabilities
  and kappa. So about half of all features are unknown, as on real data, unless there is an
  outside pool, which adds outside features;
* every feature proposes n_f ~ Uniform{1..20} structures by default, the same distribution
  for both kinds, so the list length carries no evidence and the unknown weight
  rho * n_f is exact. A molecule feature lists its true molecule plus decoys, an unknown
  feature decoys only; both draw uniformly without replacement. Decoys are molecule
  leaves, and with an outside pool (`MassSpecOptions.outside_pool`, ADR-0020) the pool's
  outside molecules too. Each outside molecule is detected in a run with probability p_o
  and each detection makes one feature that lists it;
* by default a wrong row scores x ~ Uniform(0,1) and a true row x ~ Beta(5,1), so the log
  correct density is ln 5 + 4 ln x and the log incorrect density is 0. All of these draw
  settings are the fields of `MassSpecOptions`.

rho = lambda / N, N being the number of molecule leaves plus the pool, and tau = p_o / (1 - p_o):
a consequence of the story, not a choice (the README derives them). The inference run is
given rho and tau and nothing else of this module's truth.
"""

from __future__ import annotations

import pathlib
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .indexing import TreeIndex

DATA_FILE = "simulated_ms_data.tsv"
ASSIGNMENT_FILE = "simulated_ms_assignment.tsv"

UNLISTED = "<unlisted>"
"""The reserved name the C++ writes the rest of the unknown molecule's fraction under, after
each outside structure has its share. The true assignment puts a feature with no
molecule there."""

OUTSIDE_PREFIX = "outside_"
"""The name of an outside molecule is this prefix and its number in the pool."""

CONTAMINATION = 0.01
OUTSIDE_FEATURE_RATIO = 2.0
"""Outside features per molecule feature, near the real export's ratio."""
N_FILTER_LEVELS = 255
"""Filter probabilities live on the byte grid `i / 255`, `i` in `0..254`: the C++
caps the index at 254 (`MAX_FILTER_PROBABILITY_INDEX`), so exactly 1 is not a level."""


@dataclass(frozen=True)
class MassSpecOptions:
    """The draw settings of the mass-spec half, with the values the scenarios use by default."""

    max_candidates: int = 20
    """A feature proposes `n_f ~ Uniform{1..max_candidates}` structures."""
    correct_shape: float = 5.0
    """A true row scores Beta(correct_shape, 1): its density is `correct_shape x^(correct_shape - 1)`,
    and `log_correct_density` is that density's log at every row's score."""
    incorrect_shape: float = 1.0
    """A wrong row scores Beta(incorrect_shape, 1); 1 is Uniform(0, 1)."""
    log_incorrect_density: float = 0.0
    """The `log_incorrect_density` column as written. 0 is a flat density."""
    outside_pool: int = 0
    """The number of outside molecules (ADR-0020): molecules with no leaf in the tree, that
    a feature can list. 0 is no outside structure at all, and leaves every draw as it was."""


FILTER_PARAMETER = "filter_proba"
CONTAMINATION_PARAMETER = "contamination_proba"
LAMBDA_PARAMETER = "ms_lambda"
UNKNOWN_RATE_PARAMETER = "ms_unknown_rate"
OUTSIDE_POOL_PARAMETER = "ms_outside_pool"
OUTSIDE_DETECTION_PARAMETER = "ms_outside_detection"
OUTSIDE_RATE_PARAMETER = "ms_outside_rate"


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
    """Structure numbers of the proposed structures: molecule leaves, then outside molecules."""
    log_correct: np.ndarray
    true_molecule: int | None
    """`None`: an unknown feature, which has no true row."""


@dataclass
class MassSpecTruth:
    layout: Layout
    options: MassSpecOptions
    filter_indices: np.ndarray
    contamination: float
    expected_features_per_run: float
    unknown_rate: float
    outside_detection: float
    """p_o: the probability that an outside molecule is detected in a run. 0 without a pool."""
    features: list[Feature]
    species_with_runs: np.ndarray

    @property
    def outside_rate(self) -> float:
        """tau, the outside rate under which the model is exact: the odds p_o / (1 - p_o)."""
        return self.outside_detection / (1.0 - self.outside_detection)

    @property
    def n_outside_features(self) -> int:
        """The features whose true structure is an outside molecule."""
        n_leaves = len(self.filter_indices)
        return sum(
            f.true_molecule is not None and f.true_molecule >= n_leaves
            for f in self.features
        )

    @property
    def n_runs(self) -> int:
        return len(self.species_with_runs) * self.layout.runs_per_species


def simulate(
    rng: np.random.Generator,
    field: np.ndarray,
    layout: Layout,
    options: MassSpecOptions = MassSpecOptions(),
    contamination: float = CONTAMINATION,
) -> MassSpecTruth:
    """Draw the filter probabilities, then every run's features, under `field`.

    `field` is `[species_leaf, molecule_leaf]`.
    """
    n_species, n_molecules = field.shape
    n_structures = n_molecules + options.outside_pool
    if n_structures < options.max_candidates:
        raise ValueError(
            f"A feature can propose {options.max_candidates} structures, so the molecule "
            f"leaves and the outside pool need that many together; they have {n_structures}."
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

    # Every outside molecule is detected in a run with probability p_o, and p_o is set
    # so that a run has OUTSIDE_FEATURE_RATIO times as many outside features as the
    # molecule features it expects.
    outside_detection = 0.0
    if options.outside_pool > 0:
        outside_detection = OUTSIDE_FEATURE_RATIO * expected / options.outside_pool
        if outside_detection >= 1.0:
            raise ValueError(
                f"The outside pool of {options.outside_pool} is too small for "
                f"{OUTSIDE_FEATURE_RATIO:g} times the {expected:.1f} molecule features "
                f"a run expects: every outside molecule would be detected."
            )

    features: list[Feature] = []
    for row, species in enumerate(species_with_runs):
        for run in range(layout.runs_per_species):
            detected = np.flatnonzero(rng.random(n_molecules) < p_detect[row])
            parts = [detected]
            if options.outside_pool > 0:
                parts.append(
                    n_molecules
                    + np.flatnonzero(rng.random(options.outside_pool) < outside_detection)
                )
            parts.append(np.full(rng.poisson(expected), -1))
            sources = np.concatenate(parts).astype(int)
            for number, source in enumerate(rng.permutation(sources)):
                features.append(
                    _draw_feature(
                        rng, options, int(species), f"extract_{run + 1}",
                        f"feature_{number}", int(source), n_structures,
                    )
                )
    return MassSpecTruth(
        layout=layout,
        options=options,
        filter_indices=filter_indices,
        contamination=contamination,
        expected_features_per_run=expected,
        unknown_rate=expected / n_structures,
        outside_detection=outside_detection,
        features=features,
        species_with_runs=species_with_runs,
    )


def _draw_feature(
    rng: np.random.Generator,
    options: MassSpecOptions,
    species: int,
    extract: str,
    feature: str,
    source: int,
    n_structures: int,
) -> Feature:
    """One feature of structure `source`, or an unknown one when `source` is -1.

    Structures are numbered molecule leaves first, then the outside pool. A feature draws
    its decoys uniformly from all of them.
    """
    n_candidates = int(rng.integers(1, options.max_candidates + 1))
    if source < 0:
        candidates = rng.choice(n_structures, n_candidates, replace=False)
    else:
        others = np.delete(np.arange(n_structures), source)
        candidates = np.concatenate(
            [[source], rng.choice(others, n_candidates - 1, replace=False)]
        )
    candidates = rng.permutation(candidates)

    # (0, 1], so no score is exactly 0 and no log density is infinite.
    u = 1.0 - rng.random(n_candidates)
    # Each score is its distribution's inverse CDF at u, Beta(shape, 1): u ** (1 / shape).
    x = np.where(
        candidates == source,
        u ** (1.0 / options.correct_shape),
        u ** (1.0 / options.incorrect_shape),
    )
    return Feature(
        species=species,
        extract=extract,
        feature=feature,
        candidates=candidates,
        log_correct=np.log(options.correct_shape)
        + (options.correct_shape - 1.0) * np.log(x),
        true_molecule=None if source < 0 else source,
    )


def structure_names(molecules: TreeIndex, outside_pool: int) -> np.ndarray:
    """The name of every structure number: the molecule leaves, then the outside pool.

    An outside molecule is named `outside_<i>`, which no tree node is.
    """
    return np.array(
        molecules.leaf_names() + [f"{OUTSIDE_PREFIX}{i}" for i in range(outside_pool)]
    )


def write_data(
    path: pathlib.Path,
    truth: MassSpecTruth,
    species: TreeIndex,
    molecules: TreeIndex,
) -> None:
    """The loader's six columns, the incorrect density as the options set it."""
    species_names = species.leaf_names()
    molecule_names = structure_names(molecules, truth.options.outside_pool)
    frames = [
        pd.DataFrame(
            {
                "species": species_names[f.species],
                "extract": f.extract,
                "feature": f.feature,
                "molecule": molecule_names[f.candidates],
                "log_correct_density": f.log_correct,
                "log_incorrect_density": truth.options.log_incorrect_density,
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
    `<unlisted>`, which carries the features with no molecule -- with the fraction 1
    on the true explanation and 0 elsewhere, so the files join row by row.
    """
    species_names = species.leaf_names()
    molecule_names = structure_names(molecules, truth.options.outside_pool)
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
                UNLISTED,
                int(f.true_molecule is None),
            )
        )
    pd.DataFrame(
        rows, columns=["species", "extract", "feature", "molecule", "fraction"]
    ).to_csv(path, sep="\t", index=False)


def truth_values(truth: MassSpecTruth, molecules: TreeIndex) -> dict[str, object]:
    """The true parameters, for `simulated_parameters.txt`.

    `filter_proba_<leaf>` and `contamination_proba` are the names the C++ trace
    gives its parameters. `ms_unknown_rate` is the name of the C++ option rho is
    passed by, and `ms_lambda` has no C++ counterpart. A filter probability is the
    byte index the C++ stores; the probability is the index over 255.
    """
    values: dict[str, object] = {
        f"{FILTER_PARAMETER}_{name}": int(index)
        for name, index in zip(molecules.leaf_names(), truth.filter_indices)
    }
    values[CONTAMINATION_PARAMETER] = truth.contamination
    values[LAMBDA_PARAMETER] = truth.expected_features_per_run
    values[UNKNOWN_RATE_PARAMETER] = truth.unknown_rate
    # Only a scenario with a pool names them, so that one without is as it was.
    if truth.options.outside_pool > 0:
        values[OUTSIDE_POOL_PARAMETER] = truth.options.outside_pool
        values[OUTSIDE_DETECTION_PARAMETER] = truth.outside_detection
        values[OUTSIDE_RATE_PARAMETER] = truth.outside_rate
    return values
