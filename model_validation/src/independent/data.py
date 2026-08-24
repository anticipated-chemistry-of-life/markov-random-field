"""The two observation models, derived from the C++ likelihoods.

Both take the latent field and produce a binary observation of the same shape.

**Simple error model.** Every cell is independently misreported with probability
`epsilon` (`simple_error_model::draw_D_given_Y`,
src/simple_error_model/TSimpleErrorModelMath.h:119). Dense: every cell is
observed.

**LOTUS.** A record is reported for a present pair with probability equal to the
research effort spent on it, and for an absent pair with probability
`error_rate` (`TLotus::_calculate_probability_of_L_given_x`, src/TLotus.cpp:202).
Research effort is a product across the kept dimensions of
`1 - exp(-gamma_i * papers_i[leaf])` (src/TLotus.cpp:171). Both trees are kept by
default, so both need paper counts.
"""

from __future__ import annotations

import numpy as np


def sample_paper_counts(rng: np.random.Generator, n_leaves: int) -> np.ndarray:
    """Per-leaf publication counts. Matches `Tree.generate_papers_number`."""
    return rng.poisson(3, size=n_leaves) + 1


def occurrence_count(papers: np.ndarray) -> np.ndarray:
    """The paper counts as the C++ actually uses them: `log(count + 1)`.

    `TTree::get_paper_counts` applies this transform when reading the file
    (src/tree/TTree.h:517), so research effort is driven by the log count, not
    the raw one. Using raw counts here would make the simulated LOTUS data
    inconsistent with the likelihood being fitted, which is indistinguishable
    from an inference bug.
    """
    return np.log(np.asarray(papers, dtype=float) + 1.0)


def research_effort(
    species_papers: np.ndarray,
    molecule_papers: np.ndarray,
    gamma: float,
) -> np.ndarray:
    """`[species_leaf, molecule_leaf]` probability that a present pair is reported."""
    return np.outer(
        1.0 - np.exp(-gamma * occurrence_count(species_papers)),
        1.0 - np.exp(-gamma * occurrence_count(molecule_papers)),
    )


def simulate_lotus(
    rng: np.random.Generator,
    field: np.ndarray,
    species_papers: np.ndarray,
    molecule_papers: np.ndarray,
    gamma: float,
    error_rate: float,
) -> np.ndarray:
    """Draw LOTUS records given the field."""
    probability = np.where(
        field, research_effort(species_papers, molecule_papers, gamma), error_rate
    )
    return rng.random(field.shape) < probability


def simulate_simple_error(
    rng: np.random.Generator, field: np.ndarray, epsilon: float
) -> np.ndarray:
    """Draw the simple error model's observation given the field."""
    return field ^ (rng.random(field.shape) < epsilon)
