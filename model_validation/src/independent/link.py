"""The link between the two tree fields and the field, derived from the maths.

Each tree owns a complete leaf-level view of the field -- its **tree field** --
and the field is a noisy reconciliation of the two: corrupt each tree field cell
independently with probability `omega`, then take the AND. A corrupted cell reads
1 with probability `1 - omega` when it is truly 1 and `omega` when it is truly 0,
and the two corruptions are independent, so the AND is a product of two
independent factors:

    P(Y = 1 | Z_s, Z_m) = c(Z_s) * c(Z_m),   c(1) = 1 - omega,  c(0) = omega

`omega` lives strictly inside `(0, 0.5)`. At 0 the link is the deterministic AND;
at or above 0.5 the tree fields are anti-correlated with the field. See
`docs/adr/0005-each-tree-owns-its-leaf-level-field.md`, which this module is
derived from -- not from `src/field/TFieldMath.h`, which is the implementation
this harness exists to disagree with.

The draw uses the four-cell table above and never the buckets. `buckets`,
`prob_for_bucket` and `link_counters` exist only so that `compare_fields.py` can
hold the C++'s own six traced counters against a tally of the same cells. Putting
the collapse in the draw would let one bucketing convention hide inside both
sides of that comparison.
"""

from __future__ import annotations

import numpy as np

N_BUCKETS = 3
"""The number of tree fields in state 1 at one leaf pair, so 0, 1 or 2."""


def check_error_probability(omega: float) -> float:
    """The open interval `(0, 0.5)`, as a statement about the model."""
    if not 0.0 < omega < 0.5:
        raise ValueError(
            f"An error probability of {omega} is outside the open interval (0, 0.5) that "
            "ADR-0005 puts it on."
        )
    return float(omega)


def corrupted(states: np.ndarray, omega: float) -> np.ndarray:
    """`P(a corrupted read of this cell is 1)`, cell by cell."""
    check_error_probability(omega)
    return np.where(np.asarray(states, dtype=bool), 1.0 - omega, omega)


def prob_field_is_one(
    species_field: np.ndarray, molecule_field: np.ndarray, omega: float
) -> np.ndarray:
    """`P(Y = 1 | Z_s, Z_m, omega)` for every leaf pair."""
    return corrupted(species_field, omega) * corrupted(molecule_field, omega)


def sample_field(
    rng: np.random.Generator,
    species_field: np.ndarray,
    molecule_field: np.ndarray,
    omega: float,
) -> np.ndarray:
    """Draw the field from the two tree fields. One Bernoulli per leaf pair."""
    if species_field.shape != molecule_field.shape:
        raise ValueError(
            f"The two tree fields are {species_field.shape} and {molecule_field.shape}; "
            "they are both indexed [species_leaf, molecule_leaf] and must agree."
        )
    return rng.random(species_field.shape) < prob_field_is_one(
        species_field, molecule_field, omega
    )


def buckets(species_field: np.ndarray, molecule_field: np.ndarray) -> np.ndarray:
    """The number of tree fields in state 1, which is all the link depends on."""
    return np.asarray(species_field, dtype=np.int64) + np.asarray(
        molecule_field, dtype=np.int64
    )


def prob_for_bucket(bucket: np.ndarray | int, omega: float) -> np.ndarray:
    """`P_k = (1 - omega)^k * omega^(2 - k)`, the three distinct link probabilities."""
    check_error_probability(omega)
    k = np.asarray(bucket, dtype=float)
    return (1.0 - omega) ** k * omega ** (2.0 - k)


def link_counters(
    species_field: np.ndarray, molecule_field: np.ndarray, field: np.ndarray
) -> np.ndarray:
    """The link's sufficient statistic, `n(bucket, field state)`. Shape `(3, 2)`."""
    bucket = buckets(species_field, molecule_field).reshape(-1)
    y = np.asarray(field, dtype=np.int64).reshape(-1)
    counters = np.zeros((N_BUCKETS, 2), dtype=np.int64)
    np.add.at(counters, (bucket, y), 1)
    return counters


def adjusted_rate(alpha: np.ndarray | float, omega: float) -> np.ndarray:
    """`a~ = omega + (1 - 2 omega) * alpha`: how often a corrupted tree field reads 1.

    A tree field cell at its stationary rate `alpha` reads 1 after corruption with
    this probability, so the expected density of the field is the *product* of the
    two trees' adjusted rates (ADR-0005, derivation 3). The field's own density
    therefore says nothing about how the rate splits between the two trees.
    """
    check_error_probability(omega)
    return omega + (1.0 - 2.0 * omega) * np.asarray(alpha, dtype=float)
