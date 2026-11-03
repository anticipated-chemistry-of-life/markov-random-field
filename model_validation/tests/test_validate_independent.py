"""The validator's per-clique coverage and the gates built on it."""

from __future__ import annotations

import numpy as np
import pandas as pd

import validate_independent as validate
from src.independent import recovery


def _draws_0_to_99(names: list[str]) -> pd.DataFrame:
    return pd.DataFrame({name: np.arange(100.0) for name in names})


def test_clique_coverage_uses_the_interval_the_band_is_centred_on():
    # Draws 0..99: the central 90% interval is about [4.95, 94.05] and the 95% one
    # about [2.48, 96.53]. A truth of 3 is inside the second and outside the first,
    # so scoring it as covered would test a 95% interval against a 90% band (#110).
    trace = _draws_0_to_99(["species_alpha_a", "species_alpha_b"])
    truth = {"species_alpha_a": 50.0, "species_alpha_b": 3.0}

    assert validate._coverage(trace, truth, "species_alpha_")["coverage"] == 0.5


def _summary(coverage: float, spearman: float) -> dict:
    n = 128
    block = {"n": n, "coverage": coverage, "passes": recovery.coverage_passes(coverage, n)}
    return {
        tree: {
            "alpha": block,
            "log_nu": block,
            "branch_lengths": {"n": 254, "spearman": spearman},
        }
        for tree in validate.TREES
    }


def test_a_rung_with_calibrated_cliques_passes_whatever_its_branch_length_rank():
    # Once node states are free the branch lengths are barely identified, so a
    # low Spearman is the right answer there, not a failure (#110).
    assert validate._apply_gates(_summary(coverage=0.9, spearman=0.1))


def test_clique_coverage_outside_the_band_fails_in_either_direction():
    assert not validate._apply_gates(_summary(coverage=0.80, spearman=0.9))
    assert not validate._apply_gates(_summary(coverage=0.99, spearman=0.9))
