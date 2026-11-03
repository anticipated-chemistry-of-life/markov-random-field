"""The gate on a rung's mean coverage over many simulated scenarios."""

from __future__ import annotations

import validate_replicates as replicates
from src.independent import recovery


def _summary(coverage: float) -> dict:
    n = 128
    block = {"n": n, "coverage": coverage, "passes": recovery.coverage_passes(coverage, n)}
    return {tree: {"alpha": block, "log_nu": block} for tree in replicates.TREES}


def test_scattered_scenarios_pass_when_their_mean_is_nominal():
    # The rung-2 log nu scatter of an exact posterior (#112): most of these
    # scenarios fail the single-scenario band, and the sampler is still right.
    coverages = [0.78, 0.82, 0.83, 0.90, 0.95, 0.98, 0.98, 0.96, 0.90, 0.83, 0.82, 0.82]
    summaries = [_summary(c) for c in coverages]

    assert sum(not s["species"]["log_nu"]["passes"] for s in summaries) > len(summaries) / 2
    assert replicates.apply_gate(replicates.coverage_scores(summaries), gate=3.0)


def test_a_shortfall_every_scenario_shares_fails_although_each_is_in_band():
    coverages = [0.85, 0.86, 0.85, 0.87, 0.86, 0.85, 0.86, 0.87]
    summaries = [_summary(c) for c in coverages]

    assert all(s["species"]["log_nu"]["passes"] for s in summaries)
    assert not replicates.apply_gate(replicates.coverage_scores(summaries), gate=3.0)


def test_a_mean_off_nominal_with_no_scatter_fails():
    summaries = [_summary(0.85)] * 4

    assert not replicates.apply_gate(replicates.coverage_scores(summaries), gate=3.0)
