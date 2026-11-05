"""The validator's per-clique coverage, the gates built on it, and its field false-presence score."""

from __future__ import annotations

import re

import ms_fixture
import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

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


def _summary(
    coverage: float,
    spearman: float = 0.9,
    *,
    log_nu_coverage: float | None = None,
    node_states_pinned: bool = True,
) -> dict:
    n = 128

    def block(value: float) -> dict:
        return {"n": n, "coverage": value, "passes": recovery.coverage_passes(value, n)}

    log_nu = coverage if log_nu_coverage is None else log_nu_coverage
    summary = {
        tree: {
            "alpha": block(coverage),
            "log_nu": block(log_nu),
            "branch_lengths": {"n": 254, "spearman": spearman},
        }
        for tree in validate.TREES
    }
    summary["node_states_pinned"] = node_states_pinned
    return summary


def test_a_rung_with_calibrated_cliques_passes_whatever_its_branch_length_rank():
    # Once node states are free the branch lengths are barely identified, so a
    # low Spearman is the right answer there, not a failure (#110).
    assert validate._apply_gates(_summary(coverage=0.9, spearman=0.1))


def test_clique_coverage_outside_the_band_fails_in_either_direction():
    assert not validate._apply_gates(_summary(coverage=0.80))
    assert not validate._apply_gates(_summary(coverage=0.99))


def test_log_nu_coverage_is_gated_only_where_the_node_states_are_pinned():
    # With the node states inferred, one shared var log nu sets every clique's
    # interval at once. An exact posterior then falls outside the binomial band
    # on about half of all datasets (#112).
    assert not validate._apply_gates(_summary(coverage=0.9, log_nu_coverage=0.75))
    assert validate._apply_gates(
        _summary(coverage=0.9, log_nu_coverage=0.75, node_states_pinned=False)
    )


def test_alpha_coverage_is_gated_whether_or_not_the_node_states_are_pinned():
    assert not validate._apply_gates(
        _summary(coverage=0.75, log_nu_coverage=0.9, node_states_pinned=False)
    )


def test_only_the_rung_that_pins_the_node_states_pins_them():
    assert validate._pins_node_states("rung1_pin_field_and_states")
    for rung in ("rung2_pin_field", "rung3_from_data_s", "rung3_from_data_lm"):
        assert not validate._pins_node_states(rung)


def test_field_false_presence_counts_only_absent_cells_with_no_lotus_record_in_species_with_runs(
    tmp_path,
):
    # The fixture lists five such cells: (s1,m2) 0.1, (s1,m4) 0.2, (s2,m1) 0.3,
    # (s2,m4) 0.4, and (s2,m3), which the posterior omits and so is 0. They mean 0.2.
    # Three kinds of cell sit at high posterior and must stay out of the score: a
    # present cell, an absent cell with a false LOTUS record, and an absent cell in a
    # species with no run.
    base = ms_fixture.write_scenario(tmp_path / "scenario")
    ms_fixture.write_rung(base, "rung3_from_data_lm")

    summary = validate.score(base, "rung3_from_data_lm")

    false_presence = summary["mass_spec"]["field_false_presence"]
    assert false_presence["n"] == 5
    assert false_presence["mean"] == pytest.approx(0.2)


def test_the_validator_prints_the_field_false_presence_of_a_mass_spec_rung(tmp_path):
    base = ms_fixture.write_scenario(tmp_path / "scenario")
    ms_fixture.write_rung(base, "rung3_from_data_lm")

    result = CliRunner().invoke(validate.main, [str(base), "rung3_from_data_lm"])

    assert result.exit_code == 0, result.output
    printed = re.search(r"field false presence\s+(\d+\.\d+)\s+\(n=(\d+)", result.output)
    assert printed, result.output
    assert float(printed.group(1)) == pytest.approx(0.2, abs=1e-3)
    assert int(printed.group(2)) == 5
