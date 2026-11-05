"""The gate on a rung's mean coverage over many simulated scenarios, and the difference report."""

from __future__ import annotations

import re

import ms_fixture
import pytest
from click.testing import CliRunner

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


def _table(output: str) -> dict[str, tuple[float, float, float]]:
    """Each row of the report: its mean, the baseline mean and the difference."""
    rows = {}
    for line in output.splitlines():
        row = re.match(r"\s+(.+?)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+([+-]\d+\.\d+)\s+\d+$", line)
        if row:
            label, mean, baseline, difference = row.groups()
            rows[label] = (float(mean), float(baseline), float(difference))
    return rows


def _row(output: str, label: str) -> tuple[float, float, float]:
    rows = _table(output)
    assert label in rows, output
    return rows[label]


def test_the_difference_report_prints_the_mean_over_seeds_the_baseline_mean_and_the_difference(
    tmp_path,
):
    # The field false presence of the fixture is 0.2 times the rung's field scale.
    # Over seeds 1-3 the scenario's scales average 1.0 and the baseline's 0.5, so the
    # means are 0.2 and 0.1. One seed or a different average would give other numbers.
    scenario_scales = {1: 0.5, 2: 1.0, 3: 1.5}
    baseline_scales = {1: 0.25, 2: 0.5, 3: 0.75}
    bases = []
    for seed in (1, 2, 3):
        base = ms_fixture.write_scenario(tmp_path / f"seed{seed}", seed=seed)
        ms_fixture.write_rung(base, "rung3_from_data_lm", field_scale=scenario_scales[seed])
        ms_fixture.write_rung(base, "rung2_pin_field", field_scale=baseline_scales[seed])
        bases.append(str(base))

    result = CliRunner().invoke(
        replicates.main,
        ["rung3_from_data_lm", *bases, "--baseline-rung", "rung2_pin_field"],
    )

    assert result.exit_code == 0, result.output
    assert _row(result.output, "field false presence") == pytest.approx((0.2, 0.1, 0.1))


def test_the_baseline_family_is_found_by_pattern_and_only_for_the_scenario_seeds(tmp_path):
    # False presence is 0.2 times the field scale. The scenario family averages 0.2
    # over seeds 1-3 and the baseline family 0.1. Seed 42 also matches the pattern, at
    # a false presence of 0: it is no scenario seed, so it must stay out of the baseline.
    scenario_scales = {1: 0.5, 2: 1.0, 3: 1.5}
    baseline_scales = {1: 0.25, 2: 0.5, 3: 0.75, 42: 0.0}
    scenarios = []
    for seed, scale in scenario_scales.items():
        base = ms_fixture.write_scenario(tmp_path / f"seed{seed}_ms_dense_distorted", seed=seed)
        ms_fixture.write_rung(base, "rung3_from_data_lm", field_scale=scale)
        scenarios.append(str(base))
    for seed, scale in baseline_scales.items():
        base = ms_fixture.write_scenario(tmp_path / f"seed{seed}_ms_dense", seed=seed)
        ms_fixture.write_rung(base, "rung3_from_data_lm", field_scale=scale)

    result = CliRunner().invoke(
        replicates.main,
        ["rung3_from_data_lm", *scenarios, "--baseline", str(tmp_path / "*_ms_dense")],
    )

    assert result.exit_code == 0, result.output
    assert _row(result.output, "field false presence") == pytest.approx((0.2, 0.1, 0.1))


def _scenario_with_rung(tmp_path, name: str, *, seed: int, layout: str = "dense") -> str:
    base = ms_fixture.write_scenario(tmp_path / name, seed=seed, layout=layout)
    ms_fixture.write_rung(base, "rung3_from_data_lm")
    return str(base)


def test_a_scenario_seed_that_the_baseline_family_lacks_is_an_error(tmp_path):
    scenarios = [
        _scenario_with_rung(tmp_path, f"seed{seed}_ms_dense_distorted", seed=seed)
        for seed in (1, 2)
    ]
    _scenario_with_rung(tmp_path, "seed1_ms_dense", seed=1)

    result = CliRunner().invoke(
        replicates.main,
        ["rung3_from_data_lm", *scenarios, "--baseline", str(tmp_path / "*_ms_dense")],
    )

    assert result.exit_code != 0
    assert "seed 2" in result.output


def test_a_baseline_with_another_layout_than_its_scenario_is_an_error(tmp_path):
    # The directory names carry no layout, so only the error can name the two.
    scenario = _scenario_with_rung(tmp_path, "seed1_scenario", seed=1, layout="dense")
    _scenario_with_rung(tmp_path, "seed1_baseline", seed=1, layout="sparse")

    result = CliRunner().invoke(
        replicates.main,
        ["rung3_from_data_lm", scenario, "--baseline", str(tmp_path / "*_baseline")],
    )

    assert result.exit_code != 0
    assert "dense" in result.output and "sparse" in result.output


def test_a_pattern_that_matches_two_baseline_scenarios_of_one_seed_is_an_error(tmp_path):
    # `*_ms_dense*` matches the baseline family and the family under test alike. Of the
    # two scenarios of seed 1, taking either one would be a guess.
    scenario = _scenario_with_rung(tmp_path, "seed1_ms_dense_distorted", seed=1)
    _scenario_with_rung(tmp_path, "seed1_ms_dense", seed=1)

    result = CliRunner().invoke(
        replicates.main,
        ["rung3_from_data_lm", scenario, "--baseline", str(tmp_path / "*_ms_dense*")],
    )

    assert result.exit_code != 0
    assert "seed 1" in result.output


def test_the_difference_report_takes_no_gate(tmp_path):
    scenario = _scenario_with_rung(tmp_path, "seed1", seed=1)

    result = CliRunner().invoke(
        replicates.main,
        ["rung3_from_data_lm", scenario, "--baseline-rung", "rung3_from_data_lm", "--gate", "3"],
    )

    assert result.exit_code != 0
    assert "no gate" in result.output


def test_the_difference_report_has_a_row_for_each_score_it_compares(tmp_path):
    # The fixture's numbers (see ms_fixture) are worked by hand. The assignment mass on
    # the true assignment is (0.9 + 0.7 + 0.7) / 3. No molecule of the fixture is present
    # in 2, 3-9 or 10+ runs, so those strata have no row.
    base = ms_fixture.write_scenario(tmp_path / "scenario")
    ms_fixture.write_rung(base, "rung3_from_data_lm", molecule_mass=0.9, unknown_mass=0.7)

    result = CliRunner().invoke(
        replicates.main,
        ["rung3_from_data_lm", str(base), "--baseline-rung", "rung3_from_data_lm"],
    )

    assert result.exit_code == 0, result.output
    means = {label: mean for label, (mean, _, _) in _table(result.output).items()}
    assert means == pytest.approx(
        {
            "assignment: mean posterior on the true assignment": 0.7667,
            "assignment: P(unknown), unknown features": 0.7,
            "field false presence": 0.2,
            "AUC, species with runs": 1.0,
            "AUC, species without runs": 0.6667,
            "species alpha coverage": 0.75,
            "species log nu coverage": 0.25,
            "molecules alpha coverage": 0.3333,
            "molecules log nu coverage": 1.0,
            "kappa coverage": 0.0,
            "filter probability coverage, runs 0": 0.5,
            "filter probability coverage, runs 1": 1.0,
        },
        abs=5e-5,
    )


def test_a_seed_without_a_value_in_either_family_drops_out_of_both_means(tmp_path):
    # False presence is 0.2 times the field scale. In the baseline's seed 2 every cell of
    # a species with a run has a LOTUS record, so no cell counts and the baseline has no
    # score there. Seed 2 must leave the scenario's mean too. Over seeds 1 and 3 the
    # scenario averages (0.1 + 0.2) / 2 and the baseline (0.05 + 0.15) / 2. Seed 2 would
    # move the scenario's mean to (0.1 + 0.3 + 0.2) / 3 = 0.2.
    recorded = tuple((s, m) for s in ms_fixture.RUN_SPECIES for m in ms_fixture.MOLECULES)
    scenario_scales = {1: 0.5, 2: 1.5, 3: 1.0}
    baseline_scales = {1: 0.25, 2: 1.0, 3: 0.75}
    scenarios = []
    for seed in (1, 2, 3):
        scenario = ms_fixture.write_scenario(tmp_path / f"scenario{seed}", seed=seed)
        ms_fixture.write_rung(scenario, "rung3_from_data_lm", field_scale=scenario_scales[seed])
        scenarios.append(str(scenario))
        baseline = ms_fixture.write_scenario(
            tmp_path / f"baseline{seed}",
            seed=seed,
            lotus=recorded if seed == 2 else ms_fixture.LOTUS,
        )
        ms_fixture.write_rung(baseline, "rung3_from_data_lm", field_scale=baseline_scales[seed])

    result = CliRunner().invoke(
        replicates.main,
        ["rung3_from_data_lm", *scenarios, "--baseline", str(tmp_path / "baseline*")],
    )

    assert result.exit_code == 0, result.output
    assert _row(result.output, "field false presence") == pytest.approx((0.15, 0.1, 0.05))


def test_a_directory_with_no_scenario_in_it_is_left_out_of_the_baseline_family(tmp_path):
    scenario = _scenario_with_rung(tmp_path, "seed1_ms_dense_distorted", seed=1)
    _scenario_with_rung(tmp_path, "seed1_ms_dense", seed=1)
    (tmp_path / "archive_ms_dense").mkdir()

    result = CliRunner().invoke(
        replicates.main,
        ["rung3_from_data_lm", scenario, "--baseline", str(tmp_path / "*_ms_dense")],
    )

    assert result.exit_code == 0, result.output
    assert _row(result.output, "field false presence") == pytest.approx((0.2, 0.2, 0.0))
