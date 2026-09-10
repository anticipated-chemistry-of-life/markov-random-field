"""The two instruments the mixing measurement reads.

The measurement itself is a pair of chains that take minutes to run, so what is
checked here is the arithmetic they are read with: an autocorrelation time
against series whose answer is known in closed form, and a tree field posterior
against files whose cells are counted by hand.
"""

from __future__ import annotations

import numpy as np
import pytest

from src.mixing.posteriors import (
    compare_posteriors,
    read_tree_field_posterior,
    ridge_point,
    ridge_shift,
)
from src.mixing.posteriors import RidgePoint
from src.mixing.traces import (
    autocorrelation,
    effective_sample_size,
    half_split_drift,
    integrated_autocorrelation_time,
    read_joint_density_factors,
    summarise_trace,
)


def ridge_point_of(species: float, molecules: float) -> RidgePoint:
    return RidgePoint(species=species, molecules=molecules)


def ar1(phi: float, n: int, seed: int) -> np.ndarray:
    """An AR(1) series, whose autocorrelation and autocorrelation time are known."""
    rng = np.random.default_rng(seed)
    noise = rng.normal(scale=np.sqrt(1.0 - phi**2), size=n)
    x = np.empty(n)
    x[0] = rng.normal()
    for i in range(1, n):
        x[i] = phi * x[i - 1] + noise[i]
    return x


# ---------------------------------------------------------------------------
# autocorrelation
# ---------------------------------------------------------------------------


def test_autocorrelation_starts_at_one():
    rho = autocorrelation(ar1(0.5, 512, seed=1))
    assert rho[0] == pytest.approx(1.0)


def test_autocorrelation_of_an_alternating_series_is_minus_one_at_lag_one():
    # Exactly anti-correlated, up to the edge term a finite series carries.
    x = np.resize([1.0, -1.0], 1024)
    assert autocorrelation(x, max_lag=1)[1] == pytest.approx(-1.0, abs=2e-3)


def test_autocorrelation_of_ar1_decays_geometrically():
    phi = 0.7
    rho = autocorrelation(ar1(phi, 200_000, seed=2), max_lag=5)
    for lag in range(1, 6):
        assert rho[lag] == pytest.approx(phi**lag, abs=0.02)


def test_autocorrelation_stops_at_the_requested_lag():
    assert autocorrelation(ar1(0.5, 100, seed=3), max_lag=7).size == 8


def test_autocorrelation_rejects_a_flat_trace():
    # No variance, so nothing to divide by: a trace that never moved has no
    # autocorrelation rather than a perfect one.
    with pytest.raises(ValueError, match="does not vary"):
        autocorrelation(np.full(64, 3.0))


def test_autocorrelation_rejects_a_trace_of_one_sample():
    with pytest.raises(ValueError, match="at least"):
        autocorrelation(np.array([1.0]))


# ---------------------------------------------------------------------------
# the autocorrelation time, and the sample size it leaves
# ---------------------------------------------------------------------------


def test_autocorrelation_time_of_white_noise_is_about_one():
    rng = np.random.default_rng(4)
    tau = integrated_autocorrelation_time(rng.normal(size=100_000))
    assert tau == pytest.approx(1.0, abs=0.15)


@pytest.mark.parametrize("phi", [0.5, 0.8, 0.9])
def test_autocorrelation_time_of_ar1_matches_the_closed_form(phi):
    # tau = (1 + phi) / (1 - phi) for an AR(1) process.
    expected = (1.0 + phi) / (1.0 - phi)
    tau = integrated_autocorrelation_time(ar1(phi, 500_000, seed=5))
    assert tau == pytest.approx(expected, rel=0.1)


def test_autocorrelation_time_grows_with_the_correlation():
    taus = [integrated_autocorrelation_time(ar1(phi, 100_000, seed=6)) for phi in (0.0, 0.5, 0.9)]
    assert taus[0] < taus[1] < taus[2]


def test_effective_sample_size_is_the_length_over_the_autocorrelation_time():
    x = ar1(0.8, 20_000, seed=7)
    assert effective_sample_size(x) == pytest.approx(x.size / integrated_autocorrelation_time(x))


def test_effective_sample_size_of_white_noise_is_about_the_whole_trace():
    rng = np.random.default_rng(8)
    x = rng.normal(size=50_000)
    assert effective_sample_size(x) == pytest.approx(x.size, rel=0.15)


# ---------------------------------------------------------------------------
# whether the chain had reached what its autocorrelation time describes
# ---------------------------------------------------------------------------


def test_a_stationary_trace_does_not_drift():
    rng = np.random.default_rng(10)
    assert half_split_drift(rng.normal(size=100_000)) == pytest.approx(0.0, abs=0.05)


def test_a_ramp_drifts_by_the_width_of_its_own_spread():
    # A trace rising linearly has its two half-means half the range apart and a
    # standard deviation of range / sqrt(12), so the drift is sqrt(12) / 2 --
    # whatever the range is.
    assert half_split_drift(np.linspace(0.0, 7.0, 100_000)) == pytest.approx(
        np.sqrt(12.0) / 2.0, rel=1e-3
    )


def test_the_drift_is_signed_towards_where_the_trace_went():
    rising = np.linspace(0.0, 1.0, 1000)
    assert half_split_drift(rising) > 0
    assert half_split_drift(rising[::-1]) < 0


def test_the_drift_of_a_flat_trace_has_nothing_to_divide_by():
    with pytest.raises(ValueError, match="does not vary"):
        half_split_drift(np.full(64, 2.0))


def test_summarise_trace_reports_what_the_report_prints():
    x = ar1(0.8, 4_000, seed=9)
    summary = summarise_trace(x)
    assert summary.n_samples == 4_000
    assert summary.mean == pytest.approx(float(np.mean(x)))
    assert summary.sd == pytest.approx(float(np.std(x, ddof=1)))
    assert summary.autocorrelation_time == pytest.approx(integrated_autocorrelation_time(x))
    assert summary.effective_sample_size == pytest.approx(effective_sample_size(x))
    assert summary.acf_lag_one == pytest.approx(autocorrelation(x, max_lag=1)[1])
    assert summary.drift == pytest.approx(half_split_drift(x))


# ---------------------------------------------------------------------------
# reading the joint density trace
# ---------------------------------------------------------------------------


JOINT_DENSITY = (
    "species_node_state\tmolecules_node_state\tlink\tdata\tjoint_density\n"
    "-1.0\t-2.0\t-3.0\t-4.0\t-10.0\n"
    "-1.5\t-2.5\t-3.5\t-4.5\t-12.0\n"
    "-1.25\t-2.25\t-3.25\t-4.25\t-11.0\n"
    "-1.75\t-2.75\t-3.75\t-4.75\t-13.0\n"
)


def test_read_joint_density_takes_every_row_by_default(tmp_path):
    path = tmp_path / "acol_joint_density.txt"
    path.write_text(JOINT_DENSITY)
    factors = read_joint_density_factors(path)
    assert factors["joint_density"].tolist() == [-10.0, -12.0, -11.0, -13.0]


def test_read_joint_density_drops_the_burn_in_rows(tmp_path):
    # The trace is written from the first iteration on, burn-in included: nothing
    # in the writer knows the chain has not started yet.
    path = tmp_path / "acol_joint_density.txt"
    path.write_text(JOINT_DENSITY)
    factors = read_joint_density_factors(path, burn_in_rows=2)
    assert factors["joint_density"].tolist() == [-11.0, -13.0]


def test_read_joint_density_factors_keeps_every_column(tmp_path):
    # The columns are the point of the file: a single total says a chain drifts,
    # the columns say which factor is dragging it.
    path = tmp_path / "acol_joint_density.txt"
    path.write_text(JOINT_DENSITY)
    factors = read_joint_density_factors(path, burn_in_rows=2)
    assert list(factors) == [
        "species_node_state",
        "molecules_node_state",
        "link",
        "data",
        "joint_density",
    ]
    assert factors["link"].tolist() == [-3.25, -3.75]


def test_read_joint_density_refuses_to_drop_the_whole_trace(tmp_path):
    path = tmp_path / "acol_joint_density.txt"
    path.write_text(JOINT_DENSITY)
    with pytest.raises(ValueError, match="burn-in"):
        read_joint_density_factors(path, burn_in_rows=4)


# ---------------------------------------------------------------------------
# tree field posteriors
# ---------------------------------------------------------------------------


def write_posterior(path, rows) -> None:
    """`rows` are (position, Z_state, fraction_of_one); the leaf names do not matter here."""
    lines = ["position\tZ_state\tspecies\tmolecules\tfraction_of_one"]
    for position, state, fraction in rows:
        lines.append(f"{position}\t{state}\tspecies_1\tmolecules_1\t{fraction}")
    path.write_text("\n".join(lines) + "\n")


def test_read_tree_field_posterior_scatters_by_position(tmp_path):
    path = tmp_path / "posterior.txt"
    write_posterior(path, [(3, 1, 0.5), (0, 0, 0.25)])
    posterior = read_tree_field_posterior(path, n_cells=5)
    assert posterior.fractions.tolist() == [0.25, 0.0, 0.0, 0.5, 0.0]


def test_a_cell_the_file_left_out_is_a_zero_and_not_a_gap(tmp_path):
    # The writer drops a cell that is a zero now and was never counted a one, so
    # an absent row is a posterior of zero rather than a missing measurement.
    path = tmp_path / "posterior.txt"
    write_posterior(path, [(1, 1, 1.0)])
    assert read_tree_field_posterior(path, n_cells=4).mean == pytest.approx(0.25)


def test_the_mean_divides_by_the_whole_leaf_pair_space(tmp_path):
    path = tmp_path / "posterior.txt"
    write_posterior(path, [(0, 1, 0.5), (1, 1, 0.5)])
    assert read_tree_field_posterior(path, n_cells=100).mean == pytest.approx(0.01)


def test_read_tree_field_posterior_rejects_a_position_outside_the_space(tmp_path):
    path = tmp_path / "posterior.txt"
    write_posterior(path, [(9, 1, 0.5)])
    with pytest.raises(ValueError, match="outside"):
        read_tree_field_posterior(path, n_cells=5)


def test_read_tree_field_posterior_rejects_a_repeated_position(tmp_path):
    path = tmp_path / "posterior.txt"
    write_posterior(path, [(2, 1, 0.5), (2, 1, 0.25)])
    with pytest.raises(ValueError, match="twice"):
        read_tree_field_posterior(path, n_cells=5)


def test_two_identical_posteriors_agree_exactly(tmp_path):
    path = tmp_path / "posterior.txt"
    write_posterior(path, [(0, 1, 0.5), (2, 1, 0.75)])
    posterior = read_tree_field_posterior(path, n_cells=4)
    agreement = compare_posteriors(posterior, posterior)
    assert agreement.mean_abs_diff == 0.0
    assert agreement.max_abs_diff == 0.0
    assert agreement.rms_diff == 0.0
    assert agreement.correlation == pytest.approx(1.0)


def test_agreement_reports_the_largest_cell_difference(tmp_path):
    left, right = tmp_path / "left.txt", tmp_path / "right.txt"
    write_posterior(left, [(0, 1, 0.5), (1, 1, 0.5)])
    write_posterior(right, [(0, 1, 0.5), (1, 1, 0.9)])
    agreement = compare_posteriors(
        read_tree_field_posterior(left, n_cells=4), read_tree_field_posterior(right, n_cells=4)
    )
    assert agreement.n_cells == 4
    assert agreement.max_abs_diff == pytest.approx(0.4)
    assert agreement.mean_abs_diff == pytest.approx(0.1)
    assert agreement.mean_left == pytest.approx(0.25)
    assert agreement.mean_right == pytest.approx(0.35)


def test_agreement_refuses_two_different_leaf_pair_spaces(tmp_path):
    path = tmp_path / "posterior.txt"
    write_posterior(path, [(0, 1, 0.5)])
    with pytest.raises(ValueError, match="leaf-pair space"):
        compare_posteriors(
            read_tree_field_posterior(path, n_cells=4),
            read_tree_field_posterior(path, n_cells=8),
        )


# ---------------------------------------------------------------------------
# the ridge
# ---------------------------------------------------------------------------


def test_the_ridge_point_is_the_two_tree_field_densities(tmp_path):
    species, molecules = tmp_path / "s.txt", tmp_path / "m.txt"
    write_posterior(species, [(0, 1, 0.8), (1, 1, 0.8)])
    write_posterior(molecules, [(0, 1, 0.5), (1, 1, 0.5)])
    point = ridge_point(
        read_tree_field_posterior(species, n_cells=2),
        read_tree_field_posterior(molecules, n_cells=2),
    )
    assert point.species == pytest.approx(0.8)
    assert point.molecules == pytest.approx(0.5)
    assert point.product == pytest.approx(0.4)


def test_a_move_along_the_ridge_holds_the_product():
    # ADR-0005, derivation 3: the field's density is the product of the two
    # corrupted rates, so the two trees can trade against each other and leave
    # the field where it was.
    shift = ridge_shift(ridge_point_of(0.8, 0.5), ridge_point_of(0.5, 0.8))
    assert shift.across == pytest.approx(0.0, abs=1e-12)
    assert shift.along != pytest.approx(0.0)


def test_a_move_off_the_ridge_shows_in_the_product():
    shift = ridge_shift(ridge_point_of(0.5, 0.5), ridge_point_of(0.25, 0.25))
    assert shift.along == pytest.approx(0.0, abs=1e-12)
    assert shift.across == pytest.approx(np.log(0.25) - np.log(0.0625))


def test_the_ridge_shift_is_signed_so_the_direction_reads():
    forward = ridge_shift(ridge_point_of(0.8, 0.5), ridge_point_of(0.5, 0.8))
    backward = ridge_shift(ridge_point_of(0.5, 0.8), ridge_point_of(0.8, 0.5))
    assert forward.along == pytest.approx(-backward.along)

