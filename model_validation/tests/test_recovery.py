"""Scores for how well inference recovers the field and the detection parameters."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.independent import mass_spec, recovery


def test_the_band_around_nine_tenths_of_a_hundred_covers_84_to_95():
    # Worked by hand from Binomial(100, 0.9): the central 95% mass is 84..95.
    low, high = recovery.coverage_band(100)
    assert (low, high) == (0.84, 0.95)


def test_a_single_trial_cannot_fall_outside_its_band():
    assert recovery.coverage_band(1) == (0.0, 1.0)


def test_a_coverage_passes_only_inside_the_band():
    assert recovery.coverage_passes(0.9, 100)
    assert not recovery.coverage_passes(0.60, 100)
    assert not recovery.coverage_passes(1.0, 100)  # over-wide intervals fail too


def test_a_filter_is_covered_when_its_truth_lies_in_the_90_percent_interval():
    # Draws 0..99 per column: the 5%..95% interval is about [4.95, 94.05].
    draws = np.tile(np.arange(100.0)[:, None], (1, 3))
    truth = np.array([50.0, 2.0, 99.0])
    assert recovery.is_covered(draws, truth).tolist() == [True, False, False]


def test_filter_coverage_is_stratified_by_the_runs_a_molecule_is_present_in():
    # Four molecules present in 0, 1, 1 and 3 runs. The molecule present in 3 is
    # missed; the one present in 0 and the two in 1 are covered.
    draws = np.tile(np.arange(100.0)[:, None], (1, 4))
    truth = np.array([50.0, 40.0, 60.0, 0.0])
    runs = np.array([0, 1, 1, 3])
    rows = {r.label: r for r in recovery.coverage_by_runs(draws, truth, runs)}
    assert rows["0"].n == 1 and rows["0"].coverage == 1.0
    assert rows["1"].n == 2 and rows["1"].coverage == 1.0
    assert rows["3-9"].n == 1 and rows["3-9"].coverage == 0.0
    assert "2" not in rows  # empty strata are not reported


def test_auc_matches_the_textbook_example():
    # Four pairs of (positive, negative), three ranked correctly: 0.75.
    score = np.array([0.1, 0.4, 0.35, 0.8])
    truth = np.array([0, 0, 1, 1])
    assert recovery.auc(truth, score) == 0.75


def test_auc_counts_ties_as_half():
    assert recovery.auc(np.array([0, 1]), np.array([0.5, 0.5])) == 0.5


def test_auc_is_undefined_with_one_class():
    assert np.isnan(recovery.auc(np.array([0, 0, 0]), np.array([0.1, 0.2, 0.3])))


def test_the_reliability_table_pairs_mean_prediction_with_observed_frequency():
    # Two cells in the 0.0-0.1 bin (one true), three in the 0.9-1.0 bin (all true).
    predicted = np.array([0.02, 0.08, 0.90, 0.95, 1.0])
    truth = np.array([0, 1, 1, 1, 1])
    table = recovery.reliability_table(truth, predicted, n_bins=10)
    assert [(r.n, r.observed) for r in table] == [(2, 0.5), (3, 1.0)]
    assert table[0].predicted == pytest.approx(0.05)
    assert table[1].predicted == pytest.approx(0.95)
    assert table[0].low == 0.0 and table[1].high == 1.0


def _field_frames():
    names = [(s, m) for s in ("a", "b") for m in ("x", "y")]
    truth = pd.DataFrame(
        {
            "species": [s for s, _ in names],
            "molecules": [m for _, m in names],
            "Y_state": [1, 0, 1, 0],
            # The truth file carries the same column name as the posterior.
            "fraction_of_one": [1.0, 0.0, 1.0, 0.0],
        }
    )
    # The inference file lists only cells that carry a posterior; (b, y) is absent.
    posterior = pd.DataFrame(
        {
            "species": ["a", "a", "b"],
            "molecules": ["x", "y", "x"],
            "fraction_of_one": [0.9, 0.2, 0.6],
        }
    )
    return truth, posterior


def test_field_scores_split_species_with_runs_from_species_without():
    truth, posterior = _field_frames()
    scores = recovery.field_scores(truth, posterior, species_with_runs={"a"})
    assert scores["with_runs"].n == 2 and scores["with_runs"].auc == 1.0
    # Species b: cell (b,x) truly 1 at 0.6, cell (b,y) truly 0 and absent => 0.
    assert scores["without_runs"].n == 2 and scores["without_runs"].auc == 1.0
    # A cell missing from the posterior was never a one: its probability is 0.
    assert scores["without_runs"].reliability[0].n == 1


def test_a_molecule_is_present_in_every_run_of_every_species_that_carries_it():
    # Species a carries x and has 2 runs; species b carries x too but has no runs.
    truth, _ = _field_frames()
    runs = recovery.runs_present(
        truth, species_with_runs={"a"}, runs_per_species=2, molecules=["x", "y"]
    )
    assert runs.tolist() == [2, 0]


# --------------------------------------------------------------------------
# The assignment posterior
# --------------------------------------------------------------------------

UNKNOWN = mass_spec.UNKNOWN
KEY = recovery.FEATURE_KEY


def _assignment(rows):
    """`rows`: (feature, molecule, value) for one species and extract."""
    return pd.DataFrame(
        [("sp", "e", f, m, v) for f, m, v in rows], columns=KEY + ["molecule", "fraction"]
    )


def _data(rows):
    """The data file lists candidates only, never the unknown molecule."""
    return pd.DataFrame(
        [("sp", "e", f, m) for f, m in rows], columns=KEY + ["molecule"]
    )


# Feature f1 is molecule A out of {A, B}; f2 is unknown, candidates {A, B};
# f3 is molecule B out of {A, B}.
TRUTH = _assignment(
    [
        ("f1", "A", 1), ("f1", "B", 0), ("f1", UNKNOWN, 0),
        ("f2", "A", 0), ("f2", "B", 0), ("f2", UNKNOWN, 1),
        ("f3", "A", 0), ("f3", "B", 1), ("f3", UNKNOWN, 0),
    ]
)  # fmt: skip
DATA = _data([(f, m) for f in ("f1", "f2", "f3") for m in ("A", "B")])


def _posterior(f1, f2, f3):
    rows = []
    for f, (a, b, u) in zip(("f1", "f2", "f3"), (f1, f2, f3)):
        rows += [(f, "A", a), (f, "B", b), (f, UNKNOWN, u)]
    return _assignment(rows)


def test_assignment_scores_are_worked_by_hand():
    # f1: mode A (right), mass on truth 0.6, unknown 0.1.
    # f2: mode unknown (right), mass on truth 0.5, unknown 0.5.
    # f3: mode A (wrong), mass on truth 0.25, unknown 0.25.
    posterior = _posterior((0.6, 0.3, 0.1), (0.25, 0.25, 0.5), (0.5, 0.25, 0.25))
    s = recovery.assignment_scores(TRUTH, posterior)
    assert s.n == 3
    assert s.accuracy == pytest.approx(2 / 3)
    assert s.mean_on_true_assignment == pytest.approx((0.6 + 0.5 + 0.25) / 3)
    # Unknown posterior: the unknown feature f2 alone; molecule features f1, f3.
    assert s.mean_unknown_of_unknown_features == pytest.approx(0.5)
    assert s.n_unknown_features == 1
    assert s.mean_unknown_of_molecule_features == pytest.approx((0.1 + 0.25) / 2)
    assert s.n_molecule_features == 2


def test_a_perfect_posterior_scores_one_everywhere_it_should():
    s = recovery.assignment_scores(TRUTH, TRUTH.copy())
    assert (s.accuracy, s.mean_on_true_assignment) == (1.0, 1.0)
    assert s.mean_unknown_of_unknown_features == 1.0
    assert s.mean_unknown_of_molecule_features == 0.0


def test_a_scenario_without_unknown_features_has_no_unknown_feature_score():
    truth = TRUTH[TRUTH["feature"] != "f2"]
    posterior = _posterior((0.6, 0.3, 0.1), (0, 0, 1), (0.5, 0.25, 0.25))
    posterior = posterior[posterior["feature"] != "f2"]
    s = recovery.assignment_scores(truth, posterior)
    assert s.n_unknown_features == 0
    assert np.isnan(s.mean_unknown_of_unknown_features)


GOOD = _posterior((0.6, 0.3, 0.1), (0.25, 0.25, 0.5), (0.5, 0.25, 0.25))


def test_a_well_formed_posterior_passes_the_format_check():
    recovery.check_assignment_format(DATA, GOOD)


def test_a_missing_feature_is_named():
    bad = GOOD[GOOD["feature"] != "f2"]
    with pytest.raises(recovery.AssignmentFormatError, match="f2"):
        recovery.check_assignment_format(DATA, bad)


def test_a_missing_candidate_is_named():
    bad = GOOD[~((GOOD["feature"] == "f3") & (GOOD["molecule"] == "B"))]
    with pytest.raises(recovery.AssignmentFormatError, match=r"f3.*B"):
        recovery.check_assignment_format(DATA, bad)


def test_a_missing_unknown_row_is_named():
    bad = GOOD[~((GOOD["feature"] == "f1") & (GOOD["molecule"] == UNKNOWN))]
    with pytest.raises(recovery.AssignmentFormatError, match=r"f1.*<unknown>"):
        recovery.check_assignment_format(DATA, bad)


def test_an_extra_row_is_named():
    extra = _assignment([("f1", "C", 0.0)])
    with pytest.raises(recovery.AssignmentFormatError, match=r"f1.*C"):
        recovery.check_assignment_format(DATA, pd.concat([GOOD, extra]))


def test_a_feature_that_does_not_sum_to_one_is_named():
    bad = _posterior((0.6, 0.3, 0.1), (0.25, 0.25, 0.4), (0.5, 0.25, 0.25))
    with pytest.raises(recovery.AssignmentFormatError, match=r"f2.*sum"):
        recovery.check_assignment_format(DATA, bad)


def test_a_fraction_outside_zero_one_is_named():
    bad = _posterior((0.6, 0.3, 0.1), (0.25, 0.25, 0.5), (1.5, -0.5, 0.0))
    with pytest.raises(recovery.AssignmentFormatError, match=r"f3"):
        recovery.check_assignment_format(DATA, bad)


def test_a_fraction_that_is_not_a_number_is_named():
    bad = GOOD.astype({"fraction": object})
    bad.loc[bad.index[7], "fraction"] = "abc"
    with pytest.raises(recovery.AssignmentFormatError, match=r"f3"):
        recovery.check_assignment_format(DATA, bad)


def test_a_tied_mode_goes_to_the_first_row_the_posterior_lists():
    # f1 ties A and B at 0.45; the posterior lists A first and A is true.
    tied = _posterior((0.45, 0.45, 0.1), (0, 0, 1), (0, 1, 0))
    reordered_truth = TRUTH.iloc[[1, 0, *range(2, len(TRUTH))]]
    assert recovery.assignment_scores(reordered_truth, tied).accuracy == 1.0
