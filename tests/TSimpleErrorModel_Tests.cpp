#include "TSimpleErrorModelMath.h"
#include "coretools/Main/TRandomGenerator.h"
#include "storages/y_storage/TStorageYMatrix.h"
#include "gtest/gtest.h"
#include <cmath>
#include <cstddef>
#include <vector>

namespace {

using namespace simple_error_model;

/// A matrix laid out as a single row {1, N}, so that the linear index of a cell equals its column.
/// Same convention as TStorageY_Tests.
TStorageYMatrix make_row(size_t n_cols, const std::vector<size_t> &ones) {
	TStorageYMatrix m(1000, {1, n_cols});
	for (const auto &i : ones) { m.insert_one(i); }
	return m;
}

/// The definition the closed form is a shortcut for: sum over every cell of log P(D | Y).
/// Deliberately walks all cells (including the absent ones) so it shares no code with
/// count_disagreements.
double brute_force_log_likelihood(const TStorageYMatrix &Y, const TStorageYMatrix &D, double eps) {
	double sum = 0.0;
	for (size_t i = 0; i < Y.total_size_of_container_space(); ++i) {
		sum += std::log(probability_of_D_given_Y(Y.is_one(i), D.is_one(i), eps));
	}
	return sum;
}

} // namespace

//-----------------------------------
// probability_of_D_given_Y / probabilities_for_both_Y_states
//-----------------------------------

TEST(SimpleErrorModel_Tests, probability_agrees_is_one_minus_epsilon) {
	constexpr double eps = 0.25;
	EXPECT_DOUBLE_EQ(probability_of_D_given_Y(true, true, eps), 0.75);
	EXPECT_DOUBLE_EQ(probability_of_D_given_Y(false, false, eps), 0.75);
}

TEST(SimpleErrorModel_Tests, probability_disagrees_is_epsilon) {
	constexpr double eps = 0.25;
	EXPECT_DOUBLE_EQ(probability_of_D_given_Y(true, false, eps), 0.25);
	EXPECT_DOUBLE_EQ(probability_of_D_given_Y(false, true, eps), 0.25);
}

TEST(SimpleErrorModel_Tests, probabilities_for_both_states_match_the_scalar_form) {
	constexpr double eps = 0.3;
	for (const bool d : {false, true}) {
		std::array<double, 2> prob{};
		probabilities_for_both_Y_states(d, eps, prob);
		EXPECT_DOUBLE_EQ(prob[0], probability_of_D_given_Y(false, d, eps)) << "observed D = " << d;
		EXPECT_DOUBLE_EQ(prob[1], probability_of_D_given_Y(true, d, eps)) << "observed D = " << d;
		// exactly one of the two Y states agrees with D, so the pair is always {eps, 1 - eps}
		EXPECT_DOUBLE_EQ(prob[0] + prob[1], 1.0) << "observed D = " << d;
	}
}

//-----------------------------------
// log_likelihood_from_counts
//-----------------------------------

TEST(SimpleErrorModel_Tests, log_likelihood_all_agree) {
	EXPECT_NEAR(log_likelihood_from_counts(10, 0, 0.1), 10.0 * std::log(0.9), 1e-12);
}

TEST(SimpleErrorModel_Tests, log_likelihood_all_disagree) {
	EXPECT_NEAR(log_likelihood_from_counts(10, 10, 0.1), 10.0 * std::log(0.1), 1e-12);
}

TEST(SimpleErrorModel_Tests, log_likelihood_mixed) {
	EXPECT_NEAR(log_likelihood_from_counts(10, 3, 0.1), 7.0 * std::log(0.9) + 3.0 * std::log(0.1),
	            1e-12);
}

TEST(SimpleErrorModel_Tests, log_likelihood_more_disagreements_than_cells_throws) {
	EXPECT_ANY_THROW((void)log_likelihood_from_counts(4, 5, 0.1));
}

/// Pins the O(1) closed form to the per-cell definition it is derived from. If either side is
/// refactored (log1p, a different parameterisation) this is the test that catches a drift.
TEST(SimpleErrorModel_Tests, log_likelihood_matches_cellwise_product) {
	constexpr double eps    = 0.2;
	const auto Y            = make_row(6, {0, 2, 5});
	const auto D            = make_row(6, {0, 3});
	// Y = 1 0 1 0 0 1, D = 1 0 0 1 0 0 -> disagreements at cells 2, 3 and 5
	const size_t n_disagree = count_disagreements(Y, D);
	ASSERT_EQ(n_disagree, 3u);
	EXPECT_NEAR(log_likelihood_from_counts(Y.total_size_of_container_space(), n_disagree, eps),
	            brute_force_log_likelihood(Y, D, eps), 1e-12);
}

/// The likelihood must be maximised at the empirical disagreement fraction. This is the cheapest
/// check that catches a swapped eps / (1 - eps), which is otherwise almost invisible.
TEST(SimpleErrorModel_Tests, log_likelihood_maximized_at_empirical_fraction) {
	constexpr size_t total      = 100;
	constexpr size_t n_disagree = 30;
	const double mle            = static_cast<double>(n_disagree) / static_cast<double>(total);
	const double ll_at_mle      = log_likelihood_from_counts(total, n_disagree, mle);

	for (size_t step = 1; step < 100; ++step) {
		const double eps = static_cast<double>(step) / 100.0;
		EXPECT_GE(ll_at_mle, log_likelihood_from_counts(total, n_disagree, eps) - 1e-12)
		    << "epsilon = " << eps;
	}
}

TEST(SimpleErrorModel_Tests, log_likelihood_ratio_matches_the_difference) {
	constexpr size_t total      = 50;
	constexpr size_t n_disagree = 12;
	EXPECT_NEAR(log_likelihood_ratio(total, n_disagree, 0.1, 0.25),
	            log_likelihood_from_counts(total, n_disagree, 0.25) -
	                log_likelihood_from_counts(total, n_disagree, 0.1),
	            1e-12);
}

//-----------------------------------
// count_disagreements (merge-join over two sparse matrices)
//-----------------------------------

TEST(SimpleErrorModel_Tests, count_disagreements_both_empty) {
	EXPECT_EQ(count_disagreements(make_row(5, {}), make_row(5, {})), 0u);
}

TEST(SimpleErrorModel_Tests, count_disagreements_identical_ones) {
	EXPECT_EQ(count_disagreements(make_row(5, {1, 3}), make_row(5, {1, 3})), 0u);
}

TEST(SimpleErrorModel_Tests, count_disagreements_disjoint_ones) {
	// Y ones at {0, 3}, D ones at {1, 3} -> cells 0 and 1 disagree, cell 3 agrees
	EXPECT_EQ(count_disagreements(make_row(5, {0, 3}), make_row(5, {1, 3})), 2u);
}

/// A cell stored with state 0 is not the same as a cell that is one, but it *is* the same as a cell
/// that was never stored. Confusing "is stored" with "is one" is the classic defect in a
/// merge-join over sparse matrices, so it gets its own test.
TEST(SimpleErrorModel_Tests, count_disagreements_stored_zero_is_an_agreement) {
	TStorageYMatrix Y(1000, {1, 5});
	TStorageYMatrix D(1000, {1, 5});
	D.insert_zero(2); // stored in D, absent in Y, both read as state 0
	EXPECT_EQ(count_disagreements(Y, D), 0u);
}

TEST(SimpleErrorModel_Tests, count_disagreements_D_full_Y_empty) {
	const auto D = make_row(4, {0, 1, 2, 3});
	EXPECT_EQ(count_disagreements(make_row(4, {}), D), D.total_size_of_container_space());
}

TEST(SimpleErrorModel_Tests, count_disagreements_spans_multiple_rows) {
	// A 3x4 layout exercises the cursor's row-skipping: linear index = row * 4 + col.
	TStorageYMatrix Y(1000, {3, 4});
	TStorageYMatrix D(1000, {3, 4});
	Y.insert_one(1);                          // row 0
	Y.insert_one(7);                          // row 1
	D.insert_one(7);                          // agrees
	D.insert_one(11);                         // row 2, Y absent there
	EXPECT_EQ(count_disagreements(Y, D), 2u); // cells 1 and 11
}

TEST(SimpleErrorModel_Tests, count_disagreements_dimension_mismatch_throws) {
	TStorageYMatrix Y(1000, {2, 3});
	TStorageYMatrix D(1000, {3, 2}); // same number of cells, different layout
	EXPECT_ANY_THROW((void)count_disagreements(Y, D));
}

//-----------------------------------
// draw_D_given_Y
//-----------------------------------

TEST(SimpleErrorModel_Tests, draw_D_given_Y_flips_at_rate_epsilon) {
	coretools::instances::randomGenerator().setSeed(42, true);
	constexpr size_t n_draws = 100000;
	constexpr double eps     = 0.3;

	for (const bool y : {false, true}) {
		size_t n_flipped = 0;
		for (size_t i = 0; i < n_draws; ++i) {
			if (draw_D_given_Y(y, eps) != y) { ++n_flipped; }
		}
		const double fraction = static_cast<double>(n_flipped) / static_cast<double>(n_draws);
		EXPECT_NEAR(fraction, eps, 0.01) << "latent state Y = " << y;
	}
}

TEST(SimpleErrorModel_Tests, draw_D_given_Y_is_deterministic_at_the_extremes) {
	coretools::instances::randomGenerator().setSeed(42, true);
	for (size_t i = 0; i < 100; ++i) {
		// pickOneOfTwo uses a strict less-than against the drawn uniform, so eps = 0 never flips
		EXPECT_TRUE(draw_D_given_Y(true, 0.0));
		EXPECT_FALSE(draw_D_given_Y(false, 0.0));
	}
}

/// End-to-end check that the simulation and the likelihood share one parameterisation: drawing D
/// from Y at rate eps must produce a disagreement fraction of about eps as measured by the very
/// function the likelihood uses.
TEST(SimpleErrorModel_Tests, simulated_disagreement_fraction_matches_epsilon) {
	coretools::instances::randomGenerator().setSeed(42, true);
	constexpr double eps = 0.2;

	// 50k cells keeps the sampling error at sd = sqrt(eps(1-eps)/n) ~ 0.0018, so the 0.01 bound
	// below sits at ~5.5 sd and stays comfortable even if the seed changes.
	TStorageYMatrix Y(1000, {200, 250});
	TStorageYMatrix D(1000, {200, 250});
	const size_t total = Y.total_size_of_container_space();

	for (size_t i = 0; i < total; ++i) {
		const bool y = coretools::instances::randomGenerator().pickOneOfTwo();
		if (y) { Y.insert_one(i); }
		if (draw_D_given_Y(y, eps)) { D.insert_one(i); }
	}

	const double fraction =
	    static_cast<double>(count_disagreements(Y, D)) / static_cast<double>(total);
	EXPECT_NEAR(fraction, eps, 0.01);
}
