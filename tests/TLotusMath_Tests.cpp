// Properties of the LOTUS reporting model.
//
// The first four mirror the research-effort checks in model_validation/tests/test_independent.py,
// which until now could only be asked of the Python reference. The rest cover ground the reference
// cannot reach at all: it simulates with a single scalar gamma, while the C++ infers one per
// tree.

#include "lotus/TLotusMath.h"
#include "gtest/gtest.h"

#include <cmath>
#include <vector>

namespace {

using lotus_math::TReportingModel;

IndexArray cell(size_t first, size_t second) { return IndexArray{first, second}; }

/// Two dimensions of three leaves each, with paper counts that differ within and across dimensions.
TReportingModel two_dimensions(double gamma_a, double gamma_b, double error_rate = 0.01) {
	return {{gamma_a, gamma_b}, error_rate, {{0, 1, 7}, {0, 3, 40}}};
}

} // namespace

// --------------------------------------------------------------------------
// Research effort
// --------------------------------------------------------------------------

TEST(LotusMath, research_effort_uses_log_paper_counts_not_raw_ones) {
	// The trap the validation README calls out: driving research effort from raw counts is
	// indistinguishable from an inference bug. log(count + 1) is the transform, and it now lives
	// here rather than in whatever reads the counts off disk.
	// A small gamma keeps both routes away from saturation, where they would agree at ~1 whatever
	// the transform: 40 papers give 0.17 through the log and 0.86 through the raw count.
	const double gamma = 0.05;
	const TReportingModel model({gamma}, 0.01, {{40}});

	const double from_log = 1.0 - std::exp(-gamma * std::log(41.0));
	const double from_raw = 1.0 - std::exp(-gamma * 40.0);

	EXPECT_DOUBLE_EQ(model.research_effort(IndexArray{0}), from_log);
	EXPECT_GT(std::abs(from_log - from_raw), 0.5) << "the two would be indistinguishable here";
}

TEST(LotusMath, a_leaf_with_no_papers_is_never_reported) {
	// occurrence_count(0) = log(1) = 0, so the factor is 1 - exp(0) = 0 exactly, and the product
	// across dimensions collapses. A pair nobody studied contributes no evidence of presence.
	const auto model = two_dimensions(1.5, 1.5);

	EXPECT_DOUBLE_EQ(model.research_effort(cell(0, 2)), 0.0);
	EXPECT_DOUBLE_EQ(model.research_effort(cell(2, 0)), 0.0);
	EXPECT_DOUBLE_EQ(model.research_effort(cell(0, 0)), 0.0);
	EXPECT_GT(model.research_effort(cell(2, 2)), 0.0);
}

TEST(LotusMath, research_effort_is_a_probability) {
	for (const double gamma : {0.01, 0.5, 1.1, 5.0, 50.0}) {
		const auto model = two_dimensions(gamma, gamma);
		for (size_t first = 0; first < 3; ++first) {
			for (size_t second = 0; second < 3; ++second) {
				const double effort = model.research_effort(cell(first, second));
				EXPECT_GE(effort, 0.0) << "gamma " << gamma;
				EXPECT_LE(effort, 1.0) << "gamma " << gamma;
			}
		}
	}
}

TEST(LotusMath, research_effort_rises_with_papers_and_with_gamma) {
	const auto model = two_dimensions(1.0, 1.0);
	// leaf order within a dimension is (0, 1, 7) and (0, 3, 40) papers
	EXPECT_LT(model.research_effort(cell(1, 1)), model.research_effort(cell(2, 1)));
	EXPECT_LT(model.research_effort(cell(1, 1)), model.research_effort(cell(1, 2)));

	const auto keener = two_dimensions(3.0, 3.0);
	EXPECT_LT(model.research_effort(cell(2, 2)), keener.research_effort(cell(2, 2)));
}

TEST(LotusMath, effort_saturates_as_gamma_grows) {
	const auto model = two_dimensions(500.0, 500.0);
	EXPECT_NEAR(model.research_effort(cell(2, 2)), 1.0, 1e-9);
}

// --------------------------------------------------------------------------
// One gamma per tree
// --------------------------------------------------------------------------

TEST(LotusMath, each_dimension_uses_its_own_gamma) {
	// The Python reference simulates with a single scalar gamma for both trees, so this asymmetry
	// is invisible to it. Here the two dimensions are given different rates and identical counts,
	// which must produce different factors.
	const TReportingModel model({0.2, 2.0}, 0.01, {{9}, {9}});

	const double slow = 1.0 - std::exp(-0.2 * std::log(10.0));
	const double fast = 1.0 - std::exp(-2.0 * std::log(10.0));
	ASSERT_NE(slow, fast);

	EXPECT_DOUBLE_EQ(model.research_effort(IndexArray{0, 0}), slow * fast);
}

TEST(LotusMath, effort_is_the_product_across_dimensions) {
	const auto model = two_dimensions(0.9, 1.4);
	const TReportingModel only_first({0.9}, 0.01, {{0, 1, 7}});
	const TReportingModel only_second({1.4}, 0.01, {{0, 3, 40}});

	EXPECT_DOUBLE_EQ(model.research_effort(cell(2, 1)),
	                 only_first.research_effort(IndexArray{2}) *
	                     only_second.research_effort(IndexArray{1}));
}

TEST(LotusMath, a_gamma_per_dimension_is_required) {
	EXPECT_THROW((TReportingModel({1.0}, 0.01, {{1}, {2}})), std::invalid_argument);
	EXPECT_THROW((TReportingModel({1.0, 2.0}, 0.01, {{1}})), std::invalid_argument);
}

// --------------------------------------------------------------------------
// The emission
// --------------------------------------------------------------------------

TEST(LotusMath, a_present_pair_is_reported_at_its_research_effort) {
	const auto model    = two_dimensions(1.1, 1.1);
	const auto position = cell(2, 2);
	const double effort = model.research_effort(position);

	EXPECT_DOUBLE_EQ(model.probability(true, true, position), effort);
	EXPECT_DOUBLE_EQ(model.probability(true, false, position), 1.0 - effort);
}

TEST(LotusMath, an_absent_pair_is_reported_at_the_flat_error_rate) {
	const double error_rate = 0.03;
	const auto model        = two_dimensions(1.1, 1.1, error_rate);

	EXPECT_DOUBLE_EQ(model.probability(false, true, cell(2, 2)), error_rate);
	EXPECT_DOUBLE_EQ(model.probability(false, false, cell(2, 2)), 1.0 - error_rate);
}

TEST(LotusMath, the_absent_case_does_not_depend_on_position) {
	// This is what lets the likelihood sweep answer absent cells without converting a linear index,
	// and what makes the bulk term for never-stored cells a single constant.
	const auto model = two_dimensions(1.1, 1.1, 0.03);
	for (size_t first = 0; first < 3; ++first) {
		for (size_t second = 0; second < 3; ++second) {
			EXPECT_DOUBLE_EQ(model.probability(false, true, cell(first, second)),
			                 model.probability_absent(true));
			EXPECT_DOUBLE_EQ(model.probability(false, false, cell(first, second)),
			                 model.probability_absent(false));
		}
	}
}

TEST(LotusMath, both_states_give_a_distribution_over_L) {
	const auto model = two_dimensions(1.1, 0.4, 0.03);
	for (size_t first = 0; first < 3; ++first) {
		for (size_t second = 0; second < 3; ++second) {
			for (const bool x : {false, true}) {
				const double reported = model.probability(x, true, cell(first, second));
				const double missing  = model.probability(x, false, cell(first, second));
				EXPECT_GE(reported, 0.0);
				EXPECT_GE(missing, 0.0);
				EXPECT_NEAR(reported + missing, 1.0, 1e-12)
				    << "x " << x << " at (" << first << ", " << second << ")";
			}
		}
	}
}

TEST(LotusMath, an_unstudied_present_pair_looks_absent_but_is_not_evidence) {
	// The asymmetry the source exists for: with no papers, a present pair is never reported, so a
	// missing record says nothing. A well-studied pair reports almost always, so a missing record
	// there is real evidence.
	const auto model = two_dimensions(2.0, 2.0);

	EXPECT_DOUBLE_EQ(model.probability(true, false, cell(0, 0)), 1.0);
	EXPECT_LT(model.probability(true, false, cell(2, 2)), 0.6);
}
