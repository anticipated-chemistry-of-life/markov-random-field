//
// The two-state draw.
//
// The draw is the last step of every cell update: a probability, or a pair of log probabilities,
// and the cell's own uniform decide the state. What is asserted here is that the uniform decides
// and nothing else does, that the three overloads agree, and that the tails need no branch.
//
// The last test joins the draw to the stream it is fed from (random/TCellUniforms.h), because the
// two together are what a cell update actually does.
//

#include "coretools/Math/TSumLog.h"
#include "coretools/Types/probability.h"
#include "random/TCellUniforms.h"
#include "random/two_state_draw.h"
#include "gtest/gtest.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace {

constexpr uint64_t SEED = 42;

using two_state_draw::sample;

TEST(Sample, theUniformDecidesAndNothingElseDoes) {
	// P(state 1) = 0.75, so the draw turns over between a uniform of 0.7 and one of 0.8.
	const double log_prob_0 = std::log(0.25);
	const double log_prob_1 = std::log(0.75);
	EXPECT_TRUE(sample(log_prob_0, log_prob_1, 0.70));
	EXPECT_FALSE(sample(log_prob_0, log_prob_1, 0.80));
	EXPECT_TRUE(sample(log_prob_0, log_prob_1, 0.0));
	EXPECT_FALSE(sample(log_prob_0, log_prob_1, 0.9999));
}

TEST(Sample, theSumOfLogsAndThePairOfLogsAgree) {
	std::array<coretools::TSumLogProbability, 2> sum_log;
	sum_log[0].add(0.25);
	sum_log[1].add(0.75);
	for (const double uniform : {0.0, 0.1, 0.5, 0.7499, 0.7501, 0.9999}) {
		EXPECT_EQ(sample(sum_log, uniform), sample(std::log(0.25), std::log(0.75), uniform))
		    << "at " << uniform;
	}
}

TEST(Sample, aProbabilityDrawsAtTheRateItNames) {
	EXPECT_TRUE(sample(coretools::Probability(0.75), 0.70));
	EXPECT_FALSE(sample(coretools::Probability(0.75), 0.80));
	EXPECT_FALSE(sample(coretools::Probability(0.0), 0.0));
	EXPECT_TRUE(sample(coretools::Probability(1.0), 0.9999999));
}

TEST(Sample, theTailsNeedNoBranch) {
	// Log odds far from zero: the exponential saturates, and the answer is the one the tail asks
	// for rather than an overflow.
	EXPECT_TRUE(sample(0.0, 1000.0, 0.999999));
	EXPECT_FALSE(sample(1000.0, 0.0, 0.000001));
	// Two impossible states leave the log odds undefined, which reads as state 0.
	const double minus_infinity = -std::numeric_limits<double>::infinity();
	EXPECT_FALSE(sample(minus_infinity, minus_infinity, 0.5));
}

TEST(Sample, theDrawKeepsItsProbabilityOverTheStream) {
	// The two halves together: uniforms that come from a cell's position, fed to the draw, give
	// state 1 at the rate the log odds name.
	const TCellUniforms uniforms(SEED, TCellStream::field, 0);
	const double log_prob_0 = std::log(0.3);
	const double log_prob_1 = std::log(0.7);

	size_t n_ones        = 0;
	const size_t n_draws = 200000;
	for (size_t cell = 0; cell < n_draws; ++cell) {
		n_ones += static_cast<size_t>(sample(log_prob_0, log_prob_1, uniforms.at(cell)));
	}
	const double fraction = static_cast<double>(n_ones) / static_cast<double>(n_draws);
	// Four standard errors of a binomial with p = 0.7 over 200000 draws is about 0.004.
	EXPECT_NEAR(fraction, 0.7, 0.004);
}

} // namespace
