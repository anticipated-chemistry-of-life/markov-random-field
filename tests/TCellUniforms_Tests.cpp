//
// The stream of cell uniforms.
//
// What is asserted here is not that the numbers are random -- SplitMix64's own literature says
// that -- but that this is a *stream per position*: a cell's uniform is a function of the seed, the
// stream, the tree, the iteration and the cell's linear index, and of nothing else. Two properties
// follow, and they are the whole point of the class:
//
//   - the order cells are visited in does not change what any of them gets, so a run gives one
//     chain at any thread count and under any traversal, and
//   - two cells, two iterations, two containers and two seeds do not share a number.
//
// The independence checks below are fixed-seed and therefore deterministic: they pass or they fail,
// and they never flake. The thresholds are the 0.1% critical values, so a change that moves the
// statistic to the edge is a change worth looking at rather than noise.
//

#include "TClique.h"
#include "coretools/Math/TSumLog.h"
#include "omp.h"
#include "random/TCellUniforms.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <set>
#include <vector>

namespace {

constexpr uint64_t SEED = 42;

/// How far a sample of uniforms is from flat, over ten equal bins. Chi-square on 9 degrees of
/// freedom; 27.88 is its 0.1% critical value.
constexpr double CHI_SQUARE_9_DF_AT_0_001 = 27.88;

double chi_square_over_ten_bins(const std::vector<double> &draws) {
	std::array<size_t, 10> counts{};
	for (const double u : draws) {
		const auto bin = static_cast<size_t>(u * 10.0);
		++counts[std::min<size_t>(bin, 9)];
	}
	const double expected = static_cast<double>(draws.size()) / 10.0;
	double statistic      = 0.0;
	for (const size_t count : counts) {
		const double difference = static_cast<double>(count) - expected;
		statistic += difference * difference / expected;
	}
	return statistic;
}

/// Four standard errors of the correlation between two independent samples of `n` draws. The
/// bound shrinks with the sample, so it is computed from the sample rather than picked.
double uncorrelated_bound(size_t n) { return 4.0 / std::sqrt(static_cast<double>(n)); }

/// Pearson correlation of two equally long samples.
double correlation(const std::vector<double> &left, const std::vector<double> &right) {
	const auto n         = static_cast<double>(left.size());
	const double mean_l  = std::accumulate(left.begin(), left.end(), 0.0) / n;
	const double mean_r  = std::accumulate(right.begin(), right.end(), 0.0) / n;
	double covariance    = 0.0;
	double variance_left = 0.0;
	double variance_righ = 0.0;
	for (size_t i = 0; i < left.size(); ++i) {
		const double dl = left[i] - mean_l;
		const double dr = right[i] - mean_r;
		covariance += dl * dr;
		variance_left += dl * dl;
		variance_righ += dr * dr;
	}
	return covariance / std::sqrt(variance_left * variance_righ);
}

/// The uniforms one stream gives the first `n_cells` cells.
std::vector<double> cells_of(const TCellUniforms &uniforms, size_t n_cells) {
	std::vector<double> draws;
	draws.reserve(n_cells);
	for (size_t cell = 0; cell < n_cells; ++cell) { draws.push_back(uniforms.at(cell)); }
	return draws;
}

//--------------------------------------------------------------------------------------
// The number a cell gets, and where it comes from
//--------------------------------------------------------------------------------------

TEST(TCellUniforms, everyDrawLiesInTheHalfOpenUnitInterval) {
	// A uniform of exactly 1 would make a probability of 1 reject, so the top of the range is open.
	const TCellUniforms uniforms(SEED, TCellStream::field, 0);
	for (size_t cell = 0; cell < 100000; ++cell) {
		const double u = uniforms.at(cell);
		EXPECT_GE(u, 0.0) << "cell " << cell;
		EXPECT_LT(u, 1.0) << "cell " << cell;
	}
}

TEST(TCellUniforms, theStreamIsThisFunctionAndNoOther) {
	// The uniforms are pinned to their values. Changing how a position is mixed changes every
	// chain the model has ever run, so it has to be a change somebody made on purpose.
	const TCellUniforms field(SEED, TCellStream::field, 0);
	EXPECT_DOUBLE_EQ(field.at(0), 0.9432402415420365);
	EXPECT_DOUBLE_EQ(field.at(1), 0.18269900947322548);
	EXPECT_DOUBLE_EQ(field.at(12345), 0.6123131166261317);

	const TCellUniforms node_state(SEED, TCellStream::node_state, 3, 1);
	EXPECT_DOUBLE_EQ(node_state.at(7), 0.432817051297131);
}

TEST(TCellUniforms, oneCellAlwaysDrawsTheSameNumber) {
	const TCellUniforms first(SEED, TCellStream::node_state, 7, 1);
	const TCellUniforms second(SEED, TCellStream::node_state, 7, 1);
	for (size_t cell = 0; cell < 1000; ++cell) { EXPECT_EQ(first.at(cell), second.at(cell)); }
}

TEST(TCellUniforms, theOrderCellsAreVisitedInChangesNothing) {
	// This is the property the field update and the node-state walk are built on, and the one the
	// dense and the sparse backend need in order to traverse their storage differently and still
	// run one chain.
	const size_t n_cells = 5000;
	const TCellUniforms uniforms(SEED, TCellStream::field, 3);
	const std::vector<double> in_order = cells_of(uniforms, n_cells);

	std::vector<size_t> visit(n_cells);
	std::iota(visit.begin(), visit.end(), 0);
	std::shuffle(visit.begin(), visit.end(), std::mt19937(12345));

	for (const size_t cell : visit) { EXPECT_EQ(uniforms.at(cell), in_order[cell]); }
}

TEST(TCellUniforms, theThreadCountChangesNothing) {
	const size_t n_cells = 20000;
	const TCellUniforms uniforms(SEED, TCellStream::field, 11);
	const std::vector<double> serial = cells_of(uniforms, n_cells);

	for (const int n_threads : {1, 2, 3, 8}) {
		std::vector<double> parallel(n_cells, -1.0);
#pragma omp parallel for num_threads(n_threads) default(none) shared(uniforms, parallel, n_cells)
		for (size_t cell = 0; cell < n_cells; ++cell) { parallel[cell] = uniforms.at(cell); }
		EXPECT_EQ(parallel, serial) << "with " << n_threads << " threads";
	}
}

//--------------------------------------------------------------------------------------
// What separates one draw from another
//--------------------------------------------------------------------------------------

TEST(TCellUniforms, noTwoCellsOfOneStreamShareANumber) {
	const size_t n_cells            = 200000;
	const std::vector<double> draws = cells_of(TCellUniforms(SEED, TCellStream::field, 0), n_cells);
	const std::set<double> distinct(draws.begin(), draws.end());
	EXPECT_EQ(distinct.size(), n_cells);
}

TEST(TCellUniforms, twoIterationsGiveOneCellTwoNumbers) {
	const size_t n_cells = 5000;
	const std::vector<double> iteration_0 =
	    cells_of(TCellUniforms(SEED, TCellStream::field, 0), n_cells);
	const std::vector<double> iteration_1 =
	    cells_of(TCellUniforms(SEED, TCellStream::field, 1), n_cells);
	for (size_t cell = 0; cell < n_cells; ++cell) {
		EXPECT_NE(iteration_0[cell], iteration_1[cell]) << "cell " << cell;
	}
}

TEST(TCellUniforms, theFieldAndANodeStateDoNotShareTheirCells) {
	// The two containers index their cells in spaces of their own, so cell 5 of one has nothing to
	// do with cell 5 of the other. The stream label is what keeps them apart.
	const size_t n_cells            = 5000;
	const std::vector<double> field = cells_of(TCellUniforms(SEED, TCellStream::field, 4), n_cells);
	const std::vector<double> node_state =
	    cells_of(TCellUniforms(SEED, TCellStream::node_state, 4), n_cells);
	for (size_t cell = 0; cell < n_cells; ++cell) {
		EXPECT_NE(field[cell], node_state[cell]) << "cell " << cell;
	}
	EXPECT_LT(std::abs(correlation(field, node_state)), uncorrelated_bound(n_cells));
}

TEST(TCellUniforms, twoTreesDoNotShareTheirNodeStateCells) {
	const size_t n_cells = 5000;
	const std::vector<double> species =
	    cells_of(TCellUniforms(SEED, TCellStream::node_state, 4, 0), n_cells);
	const std::vector<double> molecules =
	    cells_of(TCellUniforms(SEED, TCellStream::node_state, 4, 1), n_cells);
	for (size_t cell = 0; cell < n_cells; ++cell) {
		EXPECT_NE(species[cell], molecules[cell]) << "cell " << cell;
	}
	EXPECT_LT(std::abs(correlation(species, molecules)), uncorrelated_bound(n_cells));
}

TEST(TCellUniforms, theDrawAChainStartsFromIsNotTheDrawOfItsFirstIteration) {
	const size_t n_cells = 5000;
	const std::vector<double> at_start =
	    cells_of(TCellUniforms(SEED, TCellStream::field_at_start, 0), n_cells);
	const std::vector<double> first_update =
	    cells_of(TCellUniforms(SEED, TCellStream::field, 0), n_cells);
	for (size_t cell = 0; cell < n_cells; ++cell) {
		EXPECT_NE(at_start[cell], first_update[cell]) << "cell " << cell;
	}
}

TEST(TCellUniforms, twoSeedsGiveTwoRuns) {
	const size_t n_cells            = 5000;
	const std::vector<double> run_1 = cells_of(TCellUniforms(1, TCellStream::field, 0), n_cells);
	const std::vector<double> run_2 = cells_of(TCellUniforms(2, TCellStream::field, 0), n_cells);
	for (size_t cell = 0; cell < n_cells; ++cell) {
		EXPECT_NE(run_1[cell], run_2[cell]) << "cell " << cell;
	}
	EXPECT_LT(std::abs(correlation(run_1, run_2)), uncorrelated_bound(n_cells));
}

//--------------------------------------------------------------------------------------
// Independence across cells and across iterations
//--------------------------------------------------------------------------------------

TEST(TCellUniforms, theCellsOfOneIterationAreFlatAndUncorrelated) {
	const size_t n_cells            = 200000;
	const std::vector<double> draws = cells_of(TCellUniforms(SEED, TCellStream::field, 0), n_cells);

	EXPECT_LT(chi_square_over_ten_bins(draws), CHI_SQUARE_9_DF_AT_0_001);

	// Neighbouring cells: the one pair a linear traversal would show up.
	const std::vector<double> cell(draws.begin(), draws.end() - 1);
	const std::vector<double> next_cell(draws.begin() + 1, draws.end());
	EXPECT_LT(std::abs(correlation(cell, next_cell)), uncorrelated_bound(n_cells));

	// Cells a row apart, for a container whose rows are 421 cells wide. A stride is what a clique
	// walks, and it is the pattern a counter that only avalanches its low bits would miss.
	const std::vector<double> strided(draws.begin(), draws.end() - 421);
	const std::vector<double> a_row_on(draws.begin() + 421, draws.end());
	EXPECT_LT(std::abs(correlation(strided, a_row_on)), uncorrelated_bound(n_cells));
}

TEST(TCellUniforms, theIterationsOfOneCellAreFlatAndUncorrelated) {
	const size_t n_iterations = 200000;
	const size_t cell         = 12345;
	std::vector<double> draws;
	draws.reserve(n_iterations);
	for (size_t iteration = 0; iteration < n_iterations; ++iteration) {
		draws.push_back(TCellUniforms(SEED, TCellStream::field, iteration).at(cell));
	}

	EXPECT_LT(chi_square_over_ten_bins(draws), CHI_SQUARE_9_DF_AT_0_001);

	const std::vector<double> iteration(draws.begin(), draws.end() - 1);
	const std::vector<double> next(draws.begin() + 1, draws.end());
	EXPECT_LT(std::abs(correlation(iteration, next)), uncorrelated_bound(n_iterations));
}

//--------------------------------------------------------------------------------------
// The draw the uniform feeds
//--------------------------------------------------------------------------------------

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
