// Properties of one clique's two-state process, discretised onto the bin grid.
//
// The C++ half of the transition-matrix checks in model_validation/tests/test_independent.py. The
// Python reference derives the same quantities from the closed form; these ask the shipped code
// directly, which until TTransitionGrid was pulled out of TClique took a full MCMC run to reach.

#include "process/TBinGrid.h"
#include "process/TTransitionGrid.h"
#include "gtest/gtest.h"

#include <armadillo>
#include <cmath>
#include <vector>

namespace {

constexpr size_t N_BINS = 10;

const std::vector<double> ALPHAS = {0.05, 0.3, 0.5, 0.72, 0.95};
const std::vector<double> NUS    = {0.05, 0.6, 3.0};

/// P(child = to | parent = from) after time `t`, from a *numerical* matrix exponential of the
/// generator.
///
/// This is deliberately the other route to the same number: TTransitionGrid evaluates the closed
/// form and walks the grid by multiplication, while this calls armadillo's expmat at `t` directly.
/// Keeping expmat here rather than in the implementation is what preserves the cross-check after
/// the closed-form swap -- had both sides used the same formula, the comparison would be vacuous.
double numerical_expm(double alpha, double nu, double t, bool from, bool to) {
	arma::mat generator(2, 2);
	generator(0, 0) = (-alpha) * nu;
	generator(0, 1) = alpha * nu;
	generator(1, 0) = (1.0 - alpha) * nu;
	generator(1, 1) = (alpha - 1.0) * nu;

	// expmat returns a lazy expression; materialise it before indexing.
	const arma::mat exponential = arma::expmat(generator * t);
	return exponential(static_cast<arma::uword>(from), static_cast<arma::uword>(to));
}

} // namespace

// --------------------------------------------------------------------------
// Every bin's matrix is a transition matrix
// --------------------------------------------------------------------------

TEST(TransitionGrid, rows_are_distributions) {
	const TBinGrid bins(N_BINS);
	for (const double alpha : ALPHAS) {
		for (const double nu : NUS) {
			const TTransitionGrid grid(alpha, nu, bins);
			for (size_t bin = 0; bin < N_BINS; ++bin) {
				for (const bool from : {false, true}) {
					const double to_zero = grid.probability(bin, from, false);
					const double to_one  = grid.probability(bin, from, true);
					EXPECT_GE(to_zero, 0.0);
					EXPECT_GE(to_one, 0.0);
					EXPECT_NEAR(to_zero + to_one, 1.0, 1e-12)
					    << "alpha " << alpha << " nu " << nu << " bin " << bin;
				}
			}
		}
	}
}

TEST(TransitionGrid, the_stationary_distribution_is_preserved) {
	const TBinGrid bins(N_BINS);
	for (const double alpha : ALPHAS) {
		for (const double nu : NUS) {
			const TTransitionGrid grid(alpha, nu, bins);
			for (size_t bin = 0; bin < N_BINS; ++bin) {
				// (1-alpha, alpha) P = (1-alpha, alpha)
				const double one = (1.0 - alpha) * grid.probability(bin, false, true) +
				                   alpha * grid.probability(bin, true, true);
				EXPECT_NEAR(one, alpha, 1e-12)
				    << "alpha " << alpha << " nu " << nu << " bin " << bin;
			}
		}
	}
}

TEST(TransitionGrid, a_long_branch_forgets_the_parent) {
	// Still below the stationary short-circuit, so this exercises the real recursion.
	const TBinGrid bins(N_BINS);
	const double alpha = 0.3;
	const TTransitionGrid grid(alpha, 24.0, bins);

	const size_t longest = N_BINS - 1;
	EXPECT_NEAR(grid.probability(longest, false, true), alpha, 1e-9);
	EXPECT_NEAR(grid.probability(longest, true, true), alpha, 1e-9);
}

// --------------------------------------------------------------------------
// The grid against an independent route to the same numbers
// --------------------------------------------------------------------------

TEST(TransitionGrid, the_recursion_matches_a_numerical_matrix_exponential_at_every_bin_centre) {
	// The grid reaches bin k by multiplying its way up from bin 0, never by evaluating anything at
	// bin k's own branch length -- see docs/adr/0003. The two routes must nonetheless agree, and
	// this is what says so.
	const TBinGrid bins(N_BINS);
	for (const double alpha : ALPHAS) {
		for (const double nu : NUS) {
			const TTransitionGrid grid(alpha, nu, bins);
			for (size_t bin = 0; bin < N_BINS; ++bin) {
				const double t = bins.grid_branch_length(bin);
				for (const bool from : {false, true}) {
					for (const bool to : {false, true}) {
						EXPECT_NEAR(grid.probability(bin, from, to),
						            numerical_expm(alpha, nu, t, from, to), 1e-12)
						    << "alpha " << alpha << " nu " << nu << " bin " << bin << " from "
						    << from << " to " << to;
					}
				}
			}
		}
	}
}

// --------------------------------------------------------------------------
// Neutrality (ADR-0001)
// --------------------------------------------------------------------------

TEST(TransitionGrid, neutral_parameters_give_exactly_uninformative_rows) {
	// The assumption the whole independent-field validation rests on: past the short-circuit, with
	// alpha 0.5, every row is exactly (0.5, 0.5) -- not approximately.
	const TBinGrid bins(N_BINS);
	const double neutral_nu = std::exp(5.0);
	const TTransitionGrid grid(0.5, neutral_nu, bins);

	for (size_t bin = 0; bin < N_BINS; ++bin) {
		for (const bool from : {false, true}) {
			for (const bool to : {false, true}) {
				EXPECT_DOUBLE_EQ(grid.probability(bin, from, to), 0.5)
				    << "bin " << bin << " from " << from << " to " << to;
			}
		}
	}
}

TEST(TransitionGrid, the_short_circuit_is_stationary_for_any_alpha) {
	const TBinGrid bins(N_BINS);
	const double alpha = 0.72;
	const TTransitionGrid grid(alpha, std::exp(5.0), bins);

	for (size_t bin = 0; bin < N_BINS; ++bin) {
		EXPECT_DOUBLE_EQ(grid.probability(bin, false, true), alpha);
		EXPECT_DOUBLE_EQ(grid.probability(bin, true, true), alpha);
		EXPECT_DOUBLE_EQ(grid.probability(bin, false, false), 1.0 - alpha);
		EXPECT_DOUBLE_EQ(grid.probability(bin, true, false), 1.0 - alpha);
	}
}

TEST(TransitionGrid, the_short_circuit_is_seamless_at_the_mean_branch_but_not_at_a_short_one) {
	// The nu > 25 threshold was derived at t = 1.0, the mean grid branch length, but it is applied
	// at every bin -- and bin 0 stands for t = Delta/2, eleven times shorter. So crossing the
	// threshold is continuous at bin 5 and a visible jump at bin 0.
	//
	// This test pins that asymmetry rather than asserting the threshold is uniformly accurate,
	// which it is not. Neutrality (ADR-0001) is unaffected: the short-circuit writes exactly
	// (0.5, 0.5), so a pinned dimension is uninformative by construction, not by approximation.
	const TBinGrid bins(N_BINS);
	const double alpha = 0.4;
	const TTransitionGrid just_below(alpha, 25.0, bins);
	const TTransitionGrid just_above(alpha, 25.0 + 1e-9, bins);

	const size_t mean_branch = N_BINS / 2;
	ASSERT_NEAR(bins.grid_branch_length(mean_branch), 1.0, 0.1);
	for (const bool from : {false, true}) {
		EXPECT_NEAR(just_below.probability(mean_branch, from, true),
		            just_above.probability(mean_branch, from, true), 1e-10)
		    << "the threshold should be seamless where it was derived";
	}

	// Bin 0 is where the approximation is worst; ~6e-2 today. Asserted loosely, as a tripwire: if
	// this ever shrinks, the threshold has been made bin-aware and this comment is stale.
	const double jump_at_shortest =
	    std::abs(just_below.probability(0, false, true) - just_above.probability(0, false, true));
	EXPECT_GT(jump_at_shortest, 1e-3)
	    << "the short-circuit no longer jumps at the shortest branch -- has the threshold changed?";
}

// --------------------------------------------------------------------------
// The stationary distribution the roots are drawn from
// --------------------------------------------------------------------------

TEST(TransitionGrid, stationary_is_the_two_state_distribution) {
	const TBinGrid bins(N_BINS);
	for (const double alpha : ALPHAS) {
		const TTransitionGrid grid(alpha, 0.6, bins);
		EXPECT_DOUBLE_EQ(grid.stationary(true), alpha);
		EXPECT_DOUBLE_EQ(grid.stationary(false), 1.0 - alpha);
		EXPECT_DOUBLE_EQ(grid.stationary(true) + grid.stationary(false), 1.0);
		EXPECT_DOUBLE_EQ(grid.alpha(), alpha);
	}
}
