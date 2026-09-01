//
// The likelihood of the "simple error model" data source.
//
// The simple error model treats the observed matrix D as a direct, noisy observation of the latent
// field Y: every cell is reported correctly with probability 1 - epsilon and inverted with
// probability epsilon. D therefore has exactly the same dimensions as Y, which makes it a
// deliberately uninformative-about-structure data source -- useful for isolating
// whether a convergence problem lives in a data likelihood or in the Markov field itself.
//
// Because every cell contributes either log(1 - epsilon) or log(epsilon), the total log-likelihood
// depends on the data only through the number of cells where D and Y disagree. Keeping that single
// integer turns the epsilon_simple_model Metropolis-Hastings move into an O(1) computation instead
// of a sum over the whole container space.
//
// Everything here is free-standing: it is deliberately NOT guarded by USE_SIMPLE_ERROR_MODEL, so
// it compiles (and is unit-tested) in every build configuration.
//

#pragma once

#include "coretools/Main/TError.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Types/probability.h"
#include "storages/storage_backend.h"
#include <array>
#include <cmath>
#include <cstddef>

namespace simple_error_model {

/// P(D = d | Y = y): 1 - eps when the observation agrees with the latent state, eps when it does
/// not.
[[nodiscard]] constexpr double probability_of_D_given_Y(bool y, bool d, double eps) noexcept {
	return (y == d) ? 1.0 - eps : eps;
}

/// Fills prob[0] = P(D = d | Y = 0) and prob[1] = P(D = d | Y = 1) -- the two values the Y Gibbs
/// step needs for a single cell. Writing them by index rather than branching keeps this free of
/// conditionals on the hot path.
constexpr void probabilities_for_both_Y_states(bool d, double eps,
                                               std::array<double, 2> &prob) noexcept {
	prob[static_cast<size_t>(d)]  = 1.0 - eps; // Y == d -> agreement
	prob[static_cast<size_t>(!d)] = eps;       // Y != d -> disagreement
}

/// Total log-likelihood of the simple error model:
///     (total - n_disagree) * log(1 - eps) + n_disagree * log(eps)
/// Every cell contributes one of only two values, so the whole data set collapses into these two
/// counts. This is what makes the epsilon_simple_model MCMC move O(1) in the number of cells.
[[nodiscard]] inline double log_likelihood_from_counts(size_t total, size_t n_disagree,
                                                       double eps) {
	if (n_disagree > total) {
		throw coretools::TDevError("Number of disagreeing cells (", n_disagree,
		                           ") exceeds the total number of cells (", total, ").");
	}
	const size_t n_agree = total - n_disagree;
	// log1p(-eps) instead of log(1 - eps): eps is typically small, and 1 - eps loses precision
	// exactly there.
	return static_cast<double>(n_agree) * std::log1p(-eps) +
	       static_cast<double>(n_disagree) * std::log(eps);
}

/// Log-likelihood ratio of a proposed epsilon against the current one. The disagreement count is a
/// property of the data and of Y, not of epsilon, so it cancels out of neither term but stays
/// fixed -- no pass over the cells is needed.
[[nodiscard]] inline double log_likelihood_ratio(size_t total, size_t n_disagree, double old_eps,
                                                 double new_eps) {
	return log_likelihood_from_counts(total, n_disagree, new_eps) -
	       log_likelihood_from_counts(total, n_disagree, old_eps);
}

/// Number of cells where Y and D differ.
///
/// A cell that is zero in both always agrees and need not be visited, so the two fields are
/// merge-joined over their ones in ascending linear-index order (the same technique as
/// TLotus::calculate_log_likelihood_of_L), which costs O(ones(Y) + ones(D)) rather
/// than O(total cells).
[[nodiscard]] inline size_t count_disagreements(const TFieldStorage &Y, const TFieldStorage &D) {
	if (Y.dimensions() != D.dimensions()) {
		throw coretools::TDevError(
		    "Cannot compare Y and the simple error model data D: they have different dimensions (Y "
		    "is ",
		    Y.dimensions()[0], "x", Y.dimensions()[1], ", D is ", D.dimensions()[0], "x",
		    D.dimensions()[1], ").");
	}

	const size_t total = Y.total_size_of_container_space();
	auto y_cur         = Y.ones_cursor();
	auto d_cur         = D.ones_cursor();

	size_t n_disagree = 0;
	while (y_cur.valid() || d_cur.valid()) {
		const size_t yi = y_cur.valid() ? y_cur.linear_index() : total;
		const size_t di = d_cur.valid() ? d_cur.linear_index() : total;
		const size_t i  = std::min(yi, di);

		bool state_of_Y = false;
		bool state_of_D = false;
		if (yi == i) {
			state_of_Y = true; // the cursor yields only ones
			y_cur.advance();
		}
		if (di == i) {
			state_of_D = true;
			d_cur.advance();
		}
		if (state_of_Y != state_of_D) { ++n_disagree; }
	}
	return n_disagree;
}

/// Draw one observed cell from a latent state: report y with probability 1 - eps, its inverse with
/// probability eps. Used when simulating D from a simulated Y.
[[nodiscard]] inline bool draw_D_given_Y(bool y, double eps) {
	const coretools::Probability p_flip(eps);
	const bool flip = coretools::instances::randomGenerator().pickOneOfTwo(p_flip);
	return flip ? !y : y;
}

} // namespace simple_error_model
