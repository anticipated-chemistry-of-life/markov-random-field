//
// The two-state draw: state 1 with the probability its caller names.
//
// A cell update ends in one question: is this cell a one? Callers name the probability in three
// ways -- as a probability, as a pair of log probabilities, or as the pair of sums a walk
// accumulates -- so three overloads share one body. The node-state walk asks the third.
//
// Every caller supplies the uniform, and it comes from the cell being drawn
// (random/TCellUniforms.h) rather than from a running generator. A caller that cannot name its cell
// therefore cannot draw. ADR-0007 says why.
//

#ifndef ACOL_TWO_STATE_DRAW_H
#define ACOL_TWO_STATE_DRAW_H

#include "coretools/Math/TSumLog.h"
#include "coretools/Types/probability.h"

#include <array>
#include <cmath>

namespace two_state_draw {

namespace detail {
/// `1 / (1 + exp(-log_odds))` is the probability of state 1. Log odds of NaN make the comparison
/// false, so a pair of probabilities that says nothing reads as state 0.
///
/// The tails need no branch of their own: the exponential saturates at either end and gives the
/// answer that tail asks for. ADR-0007 says why this is written here rather than taken from the
/// odds-ratio helper it replaces.
///
/// The comparison is written out rather than handed to the probability overload, because a weak
/// probability type rejects the NaN this has to let through.
inline bool state_one(double log_odds, double uniform) {
	return uniform < 1.0 / (1.0 + std::exp(-log_odds));
}
} // namespace detail

/// State 1 at the rate `probability_of_one` names.
inline bool sample(coretools::Probability probability_of_one, double uniform) {
	return uniform < probability_of_one;
}

/// State 1 from the two sums a walk accumulates, one per state.
inline bool sample(std::array<coretools::TSumLogProbability, 2> &sum_log, double uniform) {
	return detail::state_one(sum_log[1].getSum() - sum_log[0].getSum(), uniform);
}

/// State 1 from a pair of log probabilities, which need not be normalised.
inline bool sample(double log_prob_0, double log_prob_1, double uniform) {
	return detail::state_one(log_prob_1 - log_prob_0, uniform);
}

} // namespace two_state_draw

#endif // ACOL_TWO_STATE_DRAW_H
