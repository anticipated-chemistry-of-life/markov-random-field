//
// The link between the two tree fields and the field, and the eight-state block the sampler draws
// from.
//
// Each tree carries its own leaf-level view of the field -- its tree field -- and the field is a
// noisy reconciliation of the two: corrupt each tree field cell independently with probability
// omega, then take the AND. That gives the four link probabilities
//
//     P(Y = 1 | Z_s = 1, Z_m = 1) = (1 - omega)^2
//     P(Y = 1 | Z_s = 1, Z_m = 0) = (1 - omega) * omega
//     P(Y = 1 | Z_s = 0, Z_m = 1) = omega * (1 - omega)
//     P(Y = 1 | Z_s = 0, Z_m = 0) = omega^2
//
// and, because the joint is then a directed factorisation rather than a product of two likelihoods
// over one shared variable, a normalising constant that is identically 1. That is the whole point
// of the model; see ADR-0005, which carries the derivations these functions implement.
//
// The table depends on the two tree field states only through their sum, so the link's whole
// contribution to the likelihood collapses to six integers -- n(bucket, field state) -- and the
// error-probability Metropolis move becomes O(1) in the number of cells rather than a sum over
// every cell. That
// is the same trick the simple error model's disagreement count plays for its own rate.
//
// Everything here is free-standing. It takes two transition probabilities, an error probability
// and two data-likelihood pairs, and it touches no storage, no tree, no parameter and no random
// generator -- which is what makes the closed forms below testable against brute-force definitions
// for a few lines of arithmetic, instead of through a constructed field.
//
// It is also deliberately NOT guarded by the data-source build flags, so it compiles and is unit
// tested in every build configuration, following lotus/TLotusMath.h and
// simple_error_model/TSimpleErrorModelMath.h.
//

#pragma once

#include "constants.h"
#include "coretools/Types/probability.h"
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace field_math {

/// The probability that a tree field cell is corrupted before the link reads it.
///
/// A struct rather than a bare double, so that giving each tree its own error probability later is
/// a change inside this type and the policy that reads it, rather than a signature change
/// everywhere. `is_shared` is what the six-counter collapse depends on: bucketing by the number of
/// tree fields in state 1 pools the two mixed cells, which is only valid while both trees are
/// corrupted at the same rate.
///
/// Constrained to the open interval (0, 0.5). At 0 the link is the deterministic AND and the block
/// update hits log 0; at or above 0.5 the tree fields are anti-correlated with the field, which is
/// meaningless as an error model and is a genuine second mode for a sampler to find (ADR-0005).
///
/// The check below is written out rather than delegated to one of coretools' interval types, for
/// two reasons. No such type spells (0, 0.5) -- `ZeroOneOpen` would leave the upper half of the
/// range unguarded. And coretools checks its intervals only under CHECK_INTERVALS, which
/// CMakeLists.txt defines for the unit tests and not for acol, so an interval type would stop
/// checking in exactly the binary that runs the chain. This is a statement about the model, not a
/// range on an argument, so it holds in every build.
class TErrorProbability {
private:
	double _omega = 0.0;

public:
	// No default constructor: it would have to pick a value, and 0.0 -- the natural choice -- is
	// exactly the one the constructor below rejects.
	explicit TErrorProbability(double omega) : _omega(omega) {
		if (!(omega > 0.0) || !(omega < 0.5)) {
			throw std::invalid_argument(
			    "The error probability must lie strictly inside (0, 0.5), but it is " +
			    std::to_string(omega) +
			    ". At 0 the link is deterministic and the block update takes log(0); at 0.5 and "
			    "above the tree fields are anti-correlated with the field.");
		}
	}

	/// The rate at which `tree`'s field is corrupted. One value today, one per tree later --
	/// which is the whole reason this is a type and not a double.
	[[nodiscard]] double for_tree(size_t tree) const {
		if (tree >= NUMBER_OF_TREES) {
			throw std::invalid_argument("There is no tree " + std::to_string(tree) +
			                            "; there are " + std::to_string(NUMBER_OF_TREES) + ".");
		}
		return _omega;
	}

	/// Whether both trees are corrupted at the same rate, which is what lets the link likelihood
	/// collapse to three buckets instead of four cells.
	[[nodiscard]] static constexpr bool is_shared() noexcept { return true; }
};

/// Bounds check shared by everything that takes a bucket index.
inline void check_bucket(size_t bucket, size_t n_buckets) {
	if (bucket >= n_buckets) {
		throw std::invalid_argument("There is no bucket " + std::to_string(bucket) +
		                            "; there are " + std::to_string(n_buckets) + ".");
	}
}

/// Index of one of the eight states a leaf pair can be in, over (field, species tree field,
/// molecule tree field). The field takes the high bit so that the four states sharing a field
/// value are contiguous.
[[nodiscard]] constexpr size_t state_index(bool y, bool z_s, bool z_m) noexcept {
	return (static_cast<size_t>(y) << 2U) | (static_cast<size_t>(z_s) << 1U) |
	       static_cast<size_t>(z_m);
}

/// The sufficient statistic of the link: `n(bucket, field state)`.
///
/// Six integers carry the whole link likelihood, so the error-probability move never walks the
/// cells. The block update knows its own old and new (bucket, field state), which is why this is
/// maintained by `add` and `remove` rather than recomputed.
///
/// The bucket count is the AND's; a link that bucketed differently would change this type with it.
class TLinkCounters {
public:
	static constexpr size_t n_buckets = 3;

private:
	std::array<std::array<size_t, 2>, n_buckets> _n{};

public:
	void add(size_t bucket, bool y) {
		check_bucket(bucket, n_buckets);
		++_n[bucket][static_cast<size_t>(y)];
	}

	void remove(size_t bucket, bool y) {
		check_bucket(bucket, n_buckets);
		size_t &count = _n[bucket][static_cast<size_t>(y)];
		if (count == 0) {
			throw std::invalid_argument("Cannot remove a cell from bucket " +
			                            std::to_string(bucket) + " with field state " +
			                            std::to_string(static_cast<int>(y)) + ": it is empty.");
		}
		--count;
	}

	[[nodiscard]] size_t count(size_t bucket, bool y) const {
		check_bucket(bucket, n_buckets);
		return _n[bucket][static_cast<size_t>(y)];
	}

	[[nodiscard]] size_t total() const noexcept {
		size_t sum = 0;
		for (const auto &bucket : _n) { sum += bucket[0] + bucket[1]; }
		return sum;
	}
};

/// What the sampler needs from a link: the probability of a field cell given the two tree field
/// cells, the bucket that cell falls in, and the likelihood of a whole configuration from the
/// bucket counts alone.
template<typename T>
concept LinkPolicy = requires(bool z_s, bool z_m, size_t bucket, const TErrorProbability &omega,
                              const TLinkCounters &counters) {
	{ T::n_buckets } -> std::convertible_to<size_t>;
	// the counters are sized for the link, so a link that bucketed differently would have to
	// change TLinkCounters with it rather than silently overflow it
	requires T::n_buckets == TLinkCounters::n_buckets;
	{ T::prob_y_is_one(z_s, z_m, omega) } -> std::same_as<double>;
	{ T::bucket(z_s, z_m) } -> std::same_as<size_t>;
	{ T::prob_for_bucket(bucket, omega) } -> std::same_as<double>;
	{ T::log_likelihood(counters, omega) } -> std::same_as<double>;
};

/// The field is the AND of the two independently corrupted tree fields.
///
/// The normalisation win comes from the *structure* -- each tree owning its own leaf-level field --
/// and not from this table, so a different link can be dropped in here without revisiting the
/// argument that the joint is normalised (ADR-0005).
class TAndLinkPolicy {
public:
	static constexpr size_t n_buckets = TLinkCounters::n_buckets;

	/// P(Y = 1 | Z_s = z_s, Z_m = z_m). Written as a product of the two corrupted reads, so it is
	/// already correct for a per-tree error probability.
	[[nodiscard]] static double prob_y_is_one(bool z_s, bool z_m, const TErrorProbability &omega) {
		const double omega_s = omega.for_tree(0);
		const double omega_m = omega.for_tree(1);
		return (z_s ? 1.0 - omega_s : omega_s) * (z_m ? 1.0 - omega_m : omega_m);
	}

	/// The number of tree fields in state 1, which is all the link depends on.
	[[nodiscard]] static size_t bucket(bool z_s, bool z_m) noexcept {
		return static_cast<size_t>(z_s) + static_cast<size_t>(z_m);
	}

	/// `P_k = (1 - omega)^k * omega^(2 - k)`.
	///
	/// Two parameter-free constraints follow, and they test different assumptions:
	/// `P_1^2 = P_0 * P_2` is the six-counter shadow of the four-cell identity and tests the
	/// independent-corruption-AND structure, while `sqrt(P_0) + sqrt(P_2) = 1` tests that the two
	/// trees share one error probability. Only the second has no blind spot -- pooling the two
	/// mixed cells means there is always a split at which the first holds despite unequal rates.
	/// See ADR-0005.
	[[nodiscard]] static double prob_for_bucket(size_t bucket, const TErrorProbability &omega) {
		check_bucket(bucket, n_buckets);
		if (!field_math::TErrorProbability::is_shared()) {
			throw std::invalid_argument(
			    "Bucketing by the number of tree fields in state 1 pools the two mixed cells, "
			    "which "
			    "is only the same probability while both trees are corrupted at one rate.");
		}
		const double w = omega.for_tree(0);
		switch (bucket) {
		case 0: return w * w;
		case 1: return w * (1.0 - w);
		default: return (1.0 - w) * (1.0 - w);
		}
	}

	/// `{ log(1 - P_k), log P_k }`, indexed by the field state, each from a closed form chosen so
	/// that nothing cancels.
	///
	/// Taking logs of `prob_for_bucket` would be shorter and is what this did first. It is worse in
	/// both directions, and not where one would guess. `1 - P_k` is *exact* whenever `P_k >= 0.5`
	/// (Sterbenz), so `log1p` buys nothing for the both-ones bucket; the loss is at the other end,
	/// where `P_0 = omega^2` is tiny and `1 - P_0` rounds to 1. And `log P_2` formed as
	/// `log((1 - w) * (1 - w))` loses the squaring's rounding -- 4e-11 relative at `omega = 1e-6`
	/// -- where ADR-0005's observation that `log P_k` is *affine* in `k` gives it exactly.
	///
	/// So: the logs come from the affine form, and the complements from expansions that subtract
	/// nothing near-equal -- `1 - P_0 = (1-w)(1+w)`, `1 - P_1 = 1 - w + w^2`, `1 - P_2 = w(2-w)`.
	[[nodiscard]] static std::array<double, 2> log_prob_for_bucket(size_t bucket,
	                                                               const TErrorProbability &omega) {
		check_bucket(bucket, n_buckets);
		const double w         = omega.for_tree(0);
		const double log_w     = std::log(w);
		const double log_1_m_w = std::log1p(-w);

		switch (bucket) {
		case 0: return {log_1_m_w + std::log1p(w), 2.0 * log_w};        // 1 - P_0 = (1-w)(1+w)
		case 1: return {std::log1p(-w + w * w), log_w + log_1_m_w};     // 1 - P_1 = 1 - w + w^2
		default: return {log_w + std::log1p(1.0 - w), 2.0 * log_1_m_w}; // 1 - P_2 = w(2-w)
		}
	}

	/// `sum over buckets of  n(k,1) * log P_k  +  n(k,0) * log(1 - P_k)`.
	///
	/// The counts are a property of the configuration, not of omega, so a Metropolis move on the
	/// error probability recomputes this from six integers rather than from the cells.
	[[nodiscard]] static double log_likelihood(const TLinkCounters &counters,
	                                           const TErrorProbability &omega) {
		double sum = 0.0;
		for (size_t bucket = 0; bucket < n_buckets; ++bucket) {
			const auto log_p = log_prob_for_bucket(bucket, omega);
			// omega is strictly inside (0, 0.5), so neither log is infinite and a zero count
			// contributes a clean zero.
			sum += static_cast<double>(counters.count(bucket, false)) * log_p[0];
			sum += static_cast<double>(counters.count(bucket, true)) * log_p[1];
		}
		return sum;
	}
};

/// The eight-state block at one leaf pair: the field and the two tree fields, updated together.
///
/// Their combined Markov blanket is these six scalars and the error probability. Sampling the
/// three as one block rather than as three single-site steps is what escapes the metastable state
/// the AND creates -- with a small omega, a field cell at one pins both tree fields, and given both
/// at one the field stays at one with probability near one. That triple mixes arbitrarily slowly
/// under single-site Gibbs, and it would present as a convergence problem rather than as a bug.
///
/// Returns the normalised probability of each of the eight states, addressed by `state_index`.
/// Drawing from them is the caller's business; nothing here touches a random generator.
///
/// @param prob_z_s_is_one  P(Z_s = 1 | the species parent's state), from that tree's transition
/// grid.
/// @param prob_z_m_is_one  P(Z_m = 1 | the molecule parent's state).
/// @param omega            The error probability standing between the tree fields and the field.
/// @param lotus            {P(L | Y = 0), P(L | Y = 1)} for this cell.
/// @param simple_error     {P(D | Y = 0), P(D | Y = 1)} for this cell.
template<LinkPolicy Policy>
[[nodiscard]] std::array<double, 8>
block_probabilities(coretools::Probability prob_z_s_is_one, coretools::Probability prob_z_m_is_one,
                    const TErrorProbability &omega,
                    const std::array<coretools::Probability, 2> &lotus,
                    const std::array<coretools::Probability, 2> &simple_error) {
	// These are probabilities by type, so the range check happens where the caller builds the
	// value, the way TSimpleErrorModelMath.h already does it. Note the bargain: coretools checks
	// its intervals only under CHECK_INTERVALS, which CMakeLists.txt defines for the unit tests and
	// not for acol. So this is a development-time guard, and a transition probability outside [0,1]
	// reaching a release build would still produce a negative weight. Every coretools type in this
	// codebase makes the same trade.
	const double p_s = prob_z_s_is_one.get();
	const double p_m = prob_z_m_is_one.get();

	std::array<double, 8> probability{};
	double total = 0.0;

	for (const bool z_s : {false, true}) {
		for (const bool z_m : {false, true}) {
			const double tree_factor = (z_s ? p_s : 1.0 - p_s) * (z_m ? p_m : 1.0 - p_m);
			const double link        = Policy::prob_y_is_one(z_s, z_m, omega);
			for (const bool y : {false, true}) {
				const auto data     = static_cast<size_t>(y);
				const double weight = tree_factor * (y ? link : 1.0 - link) * lotus[data].get() *
				                      simple_error[data].get();
				probability[state_index(y, z_s, z_m)] = weight;
				total += weight;
			}
		}
	}

	if (!(total > 0.0)) {
		throw std::invalid_argument(
		    "The eight-state block has no probability mass: the eight weights sum to " +
		    std::to_string(total) +
		    ". Either a data likelihood is zero for both field states, or the weights "
		    "underflowed.");
	}
	for (double &p : probability) { p /= total; }
	return probability;
}

} // namespace field_math
