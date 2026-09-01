#include "coretools/Types/probability.h"
#include "field/TFieldMath.h"
#include "field/link_backend.h"
#include "gtest/gtest.h"
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using namespace field_math;

/// The definition the link table is a shortcut for: corrupt each tree field cell independently
/// with probability omega, then AND the two corrupted values.
///
/// This enumerates both corruption events, so it shares no algebra with `prob_y_is_one` -- the
/// same closed-form-versus-independent-brute-force pattern the simple error model's test uses.
double brute_force_prob_y_is_one(bool z_s, bool z_m, double omega) {
	double p = 0.0;
	for (const bool flip_s : {false, true}) {
		for (const bool flip_m : {false, true}) {
			const double weight = (flip_s ? omega : 1.0 - omega) * (flip_m ? omega : 1.0 - omega);
			// a corrupted cell is `z != flip`; the link is the AND of the two
			if ((z_s != flip_s) && (z_m != flip_m)) { p += weight; }
		}
	}
	return p;
}

/// One cell of a configuration: the two tree field states and the field state.
struct TCell {
	bool z_s = false;
	bool z_m = false;
	bool y   = false;
};

/// The log-likelihood written the long way: walk every cell and add its own term. Shares no code
/// with the six-counter closed form.
double brute_force_log_likelihood(const std::vector<TCell> &cells, double omega) {
	double sum = 0.0;
	for (const auto &cell : cells) {
		const double p = brute_force_prob_y_is_one(cell.z_s, cell.z_m, omega);
		sum += cell.y ? std::log(p) : std::log(1.0 - p);
	}
	return sum;
}

/// The eight-state block built from its definition -- the product of the two tree factors, the
/// link, and the data, normalised -- rather than from `block_probabilities`.
///
/// Being honest about what this does and does not check: a posterior over eight states is
/// *necessarily* prior x link x data, normalised, so that shape is shared and no rewriting makes it
/// otherwise. What is independent is the link, which comes from the corruption enumeration above
/// rather than from any closed form. The relations that genuinely do not share arithmetic with
/// `block_probabilities` are the three below it -- the tree marginal, the data-tilted odds, and the
/// corrupted-rate identity.
std::array<double, 8> brute_force_block(double p_s, double p_m, double omega,
                                        const std::array<double, 2> &lotus,
                                        const std::array<double, 2> &simple_error) {
	std::array<double, 8> weight{};
	double total = 0.0;
	for (const bool y : {false, true}) {
		for (const bool z_s : {false, true}) {
			for (const bool z_m : {false, true}) {
				const double link = brute_force_prob_y_is_one(z_s, z_m, omega);
				const double w    = (z_s ? p_s : 1.0 - p_s) * (z_m ? p_m : 1.0 - p_m) *
				                 (y ? link : 1.0 - link) * lotus[static_cast<size_t>(y)] *
				                 simple_error[static_cast<size_t>(y)];
				weight[state_index(y, z_s, z_m)] = w;
				total += w;
			}
		}
	}
	for (auto &w : weight) { w /= total; }
	return weight;
}

/// The block takes probabilities by type. The brute force below stays on doubles, so that it
/// shares as little as possible with the thing it checks.
std::array<coretools::Probability, 2> as_probabilities(const std::array<double, 2> &pair) {
	return {coretools::P(pair[0]), coretools::P(pair[1])};
}

/// The range of values across the open interval the error probability is constrained to.
const std::vector<double> &omega_values() {
	static const std::vector<double> values = {1e-6, 0.001, 0.01, 0.05, 0.1, 0.2,
	                                          0.25, 0.3,   0.4,  0.45, 0.499};
	return values;
}

/// The link probability of a bucket as the counters see it: n(k,1) / (n(k,0) + n(k,1)).
/// No parameter is estimated first, which is what makes the two constraints falsifiable.
double empirical_prob(const TLinkCounters &counters, size_t bucket) {
	const double ones  = static_cast<double>(counters.count(bucket, true));
	const double zeros = static_cast<double>(counters.count(bucket, false));
	return ones / (ones + zeros);
}

/// Counters realising a chosen number of ones and zeros per bucket.
TLinkCounters counters_from(const std::array<std::pair<size_t, size_t>, 3> &ones_and_zeros) {
	TLinkCounters counters;
	for (size_t bucket = 0; bucket < ones_and_zeros.size(); ++bucket) {
		for (size_t i = 0; i < ones_and_zeros[bucket].first; ++i) { counters.add(bucket, true); }
		for (size_t i = 0; i < ones_and_zeros[bucket].second; ++i) { counters.add(bucket, false); }
	}
	return counters;
}

} // namespace

//-----------------------------------
// TErrorProbability
//-----------------------------------

TEST(FieldMath_Tests, error_probability_reports_the_same_value_for_both_trees) {
	const TErrorProbability omega(0.125);
	EXPECT_DOUBLE_EQ(omega.for_tree(0), 0.125);
	EXPECT_DOUBLE_EQ(omega.for_tree(1), 0.125);
	EXPECT_TRUE(omega.is_shared());
}

TEST(FieldMath_Tests, error_probability_rejects_values_outside_the_open_interval) {
	EXPECT_THROW(TErrorProbability(0.0), std::invalid_argument);
	EXPECT_THROW(TErrorProbability(0.5), std::invalid_argument);
	EXPECT_THROW(TErrorProbability(-0.1), std::invalid_argument);
	EXPECT_THROW(TErrorProbability(0.7), std::invalid_argument);
	EXPECT_NO_THROW(TErrorProbability(1e-9));
	EXPECT_NO_THROW(TErrorProbability(0.4999));
}

TEST(FieldMath_Tests, error_probability_rejects_a_tree_it_does_not_have) {
	const TErrorProbability omega(0.1);
	EXPECT_THROW(omega.for_tree(NUMBER_OF_TREES), std::invalid_argument);
}

//-----------------------------------
// The link table (ADR-0005, derivation 1)
//-----------------------------------

TEST(FieldMath_Tests, link_table_is_independent_corruption_followed_by_an_and) {
	for (const double w : omega_values()) {
		const TErrorProbability omega(w);
		for (const bool z_s : {false, true}) {
			for (const bool z_m : {false, true}) {
				// exactly one corruption event reaches each cell, so the enumeration reduces to
				// the same single product the closed form writes -- these agree bit for bit
				EXPECT_DOUBLE_EQ(TLinkPolicy::prob_y_is_one(z_s, z_m, omega),
				                 brute_force_prob_y_is_one(z_s, z_m, w))
				    << "omega = " << w << ", Z_s = " << z_s << ", Z_m = " << z_m;
			}
		}
	}
}

TEST(FieldMath_Tests, link_table_matches_the_four_probabilities_in_the_record) {
	const double w = 0.2;
	const TErrorProbability omega(w);
	EXPECT_NEAR(TLinkPolicy::prob_y_is_one(true, true, omega), (1.0 - w) * (1.0 - w), 1e-15);
	EXPECT_NEAR(TLinkPolicy::prob_y_is_one(true, false, omega), (1.0 - w) * w, 1e-15);
	EXPECT_NEAR(TLinkPolicy::prob_y_is_one(false, true, omega), w * (1.0 - w), 1e-15);
	EXPECT_NEAR(TLinkPolicy::prob_y_is_one(false, false, omega), w * w, 1e-15);
}

TEST(FieldMath_Tests, the_link_is_symmetric_in_the_two_trees) {
	for (const double w : omega_values()) {
		const TErrorProbability omega(w);
		EXPECT_DOUBLE_EQ(TLinkPolicy::prob_y_is_one(true, false, omega),
		                 TLinkPolicy::prob_y_is_one(false, true, omega))
		    << "omega = " << w;
	}
}

//-----------------------------------
// The bucketing (ADR-0005, derivation 2)
//-----------------------------------

TEST(FieldMath_Tests, bucket_is_the_number_of_tree_fields_in_state_one) {
	EXPECT_EQ(TLinkPolicy::bucket(false, false), 0u);
	EXPECT_EQ(TLinkPolicy::bucket(true, false), 1u);
	EXPECT_EQ(TLinkPolicy::bucket(false, true), 1u);
	EXPECT_EQ(TLinkPolicy::bucket(true, true), 2u);
}

TEST(FieldMath_Tests, cells_in_one_bucket_share_a_link_probability) {
	for (const double w : omega_values()) {
		const TErrorProbability omega(w);
		for (const bool z_s : {false, true}) {
			for (const bool z_m : {false, true}) {
				EXPECT_DOUBLE_EQ(TLinkPolicy::prob_for_bucket(TLinkPolicy::bucket(z_s, z_m), omega),
				                 TLinkPolicy::prob_y_is_one(z_s, z_m, omega))
				    << "omega = " << w << ", Z_s = " << z_s << ", Z_m = " << z_m;
			}
		}
	}
}

//-----------------------------------
// The six counters (ADR-0005, derivation 2)
//-----------------------------------

TEST(FieldMath_Tests, counters_hold_one_entry_per_bucket_and_field_state) {
	TLinkCounters counters;
	EXPECT_EQ(counters.total(), 0u);
	counters.add(TLinkPolicy::bucket(true, true), true);
	counters.add(TLinkPolicy::bucket(true, true), true);
	counters.add(TLinkPolicy::bucket(false, false), false);
	EXPECT_EQ(counters.count(2, true), 2u);
	EXPECT_EQ(counters.count(0, false), 1u);
	EXPECT_EQ(counters.count(1, true), 0u);
	EXPECT_EQ(counters.total(), 3u);
}

TEST(FieldMath_Tests, counters_are_maintained_incrementally) {
	// what the block update does: it knows its own old and new (bucket, field state)
	TLinkCounters counters;
	counters.add(1, true);
	counters.add(1, true);
	counters.remove(1, true);
	counters.add(2, false);
	EXPECT_EQ(counters.count(1, true), 1u);
	EXPECT_EQ(counters.count(2, false), 1u);
	EXPECT_EQ(counters.total(), 2u);
}

TEST(FieldMath_Tests, removing_from_an_empty_bucket_is_rejected) {
	TLinkCounters counters;
	EXPECT_THROW(counters.remove(0, true), std::invalid_argument);
}

TEST(FieldMath_Tests, six_counters_match_a_naive_per_cell_recomputation) {
	// a configuration with every bucket and both field states represented
	std::vector<TCell> cells;
	uint32_t seed = 1u; // fixed width, so the configuration does not depend on the platform
	for (size_t i = 0; i < 400; ++i) {
		seed          = seed * 1103515245u + 12345u;
		const bool zs = ((seed >> 16u) & 1u) != 0u;
		const bool zm = ((seed >> 17u) & 1u) != 0u;
		const bool y  = ((seed >> 18u) & 1u) != 0u;
		cells.push_back({zs, zm, y});
	}

	TLinkCounters counters;
	for (const auto &cell : cells) { counters.add(TLinkPolicy::bucket(cell.z_s, cell.z_m), cell.y); }
	ASSERT_EQ(counters.total(), cells.size());
	// the comparison is only worth anything if every one of the six counters is live
	for (size_t bucket = 0; bucket < TLinkCounters::n_buckets; ++bucket) {
		for (const bool y : {false, true}) {
			ASSERT_GT(counters.count(bucket, y), 0u) << "bucket " << bucket << ", field state " << y;
		}
	}

	for (const double w : omega_values()) {
		const TErrorProbability omega(w);
		const double want = brute_force_log_likelihood(cells, w);
		// Relative, not absolute. The residual here is the brute force's, not the closed form's:
		// walking cell by cell it forms P_k by squaring and then takes a log, which costs 4e-11
		// relative at omega = 1e-6, where log_prob_for_bucket uses the affine form and does not.
		EXPECT_NEAR(TLinkPolicy::log_likelihood(counters, omega), want, 1e-11 * std::abs(want))
		    << "omega = " << w;
	}
}

TEST(FieldMath_Tests, an_empty_configuration_has_zero_log_likelihood) {
	const TLinkCounters counters;
	EXPECT_DOUBLE_EQ(TLinkPolicy::log_likelihood(counters, TErrorProbability(0.1)), 0.0);
}

//-----------------------------------
// The two parameter-free constraints (ADR-0005, derivation 2)
//-----------------------------------

TEST(FieldMath_Tests, the_and_identity_holds_across_the_error_probability_range) {
	// P_1^2 = P_0 * P_2, for every omega. This is what makes the AND falsifiable from the
	// counters alone, with no parameter estimated first.
	for (const double w : omega_values()) {
		const TErrorProbability omega(w);
		const double p_0 = TLinkPolicy::prob_for_bucket(0, omega);
		const double p_1 = TLinkPolicy::prob_for_bucket(1, omega);
		const double p_2 = TLinkPolicy::prob_for_bucket(2, omega);
		// exact in real arithmetic, so the only slack allowed is rounding: EXPECT_DOUBLE_EQ is
		// 4 ULP, and an absolute tolerance here would sit below one ULP of P_0 * P_2
		EXPECT_DOUBLE_EQ(p_1 * p_1, p_0 * p_2) << "omega = " << w;
	}
}

TEST(FieldMath_Tests, the_shared_error_probability_constraint_holds_across_the_range) {
	// sqrt(P_0) + sqrt(P_2) = 1. Unlike the identity above this one has no blind spot: it fails
	// exactly when the two trees do not share one error probability (ADR-0005).
	for (const double w : omega_values()) {
		const TErrorProbability omega(w);
		EXPECT_NEAR(std::sqrt(TLinkPolicy::prob_for_bucket(0, omega)) +
		                std::sqrt(TLinkPolicy::prob_for_bucket(2, omega)),
		            1.0, 1e-15)
		    << "omega = " << w;
	}
}

TEST(FieldMath_Tests, both_constraints_hold_on_counters_from_one_shared_error_probability) {
	// omega = 0.2 gives P_0 = 0.04, P_1 = 0.16, P_2 = 0.64, each realised exactly by these counts.
	// Unlike the two tests above, this reads the constraints off n(k, y) rather than off the
	// closed form they were derived from, so it is a statement about a configuration.
	const auto counters = counters_from({{{400, 9600}, {1600, 8400}, {6400, 3600}}});
	const double p_0    = empirical_prob(counters, 0);
	const double p_1    = empirical_prob(counters, 1);
	const double p_2    = empirical_prob(counters, 2);

	EXPECT_NEAR(p_1 * p_1, p_0 * p_2, 1e-15);
	EXPECT_NEAR(std::sqrt(p_0) + std::sqrt(p_2), 1.0, 1e-15);
}

TEST(FieldMath_Tests, the_identity_has_a_blind_spot_that_the_shared_rate_constraint_does_not) {
	// The claim ADR-0005 makes, as a test. Corrupt the two trees at *different* rates -- omega_s =
	// 0.05, omega_m = 0.2 -- so the four link probabilities are P_11 = 0.76, P_10 = 0.19,
	// P_01 = 0.04, P_00 = 0.01. Bucketing pools the two mixed cells, and at the mix that lands the
	// pooled rate on their geometric mean, sqrt(0.19 * 0.04) = 0.08718, the identity P_1^2 = P_0 P_2
	// holds *exactly* on a model it should reject. Such a mix always exists, because the geometric
	// mean of the two mixed rates always lies between them.
	const auto counters = counters_from({{{1000, 99000}, {8718, 91282}, {76000, 24000}}});
	const double p_0    = empirical_prob(counters, 0);
	const double p_1    = empirical_prob(counters, 1);
	const double p_2    = empirical_prob(counters, 2);

	// the identity is fooled: it sees 3.5e-7 where it would need to see 2.8e-2 to object
	EXPECT_LT(std::abs(p_1 * p_1 - p_0 * p_2), 1e-6);
	// the shared-rate constraint is not, and misses 1 by the amount the record quotes
	EXPECT_NEAR(std::sqrt(p_0) + std::sqrt(p_2), 1.0 - 0.0282202113, 1e-9);
	EXPECT_GT(std::abs(std::sqrt(p_0) + std::sqrt(p_2) - 1.0), 1e-2);
}

//-----------------------------------
// log_prob_for_bucket
//-----------------------------------

TEST(FieldMath_Tests, the_log_probabilities_agree_with_the_logs_of_the_probabilities) {
	// They are not required to be bit-identical -- the cancellation-free forms are the more
	// accurate of the two, by up to 4e-11 relative at the small end -- but they must agree.
	for (const double w : omega_values()) {
		const TErrorProbability omega(w);
		for (size_t bucket = 0; bucket < TLinkCounters::n_buckets; ++bucket) {
			const double p   = TLinkPolicy::prob_for_bucket(bucket, omega);
			const auto log_p = TLinkPolicy::log_prob_for_bucket(bucket, omega);
			EXPECT_NEAR(log_p[1], std::log(p), 1e-10 * std::abs(std::log(p)) + 1e-15)
			    << "omega = " << w << ", bucket " << bucket;
			EXPECT_NEAR(log_p[0], std::log1p(-p), 1e-10 * std::abs(std::log1p(-p)) + 1e-15)
			    << "omega = " << w << ", bucket " << bucket;
		}
	}
}

//-----------------------------------
// The eight-state block
//-----------------------------------

TEST(FieldMath_Tests, state_index_is_a_bijection_onto_the_eight_states) {
	std::array<bool, 8> seen{};
	for (const bool y : {false, true}) {
		for (const bool z_s : {false, true}) {
			for (const bool z_m : {false, true}) {
				const size_t index = state_index(y, z_s, z_m);
				ASSERT_LT(index, 8u);
				EXPECT_FALSE(seen[index]) << "index " << index << " reached twice";
				seen[index] = true;
			}
		}
	}
}

TEST(FieldMath_Tests, block_probabilities_match_brute_force_enumeration) {
	const std::array<double, 2> lotus        = {0.7, 0.2};
	const std::array<double, 2> simple_error = {0.15, 0.85};
	for (const double w : {0.01, 0.1, 0.3, 0.45}) {
		for (const double p_s : {0.05, 0.5, 0.95}) {
			for (const double p_m : {0.1, 0.5, 0.8}) {
				const auto got = block_probabilities<TLinkPolicy>(
				    coretools::P(p_s), coretools::P(p_m), TErrorProbability(w),
				    as_probabilities(lotus), as_probabilities(simple_error));
				const auto want = brute_force_block(p_s, p_m, w, lotus, simple_error);
				for (size_t i = 0; i < 8; ++i) {
					EXPECT_NEAR(got[i], want[i], 1e-14)
					    << "state " << i << ", omega = " << w << ", p_s = " << p_s
					    << ", p_m = " << p_m;
				}
			}
		}
	}
}

TEST(FieldMath_Tests, block_probabilities_are_normalised) {
	const std::array<double, 2> lotus        = {0.9, 0.05};
	const std::array<double, 2> simple_error = {0.3, 0.7};
	for (const double w : {0.02, 0.2, 0.49}) {
		const auto block = block_probabilities<TLinkPolicy>(
		    coretools::P(0.3), coretools::P(0.6), TErrorProbability(w), as_probabilities(lotus),
		    as_probabilities(simple_error));
		double sum       = 0.0;
		for (const double p : block) {
			EXPECT_GE(p, 0.0) << "omega = " << w;
			sum += p;
		}
		EXPECT_NEAR(sum, 1.0, 1e-15) << "omega = " << w;
	}
}

TEST(FieldMath_Tests, the_block_reduces_to_the_link_when_the_data_says_nothing) {
	// With a flat data likelihood the field's conditional is the link alone, so summing the block
	// over the two tree fields must reproduce P(Y = 1) = a~_s * a~_m (ADR-0005, derivation 3).
	const std::array<double, 2> flat = {1.0, 1.0};
	for (const double w : {0.05, 0.25, 0.4}) {
		for (const double p_s : {0.2, 0.5, 0.9}) {
			for (const double p_m : {0.3, 0.7}) {
				const auto block = block_probabilities<TLinkPolicy>(
				    coretools::P(p_s), coretools::P(p_m), TErrorProbability(w),
				    as_probabilities(flat), as_probabilities(flat));
				double prob_y_is_one = 0.0;
				for (const bool z_s : {false, true}) {
					for (const bool z_m : {false, true}) {
						prob_y_is_one += block[state_index(true, z_s, z_m)];
					}
				}
				const double corrupted_s = w + (1.0 - 2.0 * w) * p_s;
				const double corrupted_m = w + (1.0 - 2.0 * w) * p_m;
				EXPECT_NEAR(prob_y_is_one, corrupted_s * corrupted_m, 1e-14)
				    << "omega = " << w << ", p_s = " << p_s << ", p_m = " << p_m;
			}
		}
	}
}

TEST(FieldMath_Tests, state_index_puts_the_field_in_the_high_bit) {
	// the documented convention: the four states sharing a field value are contiguous
	for (const bool z_s : {false, true}) {
		for (const bool z_m : {false, true}) {
			EXPECT_LT(state_index(false, z_s, z_m), 4u);
			EXPECT_GE(state_index(true, z_s, z_m), 4u);
		}
	}
}

TEST(FieldMath_Tests, the_block_leaves_the_tree_priors_alone_when_the_data_says_nothing) {
	// The link is a proper conditional density, so summing the block over the field must return
	// the two tree factors untouched. Shares no arithmetic with block_probabilities: it never
	// forms a link probability at all.
	const std::array<double, 2> flat = {1.0, 1.0};
	for (const double w : {0.03, 0.2, 0.48}) {
		for (const double p_s : {0.15, 0.5, 0.85}) {
			for (const double p_m : {0.25, 0.75}) {
				const auto block = block_probabilities<TLinkPolicy>(
				    coretools::P(p_s), coretools::P(p_m), TErrorProbability(w),
				    as_probabilities(flat), as_probabilities(flat));
				for (const bool z_s : {false, true}) {
					for (const bool z_m : {false, true}) {
						const double want = (z_s ? p_s : 1.0 - p_s) * (z_m ? p_m : 1.0 - p_m);
						const double got  = block[state_index(false, z_s, z_m)] +
						                   block[state_index(true, z_s, z_m)];
						EXPECT_NEAR(got, want, 1e-14)
						    << "omega = " << w << ", p_s = " << p_s << ", p_m = " << p_m
						    << ", Z_s = " << z_s << ", Z_m = " << z_m;
					}
				}
			}
		}
	}
}

TEST(FieldMath_Tests, the_block_tilts_the_link_by_the_data_odds) {
	// Conditioning on a tree field pair cancels both tree factors, leaving the link's odds
	// multiplied by the data's. A relation, not a recomputation -- the brute force above never
	// forms an odds ratio.
	const std::array<double, 2> lotus        = {0.6, 0.25};
	const std::array<double, 2> simple_error = {0.2, 0.9};
	for (const double w : {0.05, 0.2, 0.45}) {
		const TErrorProbability omega(w);
		const auto block = block_probabilities<TLinkPolicy>(
		    coretools::P(0.35), coretools::P(0.7), omega, as_probabilities(lotus),
		    as_probabilities(simple_error));
		for (const bool z_s : {false, true}) {
			for (const bool z_m : {false, true}) {
				const double link = TLinkPolicy::prob_y_is_one(z_s, z_m, omega);
				const double want =
				    (link * lotus[1] * simple_error[1]) / ((1.0 - link) * lotus[0] * simple_error[0]);
				const double got =
				    block[state_index(true, z_s, z_m)] / block[state_index(false, z_s, z_m)];
				EXPECT_NEAR(got, want, 1e-12 * want)
				    << "omega = " << w << ", Z_s = " << z_s << ", Z_m = " << z_m;
			}
		}
	}
}

// There is deliberately no test that the block rejects an argument outside [0, 1]. The arguments
// are coretools::Probability, so that guarantee belongs to the type and is tested in coretools.
// It is worth knowing that it holds only under CHECK_INTERVALS, which CMakeLists.txt defines for
// this test target and not for acol.

TEST(FieldMath_Tests, the_block_rejects_a_configuration_with_no_mass) {
	// zero is a perfectly good probability, so the type lets this through and the block itself
	// has to notice that nothing is left to normalise
	const std::array<double, 2> nothing = {0.0, 0.0};
	EXPECT_THROW(block_probabilities<TLinkPolicy>(coretools::P(0.5), coretools::P(0.5),
	                                              TErrorProbability(0.1), as_probabilities(nothing),
	                                              as_probabilities(nothing)),
	             std::invalid_argument);
}
