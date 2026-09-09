#include "coretools/Types/probability.h"
#include "field/TFieldMath.h"
#include "field/link_backend.h"
#include "gtest/gtest.h"
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
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

	bool operator==(const TCell &) const = default;
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

/// Where one of the eight states a leaf pair can be in sits in the fixture's own enumeration. The
/// kernel has no such index space: it draws the field cell alone, and each tree draws its own leaf.
constexpr size_t cell_index(bool y, bool z_s, bool z_m) {
	return (static_cast<size_t>(y) << 2U) | (static_cast<size_t>(z_s) << 1U) |
	       static_cast<size_t>(z_m);
}

/// The joint over one leaf pair, built from its definition -- the product of the two tree factors,
/// the link, and the data, normalised.
///
/// The sampler no longer forms this. It is kept as a fixture, because the conditional the field
/// update *does* draw from has to agree with it. The link inside it comes from the corruption
/// enumeration above rather than from any closed form, so the two sides share no arithmetic.
std::array<double, 8> brute_force_joint(double p_s, double p_m, double omega,
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
				weight[cell_index(y, z_s, z_m)] = w;
				total += w;
			}
		}
	}
	for (auto &w : weight) { w /= total; }
	return weight;
}

/// The field cell takes probabilities by type. The brute force stays on doubles, so that it shares
/// as little as possible with the thing it checks.
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

/// The bucket of a cell, written out here rather than taken from the link. A test of the
/// incremental counter arithmetic is worth nothing if it counts the way the kernel counts.
size_t bucket_written_out(const TCell &cell) {
	size_t ones = 0;
	if (cell.z_s) { ++ones; }
	if (cell.z_m) { ++ones; }
	return ones;
}

/// A configuration of `n` cells, drawn from a fixed-width generator so that it does not depend on
/// the platform.
std::vector<TCell> cells_from_seed(uint32_t seed, size_t n) {
	std::vector<TCell> cells;
	cells.reserve(n);
	for (size_t i = 0; i < n; ++i) {
		seed = seed * 1103515245u + 12345u;
		cells.push_back({.z_s = ((seed >> 16u) & 1u) != 0u,
		                 .z_m = ((seed >> 17u) & 1u) != 0u,
		                 .y   = ((seed >> 18u) & 1u) != 0u});
	}
	return cells;
}

/// The counters of a configuration, tallied the way the field update tallies them.
TLinkCounters counters_of(const std::vector<TCell> &cells) {
	TLinkCounters counters;
	for (const auto &c : cells) { counters.add(TLinkPolicy::bucket(c.z_s, c.z_m), c.y); }
	return counters;
}

/// The six counters recomputed from scratch: walk every cell and count it.
std::array<std::array<size_t, 2>, 3> recount(const std::vector<TCell> &cells) {
	std::array<std::array<size_t, 2>, 3> n{};
	for (const auto &cell : cells) { ++n[bucket_written_out(cell)][cell.y ? 1 : 0]; }
	return n;
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
	// The cast to void is what a [[nodiscard]] return needs when the call is the whole statement.
	// Every EXPECT_THROW below over a value-returning call carries it for the same reason.
	EXPECT_THROW(static_cast<void>(omega.for_tree(NUMBER_OF_TREES)), std::invalid_argument);
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
	// what a caller that moves one cell does: it knows its own old and new (bucket, field state)
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

TEST(FieldMath_Tests, merging_adds_one_tally_to_another) {
	// what the field update does after its parallel region: each thread counted its own share of
	// the cells, and the shares come back together
	TLinkCounters first;
	first.add(0, false);
	first.add(2, true);
	first.add(2, true);

	TLinkCounters second;
	second.add(2, true);
	second.add(1, false);

	first.merge(second);
	EXPECT_EQ(first.count(0, false), 1u);
	EXPECT_EQ(first.count(1, false), 1u);
	EXPECT_EQ(first.count(2, true), 3u);
	EXPECT_EQ(first.total(), 5u);
	// the tally merged in is left alone, so a thread's share can be read after it is committed
	EXPECT_EQ(second.total(), 2u);
}

TEST(FieldMath_Tests, merging_an_empty_tally_changes_nothing) {
	TLinkCounters counters;
	counters.add(1, true);
	counters.merge(TLinkCounters());
	EXPECT_EQ(counters.count(1, true), 1u);
	EXPECT_EQ(counters.total(), 1u);
}

TEST(FieldMath_Tests, six_counters_match_a_naive_per_cell_recomputation) {
	// a configuration with every bucket and both field states represented
	const auto cells    = cells_from_seed(1u, 400);
	const auto counters = counters_of(cells);
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
// The error probability's move (ADR-0005, derivation 2)
//-----------------------------------

TEST(FieldMath_Tests, the_move_reads_the_counters_and_no_cell) {
	// The whole point of the six counters: a proposal on the error probability is scored from them
	// alone. The configuration below is 400 cells, and the ratio never sees one.
	const auto cells    = cells_from_seed(7u, 400);
	const auto counters = counters_of(cells);

	for (const double old_w : omega_values()) {
		for (const double new_w : omega_values()) {
			const double want =
			    brute_force_log_likelihood(cells, new_w) - brute_force_log_likelihood(cells, old_w);
			const double got = TLinkPolicy::log_likelihood_ratio(counters, TErrorProbability(old_w),
			                                                     TErrorProbability(new_w));
			// The tolerance is the brute force's, as in the log-likelihood test above: it squares
			// and then takes a log where the closed form is affine in the bucket.
			EXPECT_NEAR(got, want, 1e-11 * std::abs(want) + 1e-9)
			    << "omega " << old_w << " -> " << new_w;
		}
	}
}

TEST(FieldMath_Tests, the_ratio_is_the_difference_of_the_two_log_likelihoods) {
	const auto counters = counters_from({{{400, 9600}, {1600, 8400}, {6400, 3600}}});
	const TErrorProbability old_omega(0.2);
	const TErrorProbability new_omega(0.05);

	EXPECT_DOUBLE_EQ(TLinkPolicy::log_likelihood_ratio(counters, old_omega, new_omega),
	                 TLinkPolicy::log_likelihood(counters, new_omega) -
	                     TLinkPolicy::log_likelihood(counters, old_omega));
}

TEST(FieldMath_Tests, a_proposal_that_does_not_move_has_a_zero_ratio) {
	const auto counters = counters_from({{{400, 9600}, {1600, 8400}, {6400, 3600}}});
	const TErrorProbability omega(0.13);
	EXPECT_DOUBLE_EQ(TLinkPolicy::log_likelihood_ratio(counters, omega, omega), 0.0);
}

TEST(FieldMath_Tests, an_empty_configuration_gives_a_zero_ratio) {
	// Before the first field update the counters are empty, so the error probability moves on its
	// prior alone.
	const TLinkCounters counters;
	EXPECT_DOUBLE_EQ(
	    TLinkPolicy::log_likelihood_ratio(counters, TErrorProbability(0.1), TErrorProbability(0.2)),
	    0.0);
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
// The diagnostic those constraints ship as
//-----------------------------------

TEST(FieldMath_Tests, the_diagnostic_reads_the_bucket_rates_off_the_counters) {
	const auto counters   = counters_from({{{400, 9600}, {1600, 8400}, {6400, 3600}}});
	const auto diagnostic = TLinkPolicy::diagnose(counters);

	EXPECT_TRUE(diagnostic.is_complete());
	for (size_t bucket = 0; bucket < TLinkCounters::n_buckets; ++bucket) {
		EXPECT_DOUBLE_EQ(diagnostic.prob[bucket], empirical_prob(counters, bucket))
		    << "bucket " << bucket;
	}
}

TEST(FieldMath_Tests, both_residuals_vanish_under_one_shared_error_probability) {
	// The same counts as the constraints test above: omega = 0.2 realised exactly.
	const auto diagnostic =
	    TLinkPolicy::diagnose(counters_from({{{400, 9600}, {1600, 8400}, {6400, 3600}}}));

	EXPECT_NEAR(diagnostic.and_identity_residual, 0.0, 1e-15);
	EXPECT_NEAR(diagnostic.shared_error_probability_residual, 0.0, 1e-15);
}

TEST(FieldMath_Tests, the_diagnostic_carries_the_identity_blind_spot_the_record_describes) {
	// Two trees corrupted at 0.05 and 0.2, mixed so that the pooled middle rate lands on the
	// geometric mean of the two mixed cells. The identity sees nothing; the shared-rate residual
	// misses 1 by the amount ADR-0005 quotes.
	const auto diagnostic =
	    TLinkPolicy::diagnose(counters_from({{{1000, 99000}, {8718, 91282}, {76000, 24000}}}));

	EXPECT_LT(std::abs(diagnostic.and_identity_residual), 1e-6);
	EXPECT_NEAR(diagnostic.shared_error_probability_residual, -0.0282202113, 1e-9);
}

TEST(FieldMath_Tests, an_empty_bucket_leaves_the_diagnostic_incomplete) {
	// A held field puts both tree fields at the field, so no leaf pair ever lands in bucket 1.
	// There is nothing to falsify then, and the diagnostic says so rather than reporting a number.
	const auto diagnostic = TLinkPolicy::diagnose(counters_from({{{0, 9600}, {0, 0}, {6400, 0}}}));

	EXPECT_FALSE(diagnostic.is_complete());
	EXPECT_TRUE(std::isnan(diagnostic.prob[1]));
	EXPECT_TRUE(std::isnan(diagnostic.and_identity_residual));
}

TEST(FieldMath_Tests, the_diagnostic_never_throws) {
	// It reports a finding about the model, so it must survive every configuration a chain can
	// reach -- an empty one included.
	EXPECT_NO_THROW(static_cast<void>(TLinkPolicy::diagnose(TLinkCounters())));
	EXPECT_NO_THROW(
	    static_cast<void>(TLinkPolicy::diagnose(counters_from({{{1, 0}, {0, 1}, {1, 1}}}))));
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
// One field cell
//-----------------------------------

TEST(FieldMath_Tests, a_field_cell_is_the_link_tilted_by_the_data) {
	// The definition: the link's probability of each field state, times what each source makes of
	// it, normalised. The link comes from the corruption enumeration, so this shares no closed form
	// with `prob_field_cell_is_one`.
	const std::array<double, 2> lotus        = {0.7, 0.2};
	const std::array<double, 2> simple_error = {0.15, 0.85};
	for (const double w : {0.01, 0.1, 0.3, 0.45}) {
		for (const bool z_s : {false, true}) {
			for (const bool z_m : {false, true}) {
				const double link     = brute_force_prob_y_is_one(z_s, z_m, w);
				const double weight_0 = (1.0 - link) * lotus[0] * simple_error[0];
				const double weight_1 = link * lotus[1] * simple_error[1];

				const auto got = prob_field_cell_is_one<TLinkPolicy>(
				    z_s, z_m, TErrorProbability(w), as_probabilities(lotus),
				    as_probabilities(simple_error));
				EXPECT_NEAR(got.get(), weight_1 / (weight_0 + weight_1), 1e-14)
				    << "omega = " << w << ", z_s = " << z_s << ", z_m = " << z_m;
			}
		}
	}
}

TEST(FieldMath_Tests, a_field_cell_with_no_data_follows_the_link) {
	// Both sources neutral, so the cell is the link and nothing else. This is the build that
	// compiled no data source in.
	const std::array<double, 2> neutral = {1.0, 1.0};
	for (const double w : {0.02, 0.2, 0.49}) {
		for (const bool z_s : {false, true}) {
			for (const bool z_m : {false, true}) {
				const auto got = prob_field_cell_is_one<TLinkPolicy>(z_s, z_m, TErrorProbability(w),
				                                                     as_probabilities(neutral),
				                                                     as_probabilities(neutral));
				EXPECT_NEAR(got.get(), TLinkPolicy::prob_y_is_one(z_s, z_m, TErrorProbability(w)),
				            1e-15)
				    << "omega = " << w << ", z_s = " << z_s << ", z_m = " << z_m;
			}
		}
	}
}

TEST(FieldMath_Tests, a_field_cell_is_the_joint_conditioned_on_the_two_tree_fields) {
	// A field cell drawn on its own must agree with the joint over the leaf pair. So the joint's
	// field marginal *given* the two tree field states is what a single-cell draw is, and the two
	// tree factors cancel from that conditional -- which is why the cell is not asked for them.
	//
	// The joint here is the fixture's own enumeration, built from the corruption events rather
	// than from any closed form, so the two sides share no arithmetic.
	const std::array<double, 2> lotus        = {0.6, 0.3};
	const std::array<double, 2> simple_error = {0.25, 0.75};
	for (const double w : {0.05, 0.25, 0.4}) {
		for (const double p_s : {0.2, 0.5, 0.9}) {
			for (const double p_m : {0.15, 0.5, 0.85}) {
				const auto joint = brute_force_joint(p_s, p_m, w, lotus, simple_error);
				for (const bool z_s : {false, true}) {
					for (const bool z_m : {false, true}) {
						const double at_zero = joint[cell_index(false, z_s, z_m)];
						const double at_one  = joint[cell_index(true, z_s, z_m)];
						const auto got       = prob_field_cell_is_one<TLinkPolicy>(
						    z_s, z_m, TErrorProbability(w), as_probabilities(lotus),
						    as_probabilities(simple_error));
						EXPECT_NEAR(got.get(), at_one / (at_zero + at_one), 1e-14)
						    << "omega = " << w << ", p_s = " << p_s << ", p_m = " << p_m
						    << ", z_s = " << z_s << ", z_m = " << z_m;
					}
				}
			}
		}
	}
}

TEST(FieldMath_Tests, a_field_cell_rejects_a_configuration_with_no_mass) {
	// zero is a perfectly good probability, so the type lets this through and the cell itself has
	// to notice that nothing is left to normalise
	const std::array<double, 2> nothing = {0.0, 0.0};
	EXPECT_THROW(static_cast<void>(prob_field_cell_is_one<TLinkPolicy>(
	                 true, true, TErrorProbability(0.1), as_probabilities(nothing),
	                 as_probabilities(nothing))),
	             std::invalid_argument);
}
