#include "coretools/Main/TRandomGenerator.h"
#include "mass_spec/feature_likelihood.h"
#include "mass_spec/msms_run.h"
#include "gtest/gtest.h"
#include <cmath>
#include <cstdint>
#include <map>
#include <vector>

// Checks that the sweep really samples from the distribution the scorer describes, which is the
// only thing that can catch a wrong Hastings ratio: the invariant tests in TAssignmentSweep_Tests
// pass just as happily with the ratios set to zero.
//
// The state space is kept small enough to enumerate. Two features with different numbers of
// candidates (2 and 3) so that the unknown-swap, whose Hastings ratio is log(c_f / c_f'), is not
// trivially 1:
//   feature 0 candidates: {10, 20}
//   feature 1 candidates: {10, 20, 30}
// The target gives every (feature, assignment) pair a weight and multiplies them, so the target
// probability of a state is the product of its two weights, normalized over the states in which no
// molecule is held twice.

namespace {

constexpr size_t NUMBER_OF_MOLECULES = 40;
constexpr uint32_t UNKNOWN           = 0; // stands in for the unknown molecule in the state key

/// Unnormalized weight of assigning `molecule` (or the unknown molecule, coded as UNKNOWN) to
/// feature `f`.
double target_weight(size_t f, uint32_t molecule) {
	static const std::map<std::pair<size_t, uint32_t>, double> weights = {
	    {{0, UNKNOWN}, 1.0}, {{0, 10}, 2.0}, {{0, 20}, 3.0}, {{1, UNKNOWN}, 1.0},
	    {{1, 10}, 4.0},      {{1, 20}, 1.5}, {{1, 30}, 2.5}};
	return weights.at({f, molecule});
}

uint32_t state_of(const TFeatureLikelihood &assignment) {
	return assignment.is_unknown_molecule() ? UNKNOWN : assignment.get_molecule_index();
}

/// Scores a proposal against `target_weight`, exactly as TMSMSData scores it against the real
/// model: log P(new state) - log P(old state), with everything that does not change cancelling.
class TWeightScorer final : public TAssignmentScorer {
public:
	[[nodiscard]] double log_likelihood_ratio(const TAssignmentProposal &p) const override {
		double log_ratio = std::log(target_weight(p.feature_a, state_of(p.new_a))) -
		                   std::log(target_weight(p.feature_a, state_of(p.old_a)));
		if (p.touches_two_features()) {
			log_ratio += std::log(target_weight(p.feature_b, state_of(p.new_b))) -
			             std::log(target_weight(p.feature_b, state_of(p.old_b)));
		}
		return log_ratio;
	}
};

using State = std::pair<uint32_t, uint32_t>;

/// All states in which no real molecule is held by both features, with their target probabilities.
std::map<State, double> analytic_distribution() {
	const std::vector<uint32_t> candidates_0 = {UNKNOWN, 10, 20};
	const std::vector<uint32_t> candidates_1 = {UNKNOWN, 10, 20, 30};

	std::map<State, double> distribution;
	double total = 0.0;
	for (const uint32_t a0 : candidates_0) {
		for (const uint32_t a1 : candidates_1) {
			if (a0 != UNKNOWN && a0 == a1) { continue; } // a molecule cannot be held twice
			const double weight    = target_weight(0, a0) * target_weight(1, a1);
			distribution[{a0, a1}] = weight;
			total += weight;
		}
	}
	for (auto &[state, weight] : distribution) { weight /= total; }
	return distribution;
}

TMassSpecRun make_run() {
	TMassSpecRun run;
	std::vector<TFeatureLikelihood> f0 = {TFeatureLikelihood(10, 100), TFeatureLikelihood(20, 120)};
	std::vector<TFeatureLikelihood> f1 = {TFeatureLikelihood(10, 90), TFeatureLikelihood(20, 80),
	                                      TFeatureLikelihood(30, 70)};
	run.add_likelihood_vector(f0);
	run.add_likelihood_vector(f1);
	run.set_probabilities_of_unknowns({50, 60});
	run.initialize_assignments_to_unknown();
	run.set_number_of_molecules(NUMBER_OF_MOLECULES);
	return run;
}

} // namespace

TEST(TAssignmentStationary_Tests, sweep_samples_the_target_distribution) {
	coretools::instances::randomGenerator().setSeed(42, true);
	auto run = make_run();
	const TWeightScorer scorer;

	constexpr size_t n_sweeps = 400000;
	constexpr size_t n_burnin = 1000;
	std::map<State, size_t> counts;
	for (size_t sweep = 0; sweep < n_sweeps + n_burnin; ++sweep) {
		run.update_all_assignments(scorer, 0.3);
		if (sweep < n_burnin) { continue; }
		++counts[{state_of(run.get_current_assignment(0)),
		          state_of(run.get_current_assignment(1))}];
	}

	const auto expected = analytic_distribution();
	ASSERT_EQ(counts.size(), expected.size()) << "the chain did not visit every reachable state";
	for (const auto &[state, probability] : expected) {
		const double observed = (double)counts.at(state) / (double)n_sweeps;
		// the smallest target probability here is ~2%, so an absolute tolerance of 0.005 is a
		// relative error of at most ~25% on the rarest state and much tighter on the common ones
		EXPECT_NEAR(observed, probability, 0.005)
		    << "state (" << state.first << ", " << state.second << ")";
	}
}
