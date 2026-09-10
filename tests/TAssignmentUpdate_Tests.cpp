#include "coretools/Main/TRandomGenerator.h"
#include "mass_spec/feature_likelihood.h"
#include "mass_spec/msms_run.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace {

/// Enough molecules for every candidate index used below.
constexpr size_t NUMBER_OF_MOLECULES = 40;

/// Records every proposal it is asked to score and returns a constant log-likelihood ratio.
/// `+ACCEPT_ALL` makes the update accept everything, `-ACCEPT_ALL` makes it reject everything while
/// still going through the real Metropolis-Hastings acceptance (the value stays finite, so the
/// non-finite shortcut in `evalLogH` is not what does the rejecting).
class TFixedScorer final : public TAssignmentScorer {
public:
	double value = 0.0;
	mutable std::vector<TAssignmentProposal> seen;

	[[nodiscard]] double log_likelihood_ratio(const TAssignmentProposal &proposal) const override {
		seen.push_back(proposal);
		return value;
	}
};

constexpr double ACCEPT_ALL = 1e9;

/// Builds a run from one candidate list per feature and one unknown probability per feature.
TMassSpecRun make_run(std::vector<std::vector<TFeatureLikelihood>> features,
                      const std::vector<uint8_t> &probabilities_of_unknowns) {
	TMassSpecRun run;
	for (auto &feature : features) { run.add_likelihood_vector(feature); }
	run.set_probabilities_of_unknowns(probabilities_of_unknowns);
	run.initialize_assignments_to_unknown();
	run.set_number_of_molecules(NUMBER_OF_MOLECULES);
	return run;
}

/// The run used by most tests, with overlapping candidate lists so every move type is reachable:
///   feature 0 candidates: {10, 20}
///   feature 1 candidates: {10, 20, 30}
///   feature 2 candidates: {20, 30}
/// All features start assigned to the unknown molecule.
TMassSpecRun make_overlapping_run() {
	return make_run(
	    {{TFeatureLikelihood(10, 100), TFeatureLikelihood(20, 120)},
	     {TFeatureLikelihood(10, 90), TFeatureLikelihood(20, 80), TFeatureLikelihood(30, 70)},
	     {TFeatureLikelihood(20, 60), TFeatureLikelihood(30, 50)}},
	    {50, 60, 70});
}

/// The invariant the moves must preserve: no real molecule is assigned to two features at once.
bool no_molecule_assigned_twice(const TMassSpecRun &run) {
	std::vector<uint32_t> assigned;
	for (size_t f = 0; f < run.number_of_assignments(); ++f) {
		const auto &a = run.get_current_assignment(f);
		if (!a.is_unknown_molecule()) { assigned.push_back(a.get_molecule_index()); }
	}
	std::sort(assigned.begin(), assigned.end());
	return std::adjacent_find(assigned.begin(), assigned.end()) == assigned.end();
}

/// The second invariant: a feature may only ever hold one of its own candidates, and it must carry
/// that candidate's binned likelihood.
bool every_assignment_is_an_own_candidate(const TMassSpecRun &run) {
	for (size_t f = 0; f < run.number_of_assignments(); ++f) {
		const auto &a = run.get_current_assignment(f);
		if (a.is_unknown_molecule()) { continue; }
		const auto hit = run.is_molecule_in_feature(f, a.get_molecule_index());
		if (!hit.found || *hit.binned_likelihood != a.get_binned_likelihood()) { return false; }
	}
	return true;
}

std::vector<TFeatureLikelihood> snapshot(const TMassSpecRun &run) {
	std::vector<TFeatureLikelihood> assignments;
	assignments.reserve(run.number_of_assignments());
	for (size_t f = 0; f < run.number_of_assignments(); ++f) {
		assignments.push_back(run.get_current_assignment(f));
	}
	return assignments;
}

bool same_assignments(const std::vector<TFeatureLikelihood> &left,
                      const std::vector<TFeatureLikelihood> &right) {
	if (left.size() != right.size()) { return false; }
	for (size_t f = 0; f < left.size(); ++f) {
		if (left[f].get_molecule_index() != right[f].get_molecule_index()) { return false; }
		if (left[f].get_binned_likelihood() != right[f].get_binned_likelihood()) { return false; }
	}
	return true;
}

size_t count_of_type(const TFixedScorer &scorer, AssignmentMoveType type) {
	return (size_t)std::count_if(scorer.seen.begin(), scorer.seen.end(),
	                             [type](const TAssignmentProposal &p) { return p.type == type; });
}

} // namespace

TEST(TAssignmentUpdate_Tests, update_reject_all_leaves_state_unchanged) {
	coretools::instances::randomGenerator().setSeed(42, true);
	auto run = make_overlapping_run();
	run.set_current_assignment(1, TFeatureLikelihood(20, 80));
	const auto before = snapshot(run);

	TFixedScorer scorer;
	scorer.value = -ACCEPT_ALL;
	run.update_all_assignments(scorer, 0.1);

	EXPECT_FALSE(scorer.seen.empty());
	EXPECT_TRUE(same_assignments(before, snapshot(run)));
}

TEST(TAssignmentUpdate_Tests, update_accept_all_preserves_invariant) {
	coretools::instances::randomGenerator().setSeed(42, true);
	auto run = make_overlapping_run();

	TFixedScorer scorer;
	scorer.value = ACCEPT_ALL;
	for (int update = 0; update < 200; ++update) {
		run.update_all_assignments(scorer, 0.1);
		ASSERT_TRUE(no_molecule_assigned_twice(run)) << "after update " << update;
		ASSERT_TRUE(every_assignment_is_an_own_candidate(run)) << "after update " << update;
	}
	// with everything accepted the run must actually have moved away from the all-unknown start
	EXPECT_FALSE(scorer.seen.empty());
}

// The features are drawn at random, so an update makes as many attempts as there are features but
// visits them in no particular order and may repeat one. Over many updates every feature must come
// up.
TEST(TAssignmentUpdate_Tests, update_makes_one_attempt_per_feature) {
	coretools::instances::randomGenerator().setSeed(42, true);
	// one candidate per feature, so every drawn feature yields exactly one proposal
	auto run = make_run({{TFeatureLikelihood(10, 100)},
	                     {TFeatureLikelihood(11, 100)},
	                     {TFeatureLikelihood(12, 100)},
	                     {TFeatureLikelihood(13, 100)}},
	                    {50, 50, 50, 50});

	TFixedScorer scorer;
	scorer.value = -ACCEPT_ALL;
	run.update_all_assignments(scorer, 0.1);
	EXPECT_EQ(scorer.seen.size(), 4u);

	for (int update = 0; update < 50; ++update) { run.update_all_assignments(scorer, 0.1); }
	EXPECT_EQ(scorer.seen.size(), 4u * 51u);

	std::vector<bool> visited(4, false);
	for (const auto &proposal : scorer.seen) {
		ASSERT_LT(proposal.feature_a, 4u);
		visited[proposal.feature_a] = true;
	}
	for (size_t f = 0; f < 4; ++f) { EXPECT_TRUE(visited[f]) << "feature " << f << " never drawn"; }
}

TEST(TAssignmentUpdate_Tests, hastings_from_unknown) {
	coretools::instances::randomGenerator().setSeed(42, true);
	// a single feature with 3 free candidates, sitting at the unknown molecule
	auto run = make_run(
	    {{TFeatureLikelihood(10, 100), TFeatureLikelihood(20, 110), TFeatureLikelihood(30, 120)}},
	    {50});

	TFixedScorer scorer;
	scorer.value          = -ACCEPT_ALL;
	constexpr double beta = 0.25;
	run.update_all_assignments(scorer, beta);

	ASSERT_EQ(scorer.seen.size(), 1u);
	EXPECT_EQ(scorer.seen[0].type, AssignmentMoveType::FromUnknown);
	EXPECT_NEAR(scorer.seen[0].log_hastings, std::log(beta * 3.0), 1e-12);
}

TEST(TAssignmentUpdate_Tests, hastings_to_unknown) {
	coretools::instances::randomGenerator().setSeed(42, true);
	auto run = make_overlapping_run();
	run.set_current_assignment(0, TFeatureLikelihood(10, 100)); // c_f = 2
	run.set_current_assignment(1, TFeatureLikelihood(20, 80));  // c_f = 3
	run.set_current_assignment(2, TFeatureLikelihood(30, 50));  // c_f = 2
	const std::vector<double> candidate_counts = {2.0, 3.0, 2.0};

	TFixedScorer scorer;
	scorer.value = -ACCEPT_ALL;
	// beta == 1 forces the to-unknown branch for every assigned feature
	run.update_all_assignments(scorer, 1.0);

	ASSERT_EQ(scorer.seen.size(), 3u);
	for (const auto &proposal : scorer.seen) {
		EXPECT_EQ(proposal.type, AssignmentMoveType::ToUnknown);
		EXPECT_NEAR(proposal.log_hastings, -std::log(candidate_counts[proposal.feature_a]), 1e-12);
	}
}

TEST(TAssignmentUpdate_Tests, hastings_zero_for_move_to_free) {
	coretools::instances::randomGenerator().setSeed(42, true);
	// one feature holding 10, with 20 free among its candidates
	auto run = make_run({{TFeatureLikelihood(10, 100), TFeatureLikelihood(20, 120)}}, {50});
	run.set_current_assignment(0, TFeatureLikelihood(10, 100));

	TFixedScorer scorer;
	scorer.value = -ACCEPT_ALL;
	// beta == 0 removes the to-unknown branch, so only the pick-a-molecule branch is exercised
	for (int update = 0; update < 50; ++update) { run.update_all_assignments(scorer, 0.0); }

	EXPECT_GT(count_of_type(scorer, AssignmentMoveType::MoveToFree), 0u);
	for (const auto &proposal : scorer.seen) {
		EXPECT_EQ(proposal.type, AssignmentMoveType::MoveToFree);
		EXPECT_EQ(proposal.log_hastings, 0.0);
	}
}

TEST(TAssignmentUpdate_Tests, hastings_zero_for_swap) {
	coretools::instances::randomGenerator().setSeed(42, true);
	// both features have {10, 20} as candidates, so the swap between them is always legal
	auto run = make_run({{TFeatureLikelihood(10, 100), TFeatureLikelihood(20, 120)},
	                     {TFeatureLikelihood(10, 90), TFeatureLikelihood(20, 80)}},
	                    {50, 60});
	run.set_current_assignment(0, TFeatureLikelihood(10, 100));
	run.set_current_assignment(1, TFeatureLikelihood(20, 80));

	TFixedScorer scorer;
	scorer.value = -ACCEPT_ALL;
	for (int update = 0; update < 50; ++update) { run.update_all_assignments(scorer, 0.0); }

	EXPECT_GT(count_of_type(scorer, AssignmentMoveType::Swap), 0u);
	for (const auto &proposal : scorer.seen) {
		EXPECT_EQ(proposal.type, AssignmentMoveType::Swap);
		EXPECT_EQ(proposal.log_hastings, 0.0);
	}
}

TEST(TAssignmentUpdate_Tests, hastings_swap_with_unknown) {
	coretools::instances::randomGenerator().setSeed(42, true);
	// feature 0 has 2 candidates and sits at the unknown molecule, feature 1 has 3 candidates and
	// holds molecule 20, which is also a candidate of feature 0. No other feature has molecule 20
	// as a candidate, so feature 0 is the only possible driver of this move.
	auto run = make_run(
	    {{TFeatureLikelihood(10, 100), TFeatureLikelihood(20, 120)},
	     {TFeatureLikelihood(10, 90), TFeatureLikelihood(20, 80), TFeatureLikelihood(30, 70)}},
	    {50, 60});
	run.set_current_assignment(1, TFeatureLikelihood(20, 80));

	TFixedScorer scorer;
	scorer.value = -ACCEPT_ALL;
	for (int update = 0; update < 50; ++update) { run.update_all_assignments(scorer, 0.0); }

	size_t n_swaps_with_unknown = 0;
	for (const auto &proposal : scorer.seen) {
		if (proposal.type != AssignmentMoveType::SwapWithUnknown) { continue; }
		++n_swaps_with_unknown;
		EXPECT_EQ(proposal.feature_a, 0u);
		EXPECT_EQ(proposal.feature_b, 1u);
		EXPECT_NEAR(proposal.log_hastings, std::log(2.0 / 3.0), 1e-12);
		// the arriving feature carries its own likelihood for molecule 20, not the other one's
		EXPECT_EQ(proposal.new_a.get_molecule_index(), 20u);
		EXPECT_EQ(proposal.new_a.get_binned_likelihood(), 120u);
		EXPECT_EQ(proposal.old_b.get_molecule_index(), 20u);
		EXPECT_EQ(proposal.old_b.get_binned_likelihood(), 80u);
		EXPECT_TRUE(proposal.new_b.is_unknown_molecule());
	}
	EXPECT_GT(n_swaps_with_unknown, 0u);
}

TEST(TAssignmentUpdate_Tests, illegal_swap_is_skipped) {
	coretools::instances::randomGenerator().setSeed(42, true);
	// molecule 20 is the only candidate the two features share, so no legal swap exists: whichever
	// feature drives it, the molecule it would give away is not a candidate of the other feature
	auto run = make_run({{TFeatureLikelihood(10, 100), TFeatureLikelihood(20, 120)},
	                     {TFeatureLikelihood(20, 80), TFeatureLikelihood(30, 70)}},
	                    {50, 60});
	run.set_current_assignment(0, TFeatureLikelihood(10, 100));
	run.set_current_assignment(1, TFeatureLikelihood(20, 80));

	TFixedScorer scorer;
	scorer.value = ACCEPT_ALL;
	for (int update = 0; update < 200; ++update) {
		run.update_all_assignments(scorer, 0.0);
		ASSERT_TRUE(no_molecule_assigned_twice(run)) << "after update " << update;
	}
	EXPECT_EQ(count_of_type(scorer, AssignmentMoveType::Swap), 0u);
}

// The holder map is rebuilt at the start of every update and throws if it finds a molecule held by
// two features, so a bookkeeping bug in one update surfaces at the start of the next one.
TEST(TAssignmentUpdate_Tests, holder_map_stays_in_sync_across_updates) {
	coretools::instances::randomGenerator().setSeed(42, true);
	auto run = make_overlapping_run();

	TFixedScorer scorer;
	scorer.value = ACCEPT_ALL;
	for (int update = 0; update < 50; ++update) {
		EXPECT_NO_THROW(run.update_all_assignments(scorer, 0.1)) << "at update " << update;
	}
}

TEST(TAssignmentUpdate_Tests, update_on_run_with_zero_features_is_noop) {
	TMassSpecRun run;
	run.set_probabilities_of_unknowns({});
	run.initialize_assignments_to_unknown();
	run.set_number_of_molecules(NUMBER_OF_MOLECULES);

	TFixedScorer scorer;
	EXPECT_NO_THROW(run.update_all_assignments(scorer, 0.1));
	EXPECT_TRUE(scorer.seen.empty());
}

TEST(TAssignmentUpdate_Tests, update_throws_when_number_of_molecules_unset) {
	TMassSpecRun run;
	std::vector<TFeatureLikelihood> f0 = {TFeatureLikelihood(10, 100)};
	run.add_likelihood_vector(f0);
	run.set_probabilities_of_unknowns({50});
	run.initialize_assignments_to_unknown();
	// set_number_of_molecules() deliberately not called

	TFixedScorer scorer;
	EXPECT_ANY_THROW(run.update_all_assignments(scorer, 0.1));
}

TEST(TAssignmentUpdate_Tests, update_throws_when_assignments_were_not_initialized) {
	TMassSpecRun run;
	std::vector<TFeatureLikelihood> f0 = {TFeatureLikelihood(10, 100)};
	run.add_likelihood_vector(f0);
	run.set_number_of_molecules(NUMBER_OF_MOLECULES);
	// initialize_assignments_to_unknown() deliberately not called

	TFixedScorer scorer;
	EXPECT_ANY_THROW(run.update_all_assignments(scorer, 0.1));
}

// The unknown molecule has a sentinel index far outside the molecule space, so the holder map must
// skip unknown assignments rather than index itself with it.
TEST(TAssignmentUpdate_Tests, update_does_not_throw_when_all_features_unknown) {
	coretools::instances::randomGenerator().setSeed(42, true);
	auto run = make_overlapping_run();
	for (size_t f = 0; f < run.number_of_assignments(); ++f) {
		ASSERT_TRUE(run.get_current_assignment(f).is_unknown_molecule());
	}

	TFixedScorer scorer;
	scorer.value = -ACCEPT_ALL;
	EXPECT_NO_THROW(run.update_all_assignments(scorer, 0.1));
}

TEST(TAssignmentUpdate_Tests, feature_with_no_candidates_is_skipped) {
	coretools::instances::randomGenerator().setSeed(42, true);
	auto run          = make_run({{}}, {42});
	const auto before = snapshot(run);

	TFixedScorer scorer;
	scorer.value = ACCEPT_ALL;
	for (int update = 0; update < 50; ++update) { run.update_all_assignments(scorer, 0.5); }

	EXPECT_TRUE(scorer.seen.empty());
	EXPECT_TRUE(same_assignments(before, snapshot(run)));
}

TEST(TAssignmentUpdate_Tests, beta_zero_never_proposes_to_unknown) {
	coretools::instances::randomGenerator().setSeed(42, true);
	auto run = make_overlapping_run();
	run.set_current_assignment(0, TFeatureLikelihood(10, 100));
	run.set_current_assignment(1, TFeatureLikelihood(20, 80));
	run.set_current_assignment(2, TFeatureLikelihood(30, 50));

	TFixedScorer scorer;
	scorer.value = -ACCEPT_ALL;
	for (int update = 0; update < 50; ++update) { run.update_all_assignments(scorer, 0.0); }

	EXPECT_FALSE(scorer.seen.empty());
	EXPECT_EQ(count_of_type(scorer, AssignmentMoveType::ToUnknown), 0u);
}

TEST(TAssignmentUpdate_Tests, beta_one_always_proposes_to_unknown) {
	coretools::instances::randomGenerator().setSeed(42, true);
	auto run = make_overlapping_run();
	run.set_current_assignment(0, TFeatureLikelihood(10, 100));
	run.set_current_assignment(1, TFeatureLikelihood(20, 80));
	run.set_current_assignment(2, TFeatureLikelihood(30, 50));

	TFixedScorer scorer;
	scorer.value = -ACCEPT_ALL;
	run.update_all_assignments(scorer, 1.0);

	ASSERT_EQ(scorer.seen.size(), 3u);
	EXPECT_EQ(count_of_type(scorer, AssignmentMoveType::ToUnknown), 3u);
}

TEST(TAssignmentUpdate_Tests, identity_pick_is_skipped) {
	coretools::instances::randomGenerator().setSeed(42, true);
	// the only candidate of the feature is the molecule it already holds
	auto run = make_run({{TFeatureLikelihood(10, 100)}}, {50});
	run.set_current_assignment(0, TFeatureLikelihood(10, 100));

	TFixedScorer scorer;
	scorer.value = ACCEPT_ALL;
	for (int update = 0; update < 50; ++update) { run.update_all_assignments(scorer, 0.0); }

	EXPECT_TRUE(scorer.seen.empty());
	EXPECT_EQ(run.get_current_assignment(0).get_molecule_index(), 10u);
}
