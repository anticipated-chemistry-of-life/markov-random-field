#include "coretools/Main/TRandomGenerator.h"
#include "mass_spec/feature_likelihood.h"
#include "mass_spec/msms_run.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <vector>

namespace {

/// A small run with overlapping candidate lists so every move type is reachable:
///   feature 0 candidates: {10, 20}
///   feature 1 candidates: {10, 20, 30}
///   feature 2 candidates: {20, 30}
/// All features start assigned to the unknown molecule.
TMassSpecRun make_run() {
	TMassSpecRun run;
	std::vector<TFeatureLikelihood> f0 = {TFeatureLikelihood(10, 100), TFeatureLikelihood(20, 120)};
	std::vector<TFeatureLikelihood> f1 = {TFeatureLikelihood(10, 90), TFeatureLikelihood(20, 80),
	                                      TFeatureLikelihood(30, 70)};
	std::vector<TFeatureLikelihood> f2 = {TFeatureLikelihood(20, 60), TFeatureLikelihood(30, 50)};
	run.add_likelihood_vector(f0);
	run.add_likelihood_vector(f1);
	run.add_likelihood_vector(f2);
	run.set_probabilities_of_unknowns({50, 60, 70});
	run.initialize_assignments_to_unknown();
	return run;
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

} // namespace

TEST(TAssignmentUpdate_Tests, initialize_to_unknown) {
	auto run = make_run();
	ASSERT_EQ(run.number_of_assignments(), 3u);
	for (size_t f = 0; f < 3; ++f) {
		EXPECT_TRUE(run.get_current_assignment(f).is_unknown_molecule());
		EXPECT_EQ(run.get_current_assignment(f).get_binned_likelihood(), run.probability_of_unknown(f));
	}
	EXPECT_FALSE(run.is_molecule_assigned(10));
	EXPECT_TRUE(no_molecule_assigned_twice(run));
}

TEST(TAssignmentUpdate_Tests, initialize_requires_matching_unknown_sizes) {
	TMassSpecRun run;
	std::vector<TFeatureLikelihood> f0 = {TFeatureLikelihood(10, 100)};
	run.add_likelihood_vector(f0);
	run.set_probabilities_of_unknowns({1, 2}); // 2 unknowns but only 1 feature
	EXPECT_ANY_THROW(run.initialize_assignments_to_unknown());
}

TEST(TAssignmentUpdate_Tests, apply_revert_to_unknown) {
	auto run = make_run();
	run.set_current_assignment(0, TFeatureLikelihood(10, 100));
	ASSERT_TRUE(run.is_molecule_assigned(10));

	TAssignmentProposal p;
	p.type      = AssignmentMoveType::ToUnknown;
	p.feature_a = 0;
	p.old_a     = run.get_current_assignment(0);
	p.new_a     = TFeatureLikelihood::new_unknown_molecule(run.probability_of_unknown(0));

	run.apply_move(p);
	EXPECT_TRUE(run.get_current_assignment(0).is_unknown_molecule());
	EXPECT_FALSE(run.is_molecule_assigned(10));

	run.revert_move(p);
	EXPECT_FALSE(run.get_current_assignment(0).is_unknown_molecule());
	EXPECT_EQ(run.get_current_assignment(0).get_molecule_index(), 10u);
	EXPECT_TRUE(run.is_molecule_assigned(10));
}

TEST(TAssignmentUpdate_Tests, swap_apply_revert) {
	auto run = make_run();
	run.set_current_assignment(0, TFeatureLikelihood(10, 100));
	run.set_current_assignment(1, TFeatureLikelihood(20, 80));

	TAssignmentProposal p;
	p.type      = AssignmentMoveType::Swap;
	p.feature_a = 0;
	p.old_a     = run.get_current_assignment(0);
	p.new_a     = TFeatureLikelihood(20, 120); // molecule 20 is candidate of feature 0 (binned 120)
	p.feature_b = 1;
	p.old_b     = run.get_current_assignment(1);
	p.new_b     = TFeatureLikelihood(10, 90); // molecule 10 is candidate of feature 1 (binned 90)

	run.apply_move(p);
	EXPECT_EQ(run.get_current_assignment(0).get_molecule_index(), 20u);
	EXPECT_EQ(run.get_current_assignment(1).get_molecule_index(), 10u);
	EXPECT_TRUE(no_molecule_assigned_twice(run));

	run.revert_move(p);
	EXPECT_EQ(run.get_current_assignment(0).get_molecule_index(), 10u);
	EXPECT_EQ(run.get_current_assignment(1).get_molecule_index(), 20u);
}

// Property test: starting from a valid state, proposing-and-applying any number of moves must never
// break the "a molecule is assigned to at most one feature" invariant.
TEST(TAssignmentUpdate_Tests, propose_move_preserves_invariant) {
	coretools::instances::randomGenerator().setSeed(42, true);
	auto run = make_run();
	run.set_current_assignment(0, TFeatureLikelihood(10, 100));
	run.set_current_assignment(1, TFeatureLikelihood(20, 80));
	ASSERT_TRUE(no_molecule_assigned_twice(run));

	int n_valid = 0;
	for (int it = 0; it < 3000; ++it) {
		const auto p = run.propose_move();
		if (!p.is_valid()) { continue; }
		++n_valid;
		run.apply_move(p);
		ASSERT_TRUE(no_molecule_assigned_twice(run))
		    << "invariant broken after move type " << static_cast<int>(p.type);
		// exercise revert on roughly half the accepted moves
		if (coretools::instances::randomGenerator().pickOneOfTwo()) {
			run.revert_move(p);
			ASSERT_TRUE(no_molecule_assigned_twice(run));
		}
	}
	EXPECT_GT(n_valid, 0); // the mixed state should make at least some moves constructible
}

TEST(TAssignmentUpdate_Tests, propose_move_invalid_when_no_candidates) {
	TMassSpecRun run;
	std::vector<TFeatureLikelihood> empty_feature; // a feature with no candidate molecules
	run.add_likelihood_vector(empty_feature);
	run.set_probabilities_of_unknowns({42});
	run.initialize_assignments_to_unknown();
	for (int i = 0; i < 50; ++i) { EXPECT_FALSE(run.propose_move().is_valid()); }
}
