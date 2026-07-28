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
		EXPECT_EQ(run.get_current_assignment(f).get_binned_likelihood(),
		          run.probability_of_unknown(f));
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

TEST(TAssignmentUpdate_Tests, swap_with_unknown_apply_revert) {
	auto run = make_run();
	run.set_current_assignment(1, TFeatureLikelihood(20, 80)); // feature 1 holds molecule 20
	ASSERT_TRUE(run.get_current_assignment(0).is_unknown_molecule());

	// feature 0 takes molecule 20 (its own binned likelihood for it is 120, not feature 1's 80) and
	// pushes feature 1 back to the unknown molecule
	TAssignmentProposal p;
	p.type      = AssignmentMoveType::SwapWithUnknown;
	p.feature_a = 0;
	p.old_a     = run.get_current_assignment(0);
	p.new_a     = TFeatureLikelihood(20, 120);
	p.feature_b = 1;
	p.old_b     = run.get_current_assignment(1);
	p.new_b     = TFeatureLikelihood::new_unknown_molecule(run.probability_of_unknown(1));

	run.apply_move(p);
	EXPECT_EQ(run.get_current_assignment(0).get_molecule_index(), 20u);
	EXPECT_EQ(run.get_current_assignment(0).get_binned_likelihood(), 120u);
	EXPECT_TRUE(run.get_current_assignment(1).is_unknown_molecule());
	EXPECT_TRUE(no_molecule_assigned_twice(run));

	run.revert_move(p);
	EXPECT_TRUE(run.get_current_assignment(0).is_unknown_molecule());
	EXPECT_EQ(run.get_current_assignment(1).get_molecule_index(), 20u);
	EXPECT_EQ(run.get_current_assignment(1).get_binned_likelihood(), 80u);
}

// Guards against a future move type that changes two features but is missed by one of the three
// gates that dispatch on this predicate (apply_move, revert_move and the likelihood ratio).
TEST(TAssignmentUpdate_Tests, touches_two_features_covers_both_swap_types) {
	TAssignmentProposal p;
	p.type = AssignmentMoveType::Swap;
	EXPECT_TRUE(p.touches_two_features());
	p.type = AssignmentMoveType::SwapWithUnknown;
	EXPECT_TRUE(p.touches_two_features());

	for (const auto type : {AssignmentMoveType::ToUnknown, AssignmentMoveType::FromUnknown,
	                        AssignmentMoveType::MoveToFree, AssignmentMoveType::Invalid}) {
		p.type = type;
		EXPECT_FALSE(p.touches_two_features()) << "for move type " << static_cast<int>(type);
	}
}

TEST(TAssignmentUpdate_Tests, add_likelihood_vector_rejects_zero_binned_likelihood) {
	TMassSpecRun run;
	// binned likelihood 0 maps to probability 0, which no assignment could ever have
	std::vector<TFeatureLikelihood> f0 = {TFeatureLikelihood(10, 100), TFeatureLikelihood(20, 0)};
	EXPECT_ANY_THROW(run.add_likelihood_vector(f0));
}

TEST(TAssignmentUpdate_Tests, add_likelihood_vector_rejects_duplicate_molecules) {
	TMassSpecRun run;
	std::vector<TFeatureLikelihood> f0 = {TFeatureLikelihood(10, 100), TFeatureLikelihood(10, 90)};
	EXPECT_ANY_THROW(run.add_likelihood_vector(f0));
}

TEST(TAssignmentUpdate_Tests, initialize_to_unknown_rejects_zero_probability) {
	TMassSpecRun run;
	std::vector<TFeatureLikelihood> f0 = {TFeatureLikelihood(10, 100)};
	run.add_likelihood_vector(f0);
	run.set_probabilities_of_unknowns({0});
	EXPECT_ANY_THROW(run.initialize_assignments_to_unknown());
}
