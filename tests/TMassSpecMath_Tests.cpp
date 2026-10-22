#include "coretools/Containers/TNestedVector.h"
#include "coretools/Main/TRandomGenerator.h"
#include "mass_spec/TMassSpecMath.h"
#include "mass_spec/feature_likelihood.h"
#include "mass_spec/msms_run.h"
#include "storages/storage_backend.h"
#include "gtest/gtest.h"
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <set>
#include <utility>
#include <vector>

namespace {

/// Every molecule index used by the runs below is smaller than this.
constexpr uint32_t NUMBER_OF_MOLECULES = 8;
constexpr size_t NUMBER_OF_SPECIES     = 3;

/// A filter probability that differs between filters and between molecules, so a factor that
/// looked up the wrong (filter, molecule) pair would show.
double filter_probability(size_t filter, uint32_t molecule) {
	return 0.15 + 0.07 * static_cast<double>(molecule) + 0.1 * static_cast<double>(filter);
}

constexpr double CONTAMINATION = 0.03;

using Cell = std::pair<size_t, uint32_t>; // (species leaf, molecule leaf)

/// A {species, molecule} field with exactly `ones` set. Spelled `TFieldStorage`, so the test runs
/// against whichever backend the suite was built with.
TFieldStorage make_field(const std::set<Cell> &ones) {
	TFieldStorage Y(1000, {NUMBER_OF_SPECIES, NUMBER_OF_MOLECULES});
	for (const auto &[species, molecule] : ones) {
		Y.insert_one(Y.get_linear_index_in_container_space({species, molecule}));
	}
	return Y;
}

/// Accepts every proposal, so each update moves the assignments as far as the kernel allows.
class TAcceptAll final : public TAssignmentScorer {
public:
	[[nodiscard]] double log_likelihood_ratio(const TAssignmentProposal &) const override {
		return 1e9;
	}
};

/// Builds a run from one candidate list per feature, every feature starting at the unknown
/// molecule.
TMassSpecRun make_run(std::vector<std::vector<TFeatureLikelihood>> features) {
	TMassSpecRun run;
	for (auto &feature : features) { run.add_likelihood_vector(feature); }
	run.set_probabilities_of_unknowns(std::vector<uint8_t>(run.number_of_features(), 50));
	run.initialize_assignments_to_unknown();
	return run;
}

/// The definition the detected query is a shortcut for: some feature is assigned to `molecule`.
bool detected_by_scan(const TMassSpecRun &run, uint32_t molecule) {
	for (size_t f = 0; f < run.number_of_assignments(); ++f) {
		const auto &a = run.get_current_assignment(f);
		if (!a.is_unknown_molecule() && a.get_molecule_index() == molecule) { return true; }
	}
	return false;
}

/// The runs of the hand-built field:
///   species 0, run A (filter 0): molecules 1 and 4 detected, a third feature at the unknown one
///   species 0, run B (filter 1): molecule 3 detected
///   species 1, one run (filter 0): molecules 0 and 5 detected
///   species 2: no runs
/// so species 0 detects molecule 1 in one run and not in the other, and vice versa for 3.
coretools::TNestedVector<TMassSpecRun> make_runs_per_species() {
	auto run_a = make_run({{TFeatureLikelihood(1, 100), TFeatureLikelihood(2, 120)},
	                       {TFeatureLikelihood(2, 80), TFeatureLikelihood(4, 70)},
	                       {TFeatureLikelihood(6, 60)}});
	run_a.set_current_assignment(0, TFeatureLikelihood(1, 100));
	run_a.set_current_assignment(1, TFeatureLikelihood(4, 70));

	auto run_b = make_run({{TFeatureLikelihood(1, 90), TFeatureLikelihood(3, 110)}});
	run_b.set_filter_index(1);
	run_b.set_current_assignment(0, TFeatureLikelihood(3, 110));

	auto run_c = make_run(
	    {{TFeatureLikelihood(0, 100), TFeatureLikelihood(1, 50)}, {TFeatureLikelihood(5, 30)}});
	run_c.set_current_assignment(0, TFeatureLikelihood(0, 100));
	run_c.set_current_assignment(1, TFeatureLikelihood(5, 30));

	coretools::TNestedVector<TMassSpecRun> runs_per_species;
	runs_per_species.push_back({run_a, run_b});
	runs_per_species.push_back({run_c});
	runs_per_species.push_back();
	return runs_per_species;
}

} // namespace

TEST(TMassSpecMath_Tests, detected_query_agrees_with_a_scan_after_accepted_moves) {
	coretools::instances::randomGenerator().setSeed(7, true);
	// overlapping candidates, so every move type, swaps included, is reachable
	auto run = make_run(
	    {{TFeatureLikelihood(1, 100), TFeatureLikelihood(2, 120)},
	     {TFeatureLikelihood(1, 90), TFeatureLikelihood(2, 80), TFeatureLikelihood(3, 70)},
	     {TFeatureLikelihood(2, 60), TFeatureLikelihood(3, 50), TFeatureLikelihood(5, 40)}});
	run.set_current_assignment(0, TFeatureLikelihood(1, 100));

	const TAcceptAll scorer;
	size_t updates_with_a_detection = 0;
	for (int update = 0; update < 300; ++update) {
		run.update_all_assignments(scorer, 0.3);
		bool any = false;
		for (uint32_t m = 0; m < NUMBER_OF_MOLECULES; ++m) {
			ASSERT_EQ(run.is_detected(m), detected_by_scan(run, m))
			    << "molecule " << m << " after update " << update;
			any = any || detected_by_scan(run, m);
		}
		updates_with_a_detection += any;
	}
	// the chain must actually have visited states with detections for the check to mean anything
	EXPECT_GT(updates_with_a_detection, 0u);
}

TEST(TMassSpecMath_Tests, species_without_runs_gets_a_neutral_factor) {
	const auto runs_per_species = make_runs_per_species();
	for (uint32_t m = 0; m < NUMBER_OF_MOLECULES; ++m) {
		const auto factor = mass_spec_math::per_cell_factor(runs_per_species.at(2), m,
		                                                    filter_probability, CONTAMINATION);
		EXPECT_EQ(factor[0], 1.0) << "molecule " << m;
		EXPECT_EQ(factor[1], 1.0) << "molecule " << m;
	}
}

TEST(TMassSpecMath_Tests, per_cell_factor_is_the_change_in_the_whole_log_likelihood) {
	const auto runs_per_species = make_runs_per_species();
	// present cells both detected and not, in runs of both filters, and in the species without
	// runs
	const std::set<Cell> field  = {{0, 1}, {0, 2}, {1, 5}, {1, 7}, {2, 3}};

	for (size_t species = 0; species < NUMBER_OF_SPECIES; ++species) {
		for (uint32_t molecule = 0; molecule < NUMBER_OF_MOLECULES; ++molecule) {
			auto absent  = field;
			auto present = field;
			absent.erase({species, molecule});
			present.insert({species, molecule});
			const double change =
			    mass_spec_math::log_likelihood(runs_per_species, make_field(present),
			                                   filter_probability, CONTAMINATION) -
			    mass_spec_math::log_likelihood(runs_per_species, make_field(absent),
			                                   filter_probability, CONTAMINATION);

			const auto factor = mass_spec_math::per_cell_factor(
			    runs_per_species.at(species), molecule, filter_probability, CONTAMINATION);
			EXPECT_NEAR(std::log(factor[1]) - std::log(factor[0]), change, 1e-12)
			    << "cell (" << species << ", " << molecule << ")";
		}
	}
}

// A worked example, so that the terms that cancel from every flip above -- the feature terms and
// the bulk of absent, undetected molecules -- are checked too: the contamination move reads them.
TEST(TMassSpecMath_Tests, whole_log_likelihood_of_a_worked_example) {
	auto run = make_run({{TFeatureLikelihood(0, 100)}, {TFeatureLikelihood(1, 120)}});
	run.set_current_assignment(0, TFeatureLikelihood(0, 100)); // feature 1 stays unknown, bin 50
	coretools::TNestedVector<TMassSpecRun> runs_per_species;
	runs_per_species.push_back({run});
	runs_per_species.push_back();
	runs_per_species.push_back();
	const auto flat_filter = [](size_t, uint32_t) { return 0.4; };

	// species 0: molecule 0 detected but absent (contamination), molecule 2 present but not
	// detected (1 - filter), the other six absent and not detected (1 - contamination). Species 1
	// has no runs, so its present cell adds nothing.
	const double expected = std::log(RELATIVE_CALIBRATION_SCORE[100]) +
	                        std::log(RELATIVE_CALIBRATION_SCORE[50]) + std::log(0.1) +
	                        std::log(0.6) + 6.0 * std::log(0.9);
	EXPECT_NEAR(mass_spec_math::log_likelihood(runs_per_species, make_field({{0, 2}, {1, 4}}),
	                                           flat_filter, 0.1),
	            expected, 1e-12);
}
