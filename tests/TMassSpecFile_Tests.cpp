#include "coretools/Main/TError.h"
#include "mass_spec/mass_spec_file.h"
#include "temp_file.h"
#include "gtest/gtest.h"
#include <cmath>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace {

const std::vector<std::string> SPECIES   = {"Ilex", "Quercus"};
const std::vector<std::string> MOLECULES = {"caffeine", "quercetin", "rutin", "theobromine"};

constexpr uint32_t CAFFEINE    = 0;
constexpr uint32_t QUERCETIN   = 1;
constexpr uint32_t RUTIN       = 2;
constexpr uint32_t THEOBROMINE = 3;

constexpr std::string_view HEADER =
    "species\textract\tfeature\tmolecule\tlog_correct_density\tlog_incorrect_density\n";

/// One row of the mass-spec file, from its two log densities as text.
std::string density_row(const std::string &species, const std::string &extract,
                        const std::string &feature, const std::string &molecule,
                        const std::string &log_correct, const std::string &log_incorrect) {
	return species + "\t" + extract + "\t" + feature + "\t" + molecule + "\t" + log_correct + "\t" +
	       log_incorrect + "\n";
}

/// One row whose likelihood ratio is e^`log_ratio`. Both densities are written, at an offset that
/// keeps either from being 0, so that a loader reading only one of them would get it wrong.
std::string row(const std::string &species, const std::string &extract, const std::string &feature,
                const std::string &molecule, double log_ratio = 0.0) {
	constexpr double LOG_INCORRECT = -2.5;
	return density_row(species, extract, feature, molecule,
	                   std::to_string(LOG_INCORRECT + log_ratio), std::to_string(LOG_INCORRECT));
}

/// Loads `rows` under the standard header, with the leaves above. The default unknown rate puts
/// the unknown weight of a feature of up to nine rows below 1, so a row of likelihood ratio 1
/// takes bin 0.
mass_spec_file::TLoaded load(const std::string &rows, double unknown_rate = 0.1) {
	const TTempFile file("ms.tsv", std::string(HEADER) + rows);
	return mass_spec_file::load(file.path(), SPECIES, MOLECULES, unknown_rate);
}

/// The bin of the unknown molecule in one feature.
uint32_t unknown_bin(const TMassSpecRun &run, size_t feature) {
	return run.get_current_assignment(feature).get_binned_likelihood();
}

/// The candidates of one feature, as molecule leaf -> bin.
std::map<uint32_t, uint32_t> candidates(const TMassSpecRun &run, size_t feature) {
	std::map<uint32_t, uint32_t> bins;
	for (const auto &candidate : run.get_likelihoods_for_feature(feature)) {
		bins[candidate.get_molecule_index()] = candidate.get_binned_likelihood();
	}
	return bins;
}

/// coretools' user error is a CTAD class template, so it cannot be named in an exception
/// declaration and EXPECT_THROW cannot be used on it. Catching the err::TError base and asserting
/// it is *not* a dev error is what "the user gets a clear error" means here.
template<typename F> void expect_user_error(F &&f, const std::string &context) {
	try {
		f();
	} catch (coretools::err::TError &e) {
		EXPECT_FALSE(e.isDevError())
		    << context << ": expected a user error, got a dev error: " << e.what();
		return;
	}
	ADD_FAILURE() << context << ": expected a user error, got none";
}

} // namespace

TEST(TMassSpecFile_Tests, a_likelihood_ratio_above_1_needs_no_special_case) {
	// ln LR: caffeine 2 (-1 - -3), quercetin 1 (-4.5 - -5.5). n_f = 2 and rho = 1 make the unknown
	// weight 2, ln 2 = 0.693. Bins against the largest weight, caffeine: 0, 16 and
	// round(16 * (2 - 0.693)) = 21.
	const auto loaded = load(density_row("Ilex", "e1", "f1", "caffeine", "-1", "-3") +
	                             density_row("Ilex", "e1", "f1", "quercetin", "-4.5", "-5.5"),
	                         1.0);

	ASSERT_EQ(loaded.runs_per_species.size(), SPECIES.size());
	ASSERT_EQ(loaded.runs_per_species.at(0).size(), 1);
	const auto &run = loaded.runs_per_species.at(0)[0];
	ASSERT_EQ(run.number_of_features(), 1);
	EXPECT_EQ(candidates(run, 0), (std::map<uint32_t, uint32_t>{{CAFFEINE, 0}, {QUERCETIN, 16}}));
	EXPECT_EQ(unknown_bin(run, 0), 21);
	EXPECT_TRUE(run.get_current_assignment(0).is_unknown_molecule());
	EXPECT_TRUE(loaded.runs_per_species.at(1).empty());
}

TEST(TMassSpecFile_Tests,
     the_unknown_weight_scales_with_every_row_of_the_feature_in_the_tree_or_not) {
	// rho = 1, so the unknown weight is n_f. f1 has one row: unknown weight 1, a tie with the
	// candidate's ratio of 1. f2 has the candidate and three rows outside the molecule tree:
	// unknown weight 4, so the candidate sits round(16 ln 4) = 22 below it.
	const auto loaded =
	    load(row("Ilex", "e1", "f1", "caffeine") + row("Ilex", "e1", "f2", "caffeine") +
		         row("Ilex", "e1", "f2", "aspirin") + row("Ilex", "e1", "f2", "morphine") +
		         row("Ilex", "e1", "f2", "ibuprofen"),
		     1.0);

	const auto &run = loaded.runs_per_species.at(0)[0];
	ASSERT_EQ(run.number_of_features(), 2);
	EXPECT_EQ(candidates(run, 0), (std::map<uint32_t, uint32_t>{{CAFFEINE, 0}}));
	EXPECT_EQ(unknown_bin(run, 0), 0);
	EXPECT_EQ(candidates(run, 1), (std::map<uint32_t, uint32_t>{{CAFFEINE, 22}}));
}

TEST(TMassSpecFile_Tests, the_unknown_weight_is_the_largest_when_it_exceeds_every_candidate) {
	// rho = 2, n_f = 1: the unknown weight is 2, ln 2 = 0.693. The candidate's ratio is e^-1, so
	// it sits round(16 * (0.693 + 1)) = 27 below it, and the unknown molecule is at bin 0.
	const auto loaded = load(row("Ilex", "e1", "f1", "caffeine", -1.0), 2.0);

	const auto &run = loaded.runs_per_species.at(0)[0];
	EXPECT_EQ(unknown_bin(run, 0), 0);
	EXPECT_EQ(candidates(run, 0), (std::map<uint32_t, uint32_t>{{CAFFEINE, 27}}));
}

TEST(TMassSpecFile_Tests, the_unknown_bin_is_per_feature) {
	// rho = 0.1. f1: the candidate's ratio is 1 and the unknown weight 0.1, so the unknown molecule
	// sits round(16 ln 10) = 37 below. f2: the candidate's ratio is 0.01, so the unknown weight
	// is the largest, at 0, and the candidate sits round(16 ln 10) = 37 below it.
	const auto loaded = load(row("Ilex", "e1", "f1", "caffeine") +
	                         row("Ilex", "e1", "f2", "quercetin", std::log(0.01)));

	const auto &run = loaded.runs_per_species.at(0)[0];
	ASSERT_EQ(run.number_of_features(), 2);
	EXPECT_EQ(unknown_bin(run, 0), 37);
	EXPECT_EQ(unknown_bin(run, 1), 0);
	EXPECT_EQ(candidates(run, 1), (std::map<uint32_t, uint32_t>{{QUERCETIN, 37}}));
}

TEST(TMassSpecFile_Tests, a_weight_below_the_floor_takes_bin_255_and_counts_as_clamped) {
	// rho n_f = 1 for four rows, so the unknown weight ties the best candidate at ln 0. The floor
	// is 255/16 = 15.94 nats below it: ln LR -20 and -30 are far below, -15.9 (bin 254) just above.
	const auto loaded = load(
	    row("Ilex", "e1", "f1", "caffeine") + row("Ilex", "e1", "f1", "quercetin", -20.0) +
	        row("Ilex", "e1", "f1", "rutin", -15.9) + row("Ilex", "e1", "f1", "theobromine", -30.0),
	    0.25);

	const auto &run = loaded.runs_per_species.at(0)[0];
	EXPECT_EQ(candidates(run, 0),
	          (std::map<uint32_t, uint32_t>{
	              {CAFFEINE, 0}, {QUERCETIN, 255}, {RUTIN, 254}, {THEOBROMINE, 255}}));
	EXPECT_EQ(unknown_bin(run, 0), 0);
	EXPECT_EQ(loaded.report.clamped_candidates, 2);
}

TEST(TMassSpecFile_Tests, an_unknown_weight_below_the_floor_takes_bin_255_and_counts_as_clamped) {
	// rho = 1e-9, n_f = 1: ln of the unknown weight is -20.7. f1's candidate has ratio 1, so the
	// unknown weight is 20.7 nats below it, past the floor. f2's candidate has ratio e^-19, so the
	// unknown weight is round(16 * 1.7) = 28 below it, inside the grid.
	const auto loaded = load(
	    row("Ilex", "e1", "f1", "caffeine") + row("Ilex", "e1", "f2", "quercetin", -19.0), 1e-9);

	const auto &run = loaded.runs_per_species.at(0)[0];
	EXPECT_EQ(unknown_bin(run, 0), 255);
	EXPECT_EQ(unknown_bin(run, 1), 28);
	EXPECT_EQ(loaded.report.clamped_unknowns, 1);
	EXPECT_EQ(loaded.report.clamped_candidates, 0);
}

TEST(TMassSpecFile_Tests, a_gap_between_two_weights_as_wide_as_a_double_allows_still_clamps) {
	// ln LR of caffeine is 1e308 and of quercetin -1e308. Each is finite, their gap is 2e308.
	const auto loaded = load(density_row("Ilex", "e1", "f1", "caffeine", "5e307", "-5e307") +
	                         density_row("Ilex", "e1", "f1", "quercetin", "-5e307", "5e307"));

	const auto &run = loaded.runs_per_species.at(0)[0];
	EXPECT_EQ(candidates(run, 0), (std::map<uint32_t, uint32_t>{{CAFFEINE, 0}, {QUERCETIN, 255}}));
	EXPECT_EQ(unknown_bin(run, 0), 255);
}

TEST(TMassSpecFile_Tests, two_finite_log_densities_whose_difference_is_not_finite_are_an_error) {
	expect_user_error(
	    [&] { (void)load(density_row("Ilex", "e1", "f1", "caffeine", "1e308", "-1e308")); },
	    "ln LR of 2e308");
}

TEST(TMassSpecFile_Tests, skipped_species_and_molecules_count_once_however_many_rows_name_them) {
	const auto loaded =
	    load(row("Pinus", "e1", "f1", "caffeine") + row("Pinus", "e2", "f1", "caffeine") +
		     row("Abies", "e1", "f1", "caffeine") + row("Ilex", "e1", "f1", "caffeine") +
		     row("Ilex", "e1", "f1", "aspirin") + row("Ilex", "e1", "f2", "aspirin") +
		     row("Ilex", "e1", "f2", "quercetin") + row("Quercus", "e1", "f1", "aspirin") +
		     row("Quercus", "e1", "f1", "morphine") + row("Quercus", "e1", "f1", "rutin"));

	EXPECT_EQ(loaded.report.skipped_species, 2);
	EXPECT_EQ(loaded.report.skipped_molecules, 2);
	EXPECT_EQ(loaded.report.features_not_loaded, 0);
	const auto &ilex = loaded.runs_per_species.at(0)[0];
	EXPECT_EQ(candidates(ilex, 0), (std::map<uint32_t, uint32_t>{{CAFFEINE, 0}}));
	EXPECT_EQ(candidates(ilex, 1), (std::map<uint32_t, uint32_t>{{QUERCETIN, 0}}));
}

TEST(TMassSpecFile_Tests, a_feature_with_only_skipped_candidates_is_not_loaded_and_is_counted) {
	const auto loaded =
	    load(row("Ilex", "e1", "f1", "aspirin") + row("Ilex", "e1", "f1", "morphine") +
		     row("Ilex", "e1", "f2", "caffeine") + row("Quercus", "e1", "f1", "aspirin"));

	EXPECT_EQ(loaded.report.features_not_loaded, 2);
	const auto &ilex = loaded.runs_per_species.at(0);
	ASSERT_EQ(ilex.size(), 1);
	ASSERT_EQ(ilex[0].number_of_features(), 1);
	EXPECT_EQ(candidates(ilex[0], 0), (std::map<uint32_t, uint32_t>{{CAFFEINE, 0}}));
}

TEST(TMassSpecFile_Tests, two_extracts_are_two_runs_and_a_feature_id_is_one_feature_per_run) {
	const auto loaded =
	    load(row("Ilex", "e1", "f1", "caffeine") + row("Quercus", "e1", "f1", "rutin") +
		     row("Ilex", "e2", "f1", "theobromine") + row("Ilex", "e1", "f1", "quercetin"));

	const auto &ilex = loaded.runs_per_species.at(0);
	ASSERT_EQ(ilex.size(), 2);
	ASSERT_EQ(ilex[0].number_of_features(), 1);
	ASSERT_EQ(ilex[1].number_of_features(), 1);
	EXPECT_EQ(candidates(ilex[0], 0),
	          (std::map<uint32_t, uint32_t>{{CAFFEINE, 0}, {QUERCETIN, 0}}));
	EXPECT_EQ(candidates(ilex[1], 0), (std::map<uint32_t, uint32_t>{{THEOBROMINE, 0}}));

	const auto &quercus = loaded.runs_per_species.at(1);
	ASSERT_EQ(quercus.size(), 1);
	EXPECT_EQ(candidates(quercus[0], 0), (std::map<uint32_t, uint32_t>{{RUTIN, 0}}));
}

TEST(TMassSpecFile_Tests, a_header_other_than_the_six_columns_in_order_is_an_error) {
	const std::string body = "Ilex\te1\tf1\tcaffeine\t-1\t-2\n";
	for (const std::string header :
	     {// the old contract
	      "species\textract\tfeature\tmolecule\tcalibration_score\n",
	      "species\textract\tfeature\tmolecule\tlog_incorrect_density\tlog_correct_density\n",
	      "species\tsample\tfeature\tmolecule\tlog_correct_density\tlog_incorrect_density\n",
	      "species\textract\tfeature\tmolecule\tlog_correct_density\n",
	      "species\textract\tfeature\tmolecule\tlog_correct_density\tlog_incorrect_density\tx\n"}) {
		const TTempFile file("ms.tsv", header + body);
		expect_user_error([&] { (void)mass_spec_file::load(file.path(), SPECIES, MOLECULES, 0.1); },
		                  "header " + header);
	}
}

TEST(TMassSpecFile_Tests, a_duplicate_row_is_an_error_even_where_it_would_be_skipped) {
	for (const std::string &duplicated :
	     {row("Ilex", "e1", "f1", "caffeine"), row("Ilex", "e1", "f1", "aspirin"),
	      row("Pinus", "e1", "f1", "caffeine")}) {
		expect_user_error(
		    [&] { (void)load(duplicated + row("Ilex", "e1", "f1", "rutin") + duplicated); },
		    "row " + duplicated);
	}
}

TEST(TMassSpecFile_Tests, a_log_density_that_is_not_a_finite_number_is_an_error) {
	for (const std::string value : {"nan", "inf", "-inf", "abc", "1x", ""}) {
		expect_user_error(
		    [&] {
			    (void)load(row("Ilex", "e1", "f1", "caffeine") +
				           density_row("Ilex", "e1", "f1", "rutin", value, "-1"));
		    },
		    "log_correct_density '" + value + "'");
		expect_user_error(
		    [&] {
			    (void)load(row("Ilex", "e1", "f1", "caffeine") +
				           density_row("Ilex", "e1", "f1", "rutin", "-1", value));
		    },
		    "log_incorrect_density '" + value + "'");
	}
}

TEST(TMassSpecFile_Tests, a_log_density_on_a_skipped_row_must_still_be_a_finite_number) {
	for (const std::string &bad : {density_row("Ilex", "e1", "f1", "aspirin", "nan", "-1"),
	                               density_row("Pinus", "e1", "f1", "caffeine", "-1", "inf")}) {
		expect_user_error([&] { (void)load(row("Ilex", "e1", "f1", "caffeine") + bad); },
		                  "row " + bad);
	}
}

TEST(TMassSpecFile_Tests, an_unknown_rate_that_is_not_positive_and_finite_is_an_error) {
	for (const double unknown_rate : {0.0, -0.1, std::nan(""), HUGE_VAL, -HUGE_VAL}) {
		expect_user_error([&] { (void)load(row("Ilex", "e1", "f1", "caffeine"), unknown_rate); },
		                  "rho = " + std::to_string(unknown_rate));
	}
}
