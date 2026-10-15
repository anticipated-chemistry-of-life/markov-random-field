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

constexpr std::string_view HEADER = "species\textract\tfeature\tmolecule\tcalibration_score\n";

/// One row of the mass-spec file.
std::string row(const std::string &species, const std::string &extract, const std::string &feature,
                const std::string &molecule, const std::string &score) {
	return species + "\t" + extract + "\t" + feature + "\t" + molecule + "\t" + score + "\n";
}

/// Loads `rows` under the standard header, with the leaves above.
mass_spec_file::TLoaded load(const std::string &rows, double unknown_fraction = 0.1) {
	const TTempFile file("ms.tsv", std::string(HEADER) + rows);
	return mass_spec_file::load(file.path(), SPECIES, MOLECULES, unknown_fraction);
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
		EXPECT_FALSE(e.isDevError()) << context << ": expected a user error, got a dev error: "
		                             << e.what();
		return;
	}
	ADD_FAILURE() << context << ": expected a user error, got none";
}

} // namespace

TEST(TMassSpecFile_Tests, the_reference_takes_bin_0_and_a_ratio_of_e_to_the_minus_1_bin_16) {
	const std::string two_over_e = std::to_string(2.0 * std::exp(-1.0));
	const auto loaded            = load(row("Ilex", "e1", "f1", "caffeine", "2.0") +
	                                    row("Ilex", "e1", "f1", "quercetin", two_over_e));

	ASSERT_EQ(loaded.runs_per_species.size(), SPECIES.size());
	ASSERT_EQ(loaded.runs_per_species.at(0).size(), 1);
	const auto &run = loaded.runs_per_species.at(0)[0];
	ASSERT_EQ(run.number_of_features(), 1);
	EXPECT_EQ(candidates(run, 0), (std::map<uint32_t, uint32_t>{{CAFFEINE, 0}, {QUERCETIN, 16}}));
	EXPECT_TRUE(loaded.runs_per_species.at(1).empty());
}

TEST(TMassSpecFile_Tests, a_ratio_below_the_floor_takes_bin_255_and_counts_as_clamped) {
	// the floor is e^(-255/16), about 1.2e-7; 1e-9 is far below it and 1e-7 just above
	const auto loaded = load(row("Ilex", "e1", "f1", "caffeine", "1") +
	                         row("Ilex", "e1", "f1", "quercetin", "1e-9") +
	                         row("Ilex", "e1", "f1", "rutin", "1e-12") +
	                         row("Ilex", "e1", "f1", "theobromine", "1.3e-7"));

	const auto &run = loaded.runs_per_species.at(0)[0];
	EXPECT_EQ(candidates(run, 0),
	          (std::map<uint32_t, uint32_t>{
	              {CAFFEINE, 0}, {QUERCETIN, 255}, {RUTIN, 255}, {THEOBROMINE, 254}}));
	EXPECT_EQ(loaded.report.clamped_candidates, 2);
}

TEST(TMassSpecFile_Tests, every_feature_starts_at_the_unknown_molecule_with_the_unknown_bin) {
	const std::string rows = row("Ilex", "e1", "f1", "caffeine", "1") +
	                         row("Ilex", "e1", "f2", "quercetin", "1");
	for (const auto &[unknown_fraction, bin] : std::map<double, uint32_t>{{0.1, 37}, {1.0, 0}}) {
		const auto loaded = load(rows, unknown_fraction);
		const auto &run   = loaded.runs_per_species.at(0)[0];
		ASSERT_EQ(run.number_of_features(), 2);
		for (size_t f = 0; f < 2; ++f) {
			EXPECT_TRUE(run.get_current_assignment(f).is_unknown_molecule());
			EXPECT_EQ(run.get_current_assignment(f).get_binned_likelihood(), bin)
			    << "at c = " << unknown_fraction;
		}
	}
}

TEST(TMassSpecFile_Tests, skipped_species_and_molecules_count_once_however_many_rows_name_them) {
	const auto loaded = load(row("Pinus", "e1", "f1", "caffeine", "1") +
	                         row("Pinus", "e2", "f1", "caffeine", "1") +
	                         row("Abies", "e1", "f1", "caffeine", "1") +
	                         row("Ilex", "e1", "f1", "caffeine", "1") +
	                         row("Ilex", "e1", "f1", "aspirin", "0.5") +
	                         row("Ilex", "e1", "f2", "aspirin", "0.5") +
	                         row("Ilex", "e1", "f2", "quercetin", "1") +
	                         row("Quercus", "e1", "f1", "aspirin", "0.5") +
	                         row("Quercus", "e1", "f1", "morphine", "0.5") +
	                         row("Quercus", "e1", "f1", "rutin", "1"));

	EXPECT_EQ(loaded.report.skipped_species, 2);
	EXPECT_EQ(loaded.report.skipped_molecules, 2);
	EXPECT_EQ(loaded.report.features_not_loaded, 0);
	const auto &ilex = loaded.runs_per_species.at(0)[0];
	EXPECT_EQ(candidates(ilex, 0), (std::map<uint32_t, uint32_t>{{CAFFEINE, 0}}));
	EXPECT_EQ(candidates(ilex, 1), (std::map<uint32_t, uint32_t>{{QUERCETIN, 0}}));
}

TEST(TMassSpecFile_Tests, a_skipped_molecule_sets_its_features_reference) {
	const auto loaded = load(row("Ilex", "e1", "f1", "aspirin", std::to_string(std::exp(1.0))) +
	                         row("Ilex", "e1", "f1", "caffeine", "1"));

	EXPECT_EQ(candidates(loaded.runs_per_species.at(0)[0], 0),
	          (std::map<uint32_t, uint32_t>{{CAFFEINE, 16}}));
}

TEST(TMassSpecFile_Tests, a_feature_with_only_skipped_candidates_is_not_loaded_and_is_counted) {
	const auto loaded = load(row("Ilex", "e1", "f1", "aspirin", "1") +
	                         row("Ilex", "e1", "f1", "morphine", "0.5") +
	                         row("Ilex", "e1", "f2", "caffeine", "1") +
	                         row("Quercus", "e1", "f1", "aspirin", "1"));

	EXPECT_EQ(loaded.report.features_not_loaded, 2);
	const auto &ilex = loaded.runs_per_species.at(0);
	ASSERT_EQ(ilex.size(), 1);
	ASSERT_EQ(ilex[0].number_of_features(), 1);
	EXPECT_EQ(candidates(ilex[0], 0), (std::map<uint32_t, uint32_t>{{CAFFEINE, 0}}));
}

TEST(TMassSpecFile_Tests, two_extracts_are_two_runs_and_a_feature_id_is_one_feature_per_run) {
	const auto loaded = load(row("Ilex", "e1", "f1", "caffeine", "1") +
	                         row("Quercus", "e1", "f1", "rutin", "1") +
	                         row("Ilex", "e2", "f1", "theobromine", "1") +
	                         row("Ilex", "e1", "f1", "quercetin", "1"));

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

TEST(TMassSpecFile_Tests, a_header_other_than_the_five_columns_in_order_is_an_error) {
	const std::string body = "Ilex\te1\tf1\tcaffeine\t1\n";
	for (const std::string header :
	     {"species\textract\tfeature\tcalibration_score\tmolecule\n",
	      "species\tsample\tfeature\tmolecule\tcalibration_score\n",
	      "species\textract\tfeature\tmolecule\n",
	      "species\textract\tfeature\tmolecule\tcalibration_score\tx\n"}) {
		const TTempFile file("ms.tsv", header + body);
		expect_user_error(
		    [&] { (void)mass_spec_file::load(file.path(), SPECIES, MOLECULES, 0.1); },
		    "header " + header);
	}
}

TEST(TMassSpecFile_Tests, a_duplicate_row_is_an_error_even_where_it_would_be_skipped) {
	for (const std::string &duplicated : {row("Ilex", "e1", "f1", "caffeine", "0.5"),
	                                      row("Ilex", "e1", "f1", "aspirin", "0.5"),
	                                      row("Pinus", "e1", "f1", "caffeine", "0.5")}) {
		expect_user_error(
		    [&] { (void)load(duplicated + row("Ilex", "e1", "f1", "rutin", "1") + duplicated); },
		    "row " + duplicated);
	}
}

TEST(TMassSpecFile_Tests, a_score_that_is_not_a_positive_finite_number_is_an_error) {
	for (const std::string score : {"0", "-1", "nan", "inf", "abc", "1x"}) {
		expect_user_error(
		    [&] {
			    (void)load(row("Ilex", "e1", "f1", "caffeine", "1") +
			               row("Ilex", "e1", "f1", "rutin", score));
		    },
		    "score '" + score + "'");
	}
}

TEST(TMassSpecFile_Tests, an_unknown_fraction_outside_0_to_1_is_an_error) {
	for (const double unknown_fraction : {0.0, -0.1, 1.0001, std::nan("")}) {
		expect_user_error(
		    [&] { (void)load(row("Ilex", "e1", "f1", "caffeine", "1"), unknown_fraction); },
		    "c = " + std::to_string(unknown_fraction));
	}
}
