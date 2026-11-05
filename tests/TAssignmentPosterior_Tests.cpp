#include "mass_spec/assignment_posterior.h"
#include "mass_spec/mass_spec_file.h"
#include "temp_file.h"
#include "gtest/gtest.h"
#include <cmath>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

const std::vector<std::string> SPECIES   = {"Ilex"};
const std::vector<std::string> MOLECULES = {"caffeine", "quercetin"};

constexpr std::string_view HEADER =
    "species\textract\tfeature\tmolecule\tlog_correct_density\tlog_incorrect_density\n";

/// One row whose likelihood ratio is e^`log_ratio`. The density is written in full precision, so
/// that the ratio the loader reads is the one the test names.
std::string row(const std::string &feature, const std::string &molecule, double log_ratio) {
	std::ostringstream log_correct;
	log_correct.precision(17);
	log_correct << -2.5 + log_ratio;
	return "Ilex\te1\t" + feature + "\t" + molecule + "\t" + log_correct.str() + "\t-2.5\n";
}

/// The rows of one feature of the written file, as molecule -> fraction.
using TFractions = std::map<std::string, double>;

/// Parses the written file into feature -> (molecule -> fraction).
std::map<std::string, TFractions> parse(const std::string &content) {
	std::map<std::string, TFractions> features;
	std::istringstream in(content);
	std::string line;
	std::getline(in, line);
	EXPECT_EQ(line, "species\textract\tfeature\tmolecule\tfraction");
	while (std::getline(in, line)) {
		std::istringstream cells(line);
		std::string species, extract, feature, molecule, fraction;
		std::getline(cells, species, '\t');
		std::getline(cells, extract, '\t');
		std::getline(cells, feature, '\t');
		std::getline(cells, molecule, '\t');
		std::getline(cells, fraction, '\t');
		EXPECT_TRUE(features[feature].emplace(molecule, std::stod(fraction)).second)
		    << "feature " << feature << " writes " << molecule << " twice";
	}
	return features;
}

/// Loads `rows`, counts `unknown_iterations` iterations with every feature on the unknown molecule
/// and one with every feature on its first candidate, then writes the posterior and parses it.
std::map<std::string, TFractions> write_and_parse(const std::string &rows, double unknown_rate,
                                                  double outside_rate, size_t unknown_iterations) {
	const TTempFile data("ms.tsv", std::string(HEADER) + rows);
	auto loaded = mass_spec_file::load(data.path(), SPECIES, MOLECULES, unknown_rate, outside_rate);

	TAssignmentPosterior posterior;
	posterior.initialize(loaded.runs_per_species, 1);
	auto &run        = loaded.runs_per_species.at(0)[0];
	size_t iteration = 0;
	for (; iteration < unknown_iterations; ++iteration) {
		posterior.add_to_counter(iteration, loaded.runs_per_species);
	}
	for (size_t f = 0; f < run.number_of_features(); ++f) {
		run.set_current_assignment(f, run.get_likelihoods_for_feature(f)[0]);
	}
	posterior.add_to_counter(iteration, loaded.runs_per_species);

	const TTempFile out("posterior.tsv");
	write_assignment_posterior(out.path(), posterior, loaded.runs_per_species, loaded.labels,
	                           SPECIES, MOLECULES);
	return parse(out.content());
}

double sum(const TFractions &fractions) {
	double total = 0.0;
	for (const auto &[molecule, fraction] : fractions) { total += fraction; }
	return total;
}

} // namespace

TEST(
    TAssignmentPosterior_Tests,
    the_unknown_molecules_fraction_is_divided_between_the_outside_structures_and_the_unlisted_row) {
	// rho = 1 and n_f = 3, so rho n_f is 3. The outside structures have likelihood
	// ratios 4 and 2 and tau = 0.5, so their weights are 2 and 1. The unknown weight is 6. The
	// feature is on the unknown molecule in 3 of 4 counted iterations, P(U) = 3/4, so the rows
	// are 3/4 * 2/6, 3/4 * 1/6 and 3/4 * 3/6. The fourth iteration is on caffeine, 1/4.
	const auto features =
	    write_and_parse(row("f1", "caffeine", 1.0) + row("f1", "aspirin", std::log(4.0)) +
		                    row("f1", "ibuprofen", std::log(2.0)),
		                1.0, 0.5, 3);

	const auto &f1 = features.at("f1");
	EXPECT_EQ(f1.size(), 4);
	EXPECT_DOUBLE_EQ(f1.at("caffeine"), 0.25);
	EXPECT_NEAR(f1.at("aspirin"), 0.25, 1e-12);
	EXPECT_NEAR(f1.at("ibuprofen"), 0.125, 1e-12);
	EXPECT_NEAR(f1.at("<unlisted>"), 0.375, 1e-12);
	EXPECT_NEAR(sum(f1), 1.0, 1e-12);
	EXPECT_EQ(f1.count("<unknown>"), 0);
}

TEST(
    TAssignmentPosterior_Tests,
    a_feature_with_no_outside_structure_puts_the_unknown_molecules_whole_fraction_on_the_unlisted_row) {
	const auto features =
	    write_and_parse(row("f1", "caffeine", 0.0) + row("f2", "quercetin", 0.0), 1.0, 0.5, 1);

	for (const auto &[feature, fractions] : features) {
		EXPECT_EQ(fractions.size(), 2) << feature;
		EXPECT_DOUBLE_EQ(fractions.at("<unlisted>"), 0.5) << feature;
		EXPECT_NEAR(sum(fractions), 1.0, 1e-12) << feature;
	}
}

TEST(TAssignmentPosterior_Tests, an_outside_rate_of_zero_writes_zero_for_every_outside_structure) {
	const auto features =
	    write_and_parse(row("f1", "caffeine", 0.0) + row("f1", "aspirin", 5.0), 1.0, 0.0, 1);

	const auto &f1 = features.at("f1");
	EXPECT_EQ(f1.at("aspirin"), 0.0);
	EXPECT_DOUBLE_EQ(f1.at("<unlisted>"), 0.5);
	EXPECT_NEAR(sum(f1), 1.0, 1e-12);
}

TEST(TAssignmentPosterior_Tests,
     a_dominant_outside_structure_takes_the_unknown_molecules_fraction) {
	// e^800 overflows a double. The split stays on the log scale: tau LR = e^800, rho n_f = 2.
	const auto features =
	    write_and_parse(row("f1", "caffeine", 0.0) + row("f1", "aspirin", 800.0), 1.0, 1.0, 1);

	const auto &f1 = features.at("f1");
	EXPECT_NEAR(f1.at("aspirin"), 0.5, 1e-12);
	EXPECT_NEAR(f1.at("<unlisted>"), 0.0, 1e-12);
	EXPECT_NEAR(sum(f1), 1.0, 1e-12);
}
