#include "mass_spec/assignment_posterior.h"
#include "coretools/Files/TOutputFile.h"
#include <algorithm>
#include <charconv>
#include <cmath>
#include <string>

namespace {

/// The shortest text that reads back as exactly `value`, so that a feature's fractions sum to 1
/// as closely as the counters do.
[[nodiscard]] std::string fraction_text(double value) {
	char buffer[32];
	const auto [end, ec] = std::to_chars(buffer, buffer + sizeof buffer, value);
	if (ec != std::errc{}) {
		throw coretools::TDevError("Could not write the fraction ", value, " as text.");
	}
	return std::string(buffer, end);
}

} // namespace

void write_assignment_posterior(const std::string &filename, const TAssignmentPosterior &posterior,
                                const coretools::TNestedVector<TMassSpecRun> &runs,
                                const coretools::TNestedVector<mass_spec_file::TRunLabels> &labels,
                                const std::vector<std::string> &species_leaves,
                                const std::vector<std::string> &molecule_leaves) {
	const size_t total = posterior.get_total_counts();
	if (total == 0) {
		throw coretools::TDevError("The assignment posterior counted no iteration, so it has no "
		                           "fractions to write.");
	}
	const std::string unlisted(mass_spec_file::UNLISTED_MOLECULE_NAME);
	const std::vector<std::string> header = {"species", "extract", "feature", "molecule",
	                                         "fraction"};
	coretools::TOutputFile file(filename, header, "\t");

	for (size_t species = 0; species < runs.size(); ++species) {
		const auto species_runs = runs.at(species);
		const auto run_labels   = labels.at(species);
		for (size_t r = 0; r < species_runs.size(); ++r) {
			const auto &run    = species_runs[r];
			const auto &label  = run_labels[r];
			const auto &counts = posterior.counts_of_run(species, r);
			size_t offset      = 0;
			for (size_t f = 0; f < run.number_of_features(); ++f) {
				const auto candidates = run.get_likelihoods_for_feature(f);
				const auto row        = [&](const std::string &molecule, double fraction) {
					file.writeln(species_leaves[species], label.extract, label.feature_ids[f],
					             molecule, fraction_text(fraction));
				};
				const auto fraction_of = [&](uint32_t count) {
					return static_cast<double>(count) / static_cast<double>(total);
				};
				for (size_t c = 0; c < candidates.size(); ++c) {
					row(molecule_leaves[candidates[c].get_molecule_index()],
					    fraction_of(counts[offset + c]));
				}

				// P(U) is the unknown molecule's fraction. The split of P(U) follows the exact
				// weights, not the bins. Each share is exp(log weight - log W_U), so a weight
				// too large for a double still divides.
				const auto &outside = label.outside[f];
				double log_largest  = outside.log_rho_n_f;
				for (const double w : outside.log_weights) {
					log_largest = std::max(log_largest, w);
				}
				double sum = std::exp(outside.log_rho_n_f - log_largest);
				for (const double w : outside.log_weights) { sum += std::exp(w - log_largest); }
				const double log_unknown_weight = log_largest + std::log(sum);
				const double p_unknown          = fraction_of(counts[offset + candidates.size()]);
				for (size_t j = 0; j < outside.molecules.size(); ++j) {
					row(outside.molecules[j],
					    p_unknown * std::exp(outside.log_weights[j] - log_unknown_weight));
				}
				row(unlisted, p_unknown * std::exp(outside.log_rho_n_f - log_unknown_weight));
				offset += candidates.size() + 1;
			}
		}
	}
}
