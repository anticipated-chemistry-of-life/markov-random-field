#include "mass_spec/assignment_posterior.h"
#include "coretools/Files/TOutputFile.h"
#include <charconv>
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
	const std::string unknown(mass_spec_file::UNKNOWN_MOLECULE_NAME);
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
				const auto row        = [&](const std::string &molecule, uint32_t count) {
					file.writeln(
					    species_leaves[species], label.extract, label.feature_ids[f], molecule,
					    fraction_text(static_cast<double>(count) / static_cast<double>(total)));
				};
				for (size_t c = 0; c < candidates.size(); ++c) {
					row(molecule_leaves[candidates[c].get_molecule_index()], counts[offset + c]);
				}
				row(unknown, counts[offset + candidates.size()]);
				offset += candidates.size() + 1;
			}
		}
	}
}
