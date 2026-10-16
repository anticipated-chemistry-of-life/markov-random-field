#include "mass_spec/mass_spec_file.h"
#include "constants.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Main/TError.h"
#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <map>
#include <string>
#include <string_view>
#include <system_error>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace mass_spec_file {
namespace {

const std::vector<std::string> HEADER = {"species", "extract", "feature", "molecule",
                                         "calibration_score"};

/// A calibration score is a density, so it must be a positive finite number.
[[nodiscard]] double parse_score(std::string_view text, const std::string &filename, size_t line) {
	double score         = 0.0;
	const auto end       = text.data() + text.size();
	const auto [ptr, ec] = std::from_chars(text.data(), end, score);
	if (ec != std::errc{} || ptr != end || !std::isfinite(score) || score <= 0.0) {
		throw coretools::TUserError("Line ", line, " of the mass-spec file '", filename,
		                            "' has the calibration score '", text,
		                            "'. A calibration score must be a positive finite number.");
	}
	return score;
}

/// The rows of one feature, as read.
struct TFeatureRows {
	/// The best calibration score over all of the feature's rows, candidates or not.
	double reference = 0.0;
	/// (molecule leaf, calibration score) of every kept candidate.
	std::vector<std::pair<uint32_t, double>> candidates;
};

/// The rows of one run, its features in the order the file first names them.
struct TRunRows {
	size_t species = 0;
	std::unordered_map<std::string, size_t> feature_of_id;
	std::vector<TFeatureRows> features;
};

[[nodiscard]] std::unordered_map<std::string, uint32_t>
index_of_name(const std::vector<std::string> &names) {
	std::unordered_map<std::string, uint32_t> index;
	index.reserve(names.size());
	for (size_t i = 0; i < names.size(); ++i) { index.emplace(names[i], static_cast<uint32_t>(i)); }
	return index;
}

/// The last bin of the relative-calibration-score grid. A score below its floor takes it.
constexpr long LAST_BIN = RELATIVE_CALIBRATION_SCORE.size() - 1;

/// The bin of a relative calibration score whose log is `log_ratio`, at most 0:
/// `round(-16 log_ratio)`, unclamped. Callers pass a difference of logs rather than the log of a
/// quotient, because the quotient of two valid scores can underflow to 0.
[[nodiscard]] long bin_of(double log_ratio) {
	return std::lround(-CALIBRATION_BINS_PER_NAT * log_ratio);
}

[[nodiscard]] TMassSpecRun build_run(const TRunRows &rows, uint8_t unknown_bin,
                                     TLoadReport &report) {
	TMassSpecRun run;
	run.set_filter_index(0);
	for (const auto &feature : rows.features) {
		if (feature.candidates.empty()) {
			++report.features_not_loaded;
			continue;
		}
		std::vector<TFeatureLikelihood> likelihoods;
		likelihoods.reserve(feature.candidates.size());
		for (const auto &[molecule, score] : feature.candidates) {
			long bin = bin_of(std::log(score) - std::log(feature.reference));
			if (bin > LAST_BIN) {
				bin = LAST_BIN;
				++report.clamped_candidates;
			}
			likelihoods.emplace_back(molecule, static_cast<uint8_t>(bin));
		}
		run.add_likelihood_vector(likelihoods);
	}
	run.set_probabilities_of_unknowns(std::vector<uint8_t>(run.number_of_features(), unknown_bin));
	run.initialize_assignments_to_unknown();
	return run;
}

} // namespace

TLoaded load(const std::string &filename, const std::vector<std::string> &species_leaves,
             const std::vector<std::string> &molecule_leaves, double unknown_fraction) {
	if (!(unknown_fraction > 0.0 && unknown_fraction <= 1.0)) {
		throw coretools::TUserError("The unknown fraction (--ms_unknown_fraction) must be in "
		                            "(0, 1], but got ",
		                            unknown_fraction, ".");
	}
	const auto species_index  = index_of_name(species_leaves);
	const auto molecule_index = index_of_name(molecule_leaves);
	// an unknown fraction below the grid's floor takes the last bin, as a candidate would
	const auto unknown_bin =
	    static_cast<uint8_t>(std::min(bin_of(std::log(unknown_fraction)), LAST_BIN));

	coretools::TInputFile file(filename, coretools::FileType::Header, "\t");
	if (file.header() != HEADER) {
		throw coretools::TUserError("The mass-spec file '", filename,
		                            "' must have the tab-separated header 'species extract feature "
		                            "molecule calibration_score'.");
	}

	std::map<std::pair<size_t, std::string>, size_t> run_of_key;
	std::vector<TRunRows> runs;
	std::unordered_set<std::string> seen_rows;
	std::unordered_set<std::string> skipped_species;
	std::unordered_set<std::string> skipped_molecules;
	for (; !file.empty(); file.popFront()) {
		const std::string species_name(file.get(0));
		const std::string molecule_name(file.get(3));
		const double score  = parse_score(file.get(4), filename, file.curLine());
		std::string row_key = species_name;
		for (const size_t column : {1, 2, 3}) {
			row_key.push_back('\t');
			row_key.append(file.get(column));
		}
		if (!seen_rows.insert(std::move(row_key)).second) {
			throw coretools::TUserError("Line ", file.curLine(), " of the mass-spec file '",
			                            filename,
			                            "' repeats the (species, extract, feature, molecule) of "
			                            "an earlier row.");
		}

		const auto species_it = species_index.find(species_name);
		if (species_it == species_index.end()) {
			skipped_species.insert(species_name);
			continue;
		}
		const auto [run_it, new_run] =
		    run_of_key.try_emplace({species_it->second, std::string(file.get(1))}, runs.size());
		if (new_run) { runs.push_back(TRunRows{species_it->second, {}, {}}); }
		auto &run = runs[run_it->second];

		const auto [feature_it, new_feature] =
		    run.feature_of_id.try_emplace(std::string(file.get(2)), run.features.size());
		if (new_feature) { run.features.emplace_back(); }
		auto &feature = run.features[feature_it->second];

		// a skipped molecule is no candidate, but its score still sets the reference
		feature.reference      = std::max(feature.reference, score);
		const auto molecule_it = molecule_index.find(molecule_name);
		if (molecule_it == molecule_index.end()) {
			skipped_molecules.insert(molecule_name);
			continue;
		}
		feature.candidates.emplace_back(molecule_it->second, score);
	}

	TLoaded loaded;
	loaded.report.skipped_species   = skipped_species.size();
	loaded.report.skipped_molecules = skipped_molecules.size();
	// a stable sort keeps the runs of one species in the order the file first names them
	std::stable_sort(runs.begin(), runs.end(),
	                 [](const auto &a, const auto &b) { return a.species < b.species; });
	auto next_run = runs.begin();
	for (size_t species = 0; species < species_leaves.size(); ++species) {
		loaded.runs_per_species.push_back();
		for (; next_run != runs.end() && next_run->species == species; ++next_run) {
			auto run = build_run(*next_run, unknown_bin, loaded.report);
			// a run whose every feature went unloaded observes nothing
			if (!run.empty()) { loaded.runs_per_species.pushToLast(std::move(run)); }
		}
	}
	return loaded;
}

} // namespace mass_spec_file
