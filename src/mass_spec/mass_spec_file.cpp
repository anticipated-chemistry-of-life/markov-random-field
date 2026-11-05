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

const std::vector<std::string> HEADER = {
    "species", "extract", "feature", "molecule", "log_correct_density", "log_incorrect_density"};

/// A log density is any finite number: the densities are on the log scale, so there is no sign
/// to enforce. A missing value, which the export writes for a structure its model never
/// calibrated, fails the same way.
[[nodiscard]] double parse_log_density(std::string_view text, std::string_view column,
                                       const std::string &filename, size_t line) {
	double value         = 0.0;
	const auto end       = text.data() + text.size();
	const auto [ptr, ec] = std::from_chars(text.data(), end, value);
	if (ec != std::errc{} || ptr != end || !std::isfinite(value)) {
		throw coretools::TUserError("Line ", line, " of the mass-spec file '", filename,
		                            "' has the ", column, " '", text,
		                            "'. A log density must be a "
		                            "finite number.");
	}
	return value;
}

/// The rows of one feature, as read.
struct TFeatureRows {
	/// Every row of the feature, molecule leaf or not. The unknown weight scales with it.
	size_t number_of_rows = 0;
	/// (molecule leaf, log likelihood ratio) of every kept candidate.
	std::vector<std::pair<uint32_t, double>> candidates;
	/// The log likelihood ratio of every outside structure of the feature.
	std::vector<double> outside_log_ratios;
	/// The molecule name of each entry of `outside_log_ratios`.
	std::vector<std::string> outside_molecules;
};

/// The rows of one run, its features in the order the file first names them.
struct TRunRows {
	size_t species = 0;
	std::string extract;
	std::unordered_map<std::string, size_t> feature_of_id;
	/// The id of each entry of `features`.
	std::vector<std::string> feature_ids;
	std::vector<TFeatureRows> features;
};

[[nodiscard]] std::unordered_map<std::string, uint32_t>
index_of_name(const std::vector<std::string> &names) {
	std::unordered_map<std::string, uint32_t> index;
	index.reserve(names.size());
	for (size_t i = 0; i < names.size(); ++i) { index.emplace(names[i], static_cast<uint32_t>(i)); }
	return index;
}

/// The last bin of the grid. A weight further below its feature's largest than the grid reaches
/// takes it.
constexpr long LAST_BIN = RELATIVE_CALIBRATION_SCORE.size() - 1;

struct TBin {
	uint8_t bin;
	/// The weight was further below the largest than the grid reaches, so `bin` is the last bin.
	bool clamped;
};

/// The bin of a weight `nats_below` (at least 0) under its feature's largest weight:
/// `round(16 nats_below)`, clamped to the last bin. Callers pass a difference of logs rather than
/// the log of a quotient, because the quotient of two valid weights can overflow or underflow. The
/// difference is capped before rounding, so that one as large as a double allows cannot overflow
/// `lround`.
[[nodiscard]] TBin bin_of(double nats_below) {
	constexpr double FLOOR_NATS = (LAST_BIN + 1) / CALIBRATION_BINS_PER_NAT;
	const long bin = std::lround(CALIBRATION_BINS_PER_NAT * std::min(nats_below, FLOOR_NATS));
	if (bin > LAST_BIN) { return {static_cast<uint8_t>(LAST_BIN), true}; }
	return {static_cast<uint8_t>(bin), false};
}

/// log(exp(a) + exp(b)), without exp of a number above 0.
[[nodiscard]] double log_add_exp(double a, double b) {
	const double larger = std::max(a, b);
	return larger + std::log1p(std::exp(std::min(a, b) - larger));
}

/// log of the sum of exp(x) over `xs`, which must not be empty.
[[nodiscard]] double log_sum_exp(const std::vector<double> &xs) {
	const double largest = *std::max_element(xs.begin(), xs.end());
	double sum           = 0.0;
	for (const double x : xs) { sum += std::exp(x - largest); }
	return largest + std::log(sum);
}

/// Moves the outside structure names of `rows` into `labels`.
[[nodiscard]] TMassSpecRun build_run(TRunRows &rows, double log_unknown_rate,
                                     double log_outside_rate, TLoadReport &report,
                                     TRunLabels &labels) {
	TMassSpecRun run;
	run.set_filter_index(0);
	std::vector<uint8_t> unknown_bins;
	for (size_t f = 0; f < rows.features.size(); ++f) {
		auto &feature = rows.features[f];
		if (feature.candidates.empty()) {
			++report.features_not_loaded;
			continue;
		}
		labels.feature_ids.push_back(rows.feature_ids[f]);
		// The unknown weight is rho n_f + tau sum LR_outside. n_f counts every row of the
		// feature. The sum runs in log space, so a large ratio cannot overflow. A tau of 0 has
		// the log -inf, and log_add_exp(x, -inf) is exactly x.
		TFeatureOutside outside;
		outside.log_rho_n_f =
		    log_unknown_rate + std::log(static_cast<double>(feature.number_of_rows));
		outside.molecules = std::move(feature.outside_molecules);
		outside.log_weights.reserve(feature.outside_log_ratios.size());
		for (const double log_ratio : feature.outside_log_ratios) {
			outside.log_weights.push_back(log_outside_rate + log_ratio);
		}
		double log_unknown = outside.log_rho_n_f;
		if (!feature.outside_log_ratios.empty()) {
			log_unknown = log_add_exp(log_unknown,
			                          log_outside_rate + log_sum_exp(feature.outside_log_ratios));
		}
		labels.outside.push_back(std::move(outside));
		report.outside_structures += feature.outside_log_ratios.size();
		double log_largest = log_unknown;
		for (const auto &[molecule, log_ratio] : feature.candidates) {
			log_largest = std::max(log_largest, log_ratio);
		}

		std::vector<TFeatureLikelihood> likelihoods;
		likelihoods.reserve(feature.candidates.size());
		for (const auto &[molecule, log_ratio] : feature.candidates) {
			const auto [bin, clamped] = bin_of(log_largest - log_ratio);
			report.clamped_candidates += clamped;
			likelihoods.emplace_back(molecule, bin);
		}
		run.add_likelihood_vector(likelihoods);
		const auto [unknown_bin, unknown_clamped] = bin_of(log_largest - log_unknown);
		report.clamped_unknowns += unknown_clamped;
		unknown_bins.push_back(unknown_bin);
	}
	run.set_probabilities_of_unknowns(unknown_bins);
	run.initialize_assignments_to_unknown();
	return run;
}

} // namespace

TLoaded load(const std::string &filename, const std::vector<std::string> &species_leaves,
             const std::vector<std::string> &molecule_leaves, double unknown_rate,
             double outside_rate) {
	if (!(unknown_rate > 0.0 && std::isfinite(unknown_rate))) {
		throw coretools::TUserError("The unknown rate (--ms_unknown_rate) must be positive and "
		                            "finite, but got ",
		                            unknown_rate, ".");
	}
	if (!(outside_rate >= 0.0 && std::isfinite(outside_rate))) {
		throw coretools::TUserError("The outside rate (--ms_outside_rate) must be at least 0 and "
		                            "finite, but got ",
		                            outside_rate, ".");
	}
	for (const auto &molecule : molecule_leaves) {
		if (molecule == UNKNOWN_MOLECULE_NAME || molecule == UNLISTED_MOLECULE_NAME) {
			throw coretools::TUserError("The molecule tree has a leaf named '", molecule,
			                            "', which is a name reserved for the unknown molecule in "
			                            "the assignment posterior. Rename the leaf.");
		}
	}
	const auto species_index      = index_of_name(species_leaves);
	const auto molecule_index     = index_of_name(molecule_leaves);
	const double log_unknown_rate = std::log(unknown_rate);
	const double log_outside_rate = std::log(outside_rate);

	coretools::TInputFile file(filename, coretools::FileType::Header, "\t");
	if (file.header() != HEADER) {
		throw coretools::TUserError("The mass-spec file '", filename,
		                            "' must have the tab-separated header 'species extract feature "
		                            "molecule log_correct_density log_incorrect_density'.");
	}

	std::map<std::pair<size_t, std::string>, size_t> run_of_key;
	std::vector<TRunRows> runs;
	std::unordered_set<std::string> seen_rows;
	std::unordered_set<std::string> skipped_species;
	std::unordered_set<std::string> skipped_molecules;
	for (; !file.empty(); file.popFront()) {
		const std::string species_name(file.get(0));
		const std::string molecule_name(file.get(3));
		const double log_correct =
		    parse_log_density(file.get(4), HEADER[4], filename, file.curLine());
		const double log_incorrect =
		    parse_log_density(file.get(5), HEADER[5], filename, file.curLine());
		// two finite densities can still differ by more than a double holds
		const double log_ratio = log_correct - log_incorrect;
		if (!std::isfinite(log_ratio)) {
			throw coretools::TUserError("Line ", file.curLine(), " of the mass-spec file '",
			                            filename,
			                            "' has a likelihood ratio that is not a finite "
			                            "number on the log scale.");
		}
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
		if (new_run) {
			runs.push_back(TRunRows{species_it->second, std::string(file.get(1)), {}, {}, {}});
		}
		auto &run = runs[run_it->second];

		const auto [feature_it, new_feature] =
		    run.feature_of_id.try_emplace(std::string(file.get(2)), run.features.size());
		if (new_feature) {
			run.features.emplace_back();
			run.feature_ids.emplace_back(file.get(2));
		}
		auto &feature = run.features[feature_it->second];

		// Every proposed structure counts toward n_f. An outside structure is no candidate, but
		// its likelihood ratio enters the unknown weight.
		++feature.number_of_rows;
		const auto molecule_it = molecule_index.find(molecule_name);
		if (molecule_it == molecule_index.end()) {
			skipped_molecules.insert(molecule_name);
			feature.outside_log_ratios.push_back(log_ratio);
			feature.outside_molecules.push_back(molecule_name);
			continue;
		}
		feature.candidates.emplace_back(molecule_it->second, log_ratio);
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
		loaded.labels.push_back();
		for (; next_run != runs.end() && next_run->species == species; ++next_run) {
			TRunLabels labels{next_run->extract, {}, {}};
			auto run =
			    build_run(*next_run, log_unknown_rate, log_outside_rate, loaded.report, labels);
			// a run whose every feature went unloaded observes nothing
			if (!run.empty()) {
				loaded.runs_per_species.pushToLast(std::move(run));
				loaded.labels.pushToLast(std::move(labels));
			}
		}
	}
	return loaded;
}

} // namespace mass_spec_file
