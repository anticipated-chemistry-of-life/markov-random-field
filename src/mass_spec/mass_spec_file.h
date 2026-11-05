//
// Reads the mass-spec file into runs (CONTEXT.md, "Mass spectrometry").
//
// The file is tab-separated. Its header is exactly
//
//     species  extract  feature  molecule  log_correct_density  log_incorrect_density
//
// and every row names one structure proposed for one feature. A run is a (species leaf, extract)
// pair, and every run gets the single filter. A feature is unique within its run, not across runs.
//
// A candidate's weight is its likelihood ratio, exp(log correct density - log incorrect density).
// A feature's unknown weight is rho n_f + tau sum LR. Here rho is the unknown rate, and n_f
// counts every proposed structure of the feature, candidate or not. The outside rate is tau. The
// sum runs over the outside structures of the feature: its proposed structures that are not
// molecule leaves. At tau = 0 the unknown weight is rho n_f.
//
// Each weight enters as a bin on a 1/16-nat grid, relative to the largest weight of its feature,
// the unknown weight included. See ADR-0018 and ADR-0020.
//
// No USE_MS_DATA guard covers this file. It compiles, and its tests run, in every build
// configuration.
//

#pragma once

#include "coretools/Containers/TNestedVector.h"
#include "mass_spec/msms_run.h"
#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace mass_spec_file {

/// The name the unknown molecule is written under in the assignment posterior, and the one the
/// simulator's true assignment will use (a later ticket of #99). `<` cannot occur in an InChIKey
/// or a SMILES, and the loader rejects a molecule leaf that carries it, so no molecule leaf can be
/// mistaken for it.
inline constexpr std::string_view UNKNOWN_MOLECULE_NAME = "<unknown>";

/// What the loader left out, for the log.
struct TLoadReport {
	/// Distinct species that are not species leaves. Their rows are skipped whole.
	size_t skipped_species     = 0;
	/// Distinct molecules that are not molecule leaves, among the rows of kept species.
	size_t skipped_molecules   = 0;
	/// Outside structures in the loaded features, counted by row. Each one's likelihood ratio
	/// enters the unknown weight of its feature, times tau.
	size_t outside_structures  = 0;
	/// Features of kept species with no candidate left.
	size_t features_not_loaded = 0;
	/// Candidates whose bin fell past 255 and was clamped to it.
	size_t clamped_candidates  = 0;
	/// Features whose unknown molecule's bin fell past 255 and was clamped to it.
	size_t clamped_unknowns    = 0;
};

/// The names the file gave one loaded run and its features.
struct TRunLabels {
	std::string extract;
	/// One id per loaded feature, in the feature order of the run.
	std::vector<std::string> feature_ids;
};

struct TLoaded {
	/// The runs of each species leaf, in leaf-index order. Every assignment starts at the unknown
	/// molecule, whose bin is per feature.
	coretools::TNestedVector<TMassSpecRun> runs_per_species;
	/// The labels of `runs_per_species`, run for run.
	coretools::TNestedVector<TRunLabels> labels;
	TLoadReport report;
};

/// Reads `filename` into runs. `species_leaves` and `molecule_leaves` are the leaf names in
/// leaf-index order. `unknown_rate` is the unknown rate rho and `outside_rate` the outside rate
/// tau.
///
/// Throws TUserError on a bad header, on a duplicate (species, extract, feature, molecule) row, on
/// a log density that is not a finite number (a missing one included), on a likelihood ratio
/// whose log is not finite, on a rho that is not positive and finite, on a tau that is negative or
/// not finite, and on a molecule leaf named `UNKNOWN_MOLECULE_NAME`.
[[nodiscard]] TLoaded load(const std::string &filename,
                           const std::vector<std::string> &species_leaves,
                           const std::vector<std::string> &molecule_leaves, double unknown_rate,
                           double outside_rate);

} // namespace mass_spec_file
