//
// Reads the mass-spec file into runs (CONTEXT.md, "Mass spectrometry").
//
// The file is tab-separated. Its header is exactly
//
//     species  extract  feature  molecule  calibration_score
//
// and every row names one candidate of one feature. A run is a (species leaf, extract) pair, and
// every run gets the single filter. A feature is unique within its run, not across runs.
//
// A calibration score enters only relative to its feature's best row, as a bin on a 1/16-nat grid.
// See ADR-0012.
//
// No USE_MS_DATA guard covers this file. It compiles, and its tests run, in every build
// configuration.
//

#pragma once

#include "coretools/Containers/TNestedVector.h"
#include "mass_spec/msms_run.h"
#include <cstddef>
#include <string>
#include <vector>

namespace mass_spec_file {

/// What the loader left out, for the log.
struct TLoadReport {
	/// Distinct species that are not species leaves. Their rows are skipped whole.
	size_t skipped_species     = 0;
	/// Distinct molecules that are not molecule leaves, among the rows of kept species.
	size_t skipped_molecules   = 0;
	/// Features of kept species with no candidate left.
	size_t features_not_loaded = 0;
	/// Candidates whose bin fell past 255 and was clamped to it.
	size_t clamped_candidates  = 0;
};

struct TLoaded {
	/// The runs of each species leaf, in leaf-index order. Every assignment starts at the unknown
	/// molecule.
	coretools::TNestedVector<TMassSpecRun> runs_per_species;
	TLoadReport report;
};

/// Reads `filename` into runs. `species_leaves` and `molecule_leaves` are the leaf names in
/// leaf-index order. `unknown_fraction` is the unknown fraction `c`, in `(0, 1]`.
///
/// Throws TUserError on a bad header, on a duplicate (species, extract, feature, molecule) row, on
/// a calibration score that is not a positive finite number, and on `c` outside `(0, 1]`.
[[nodiscard]] TLoaded load(const std::string &filename,
                           const std::vector<std::string> &species_leaves,
                           const std::vector<std::string> &molecule_leaves,
                           double unknown_fraction);

} // namespace mass_spec_file
