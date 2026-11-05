//
// How often each feature of each run was explained by each of its candidates and by the unknown
// molecule, and what the unknown molecule's share is made of.
//
// The posterior of Y counts a cell on one iteration in its thinning factor and starts again
// when the burn-in ends. This counts the assignments on the same iterations and clears at the same
// point, so that both report the same chain and their denominators agree.
//
// A feature's slots are its candidates, in the order the run lists them (ascending molecule leaf),
// followed by one slot for the unknown molecule. The writer divides the unknown molecule's
// fraction P(U) with the exact weights the loader kept (ADR-0020): each outside structure j gets
// P(U) tau LR_j / W_U, and the row `<unlisted>` gets P(U) rho n_f / W_U. Every row is written,
// including the ones that were never counted, so a file holds a row for each explanation a
// feature could take and its fractions sum to 1. The chain counts nothing more than the slots:
// the split costs no state.
//

#pragma once

#include "coretools/Containers/TNestedVector.h"
#include "coretools/Main/TError.h"
#include "mass_spec/mass_spec_file.h"
#include "mass_spec/msms_run.h"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

class TAssignmentPosterior {
private:
	/// Counters by species leaf, then run. A run's counters hold, feature after feature, one slot
	/// per candidate and one for the unknown molecule.
	std::vector<std::vector<std::vector<uint32_t>>> _counts;
	/// The number of iterations actually counted.
	size_t _total_counts    = 0;
	size_t _thinning_factor = 1;

public:
	/// @param runs            the runs of every species leaf
	/// @param thinning_factor the field's own, so that both posteriors count the same iterations
	void initialize(const coretools::TNestedVector<TMassSpecRun> &runs, size_t thinning_factor) {
		if (thinning_factor == 0) {
			throw coretools::TDevError("An assignment posterior needs a thinning factor of at "
			                           "least one.");
		}
		_thinning_factor = thinning_factor;
		_total_counts    = 0;
		_counts.assign(runs.size(), {});
		for (size_t species = 0; species < runs.size(); ++species) {
			for (const auto &run : runs.at(species)) {
				size_t slots = 0;
				for (size_t f = 0; f < run.number_of_features(); ++f) {
					slots += run.number_of_candidates(f) + 1;
				}
				_counts[species].emplace_back(slots, 0);
			}
		}
	}

	/// Counts the current assignment of every feature, on one iteration in the thinning factor.
	/// `runs` must be the runs `initialize` was given.
	void add_to_counter(size_t iteration, const coretools::TNestedVector<TMassSpecRun> &runs) {
		if (iteration % _thinning_factor != 0) { return; }
		for (size_t species = 0; species < _counts.size(); ++species) {
			const auto species_runs = runs.at(species);
			for (size_t r = 0; r < species_runs.size(); ++r) {
				const auto &run = species_runs[r];
				auto &counts    = _counts[species][r];
				size_t offset   = 0;
				for (size_t f = 0; f < run.number_of_features(); ++f) {
					const size_t n_candidates = run.number_of_candidates(f);
					const auto &assignment    = run.get_current_assignment(f);
					if (assignment.is_unknown_molecule()) {
						++counts[offset + n_candidates];
					} else {
						const auto found =
						    run.is_molecule_in_feature(f, assignment.get_molecule_index());
						if (!found.found) {
							throw coretools::TDevError("Feature ", f,
							                           " is assigned to a molecule it does not "
							                           "list as a candidate.");
						}
						++counts[offset + found.index];
					}
					offset += n_candidates + 1;
				}
			}
		}
		++_total_counts;
	}

	/// Forgets what was counted, as the burn-in ends.
	void reset_counts() {
		for (auto &species : _counts) {
			for (auto &run : species) { std::fill(run.begin(), run.end(), 0); }
		}
		_total_counts = 0;
	}

	[[nodiscard]] size_t get_total_counts() const { return _total_counts; }
	[[nodiscard]] size_t get_thinning_factor() const { return _thinning_factor; }

	/// The counters of one run, feature after feature: one slot per candidate in the order the run
	/// lists them, then one for the unknown molecule.
	[[nodiscard]] const std::vector<uint32_t> &counts_of_run(size_t species, size_t run) const {
		return _counts.at(species).at(run);
	}
};

/// Writes `filename`: one row per (species, extract, feature, explanation), with the fraction of
/// the counted iterations the feature spent on that explanation. Every candidate, every outside
/// structure and the row `mass_spec_file::UNLISTED_MOLECULE_NAME` have a row, zero fractions
/// included. The unknown molecule has no row of its own.
///
/// @param runs           the runs the posterior counted
/// @param labels         the names of those runs and their features (`mass_spec_file::TLoaded`)
/// @param species_leaves the species leaf names, in leaf-index order
/// @param molecule_leaves the molecule leaf names, in leaf-index order
void write_assignment_posterior(const std::string &filename, const TAssignmentPosterior &posterior,
                                const coretools::TNestedVector<TMassSpecRun> &runs,
                                const coretools::TNestedVector<mass_spec_file::TRunLabels> &labels,
                                const std::vector<std::string> &species_leaves,
                                const std::vector<std::string> &molecule_leaves);
