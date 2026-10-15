//
// The likelihood of the mass-spectrometry data source.
//
// A species' runs observe the field through their assignments (CONTEXT.md, "Mass spectrometry").
// A molecule is detected in a run when some feature is assigned to it. Whether it is detected
// depends on whether it is present. Take the field state `y` of the (species, molecule) cell.
// Take the filter probability `f` of (run's filter, molecule) and the contamination probability
// `c`. One run then contributes
//
//     P(detected | y = 1) = f          P(detected | y = 0) = c
//     P(not detected | y = 1) = 1 - f  P(not detected | y = 0) = 1 - c
//
// which is `TMassSpecRun::probability_of_assignment`. Every run also carries a feature term. It is
// the relative calibration score of the molecule each feature is assigned to (ADR-0012). It does
// not read the field.
//
// Conditional on the assignments, the source is therefore per-cell shaped. A cell's factor is the
// product over the species' runs of its detection term. The whole log-likelihood is the sum of
// every feature term and every detection term. The two agree exactly. A field draw that reads the
// factor and a parameter move that reads the whole log-likelihood therefore score the same model.
//
// Both skip a run without features. Such a run contributes no feature term and no detection term.
//
// Nothing here is guarded by USE_MS_DATA. It compiles, and is unit-tested, in every build
// configuration. The caller passes the filter probability as a callable
// `(size_t filter, uint32_t molecule) -> double`, so nothing here reads a parameter.
//

#pragma once

#include "coretools/Containers/TNestedVector.h"
#include "coretools/Containers/TView.h"
#include "coretools/Main/TError.h"
#include "coretools/Math/TSumLog.h"
#include "mass_spec/feature_likelihood.h"
#include "mass_spec/msms_run.h"
#include "storages/storage_concepts.h"
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace mass_spec_math {

/// The feature term of one feature: the relative calibration score of the molecule it is assigned
/// to. It reads the assignment's bin. The unknown molecule carries its bin in the same place, so it
/// reads the same way.
[[nodiscard]] inline double feature_term(const TFeatureLikelihood &assignment) {
	return RELATIVE_CALIBRATION_SCORE[assignment.get_binned_likelihood()];
}

/// `{P(runs | Y = 0), P(runs | Y = 1)}` for the cell of `molecule` in a species whose runs are
/// `runs`: the product over the runs of their detection terms. `{1, 1}` for a species without
/// runs, so mass spec changes nothing where there is no data.
template<typename FilterProbability>
[[nodiscard]] std::array<double, 2> per_cell_factor(coretools::TConstView<TMassSpecRun> runs,
                                                    uint32_t molecule,
                                                    const FilterProbability &filter_probability,
                                                    double contamination) {
	std::array<double, 2> factor{1.0, 1.0};
	for (const auto &run : runs) {
		if (run.empty()) { continue; }
		const bool detected = run.is_detected(molecule);
		const double pass   = filter_probability(run.filter_index(), molecule);
		factor[0] *= TMassSpecRun::probability_of_assignment(false, detected, pass, contamination);
		factor[1] *= TMassSpecRun::probability_of_assignment(true, detected, pass, contamination);
	}
	return factor;
}

/// The whole mass-spec log-likelihood: every feature term plus, for every run, the detection term
/// of every molecule leaf. `runs_per_species[s]` are the runs of species leaf `s`, and `Y` is laid
/// out `{species leaf, molecule leaf}`.
///
/// The detection term spans the whole molecule space, but almost every molecule is absent and not
/// detected, so that case, `log(1 - contamination)`, is added in bulk; only the present and the
/// detected molecules are scored one by one. That makes the cost O(ones(Y) + features) rather than
/// O(runs x molecule leaves).
template<BinaryStorage Field, typename FilterProbability>
[[nodiscard]] double log_likelihood(const coretools::TNestedVector<TMassSpecRun> &runs_per_species,
                                    const Field &Y, const FilterProbability &filter_probability,
                                    double contamination) {
	const size_t n_species   = runs_per_species.size();
	const size_t n_molecules = Y.dimensions()[1];
	if (Y.dimensions()[0] != n_species) {
		throw coretools::TDevError("The field has ", Y.dimensions()[0], " species leaves but ",
		                           n_species, " species carry mass-spec runs.");
	}

	// the present molecules of each species, walked through the ones cursor in O(ones(Y))
	std::vector<std::vector<uint32_t>> present_per_species(n_species);
	for (auto cursor = Y.ones_cursor(); cursor.valid(); cursor.advance()) {
		const auto cell = Y.get_multi_dimensional_index(cursor.linear_index());
		present_per_species[cell[0]].push_back(static_cast<uint32_t>(cell[1]));
	}

	const double log_absent_not_detected = std::log(1.0 - contamination);
	coretools::TSumLogProbability sum_log{};
	double bulk_log = 0.0;
	for (size_t species = 0; species < n_species; ++species) {
		const auto &present = present_per_species[species];
		for (const auto &run : runs_per_species.at(species)) {
			if (run.empty()) { continue; }
			const size_t filter = run.filter_index();

			// feature terms, and the absent molecules this run detects
			size_t n_scored = 0;
			for (size_t f = 0; f < run.number_of_assignments(); ++f) {
				const auto &assignment = run.get_current_assignment(f);
				sum_log.add(feature_term(assignment));
				if (assignment.is_unknown_molecule()) { continue; }
				const uint32_t m = assignment.get_molecule_index();
				if (Y.is_one(Y.get_linear_index_in_container_space({species, m}))) { continue; }
				sum_log.add(TMassSpecRun::probability_of_assignment(
				    false, true, filter_probability(filter, m), contamination));
				++n_scored;
			}
			// the present molecules, detected or not
			for (const uint32_t m : present) {
				sum_log.add(TMassSpecRun::probability_of_assignment(
				    true, run.is_detected(m), filter_probability(filter, m), contamination));
				++n_scored;
			}
			bulk_log += static_cast<double>(n_molecules - n_scored) * log_absent_not_detected;
		}
	}
	return sum_log.getSum() + bulk_log;
}

} // namespace mass_spec_math
