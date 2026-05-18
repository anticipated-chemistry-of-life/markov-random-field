#pragma once

#include "./msms_run.h"
#include "coretools/Containers/TNestedVector.h"
#include "coretools/Math/TSumLog.h"
#include <array>
#include <cstddef>
#include <cstdint>

/// Class TRun that stores 2 vectors :
/// - vector of 32 bit with value uin8t which is the binned likelihood and molecule index
/// - vector of size_t that tells where the feature starts
/// Use TNestedVector from coretools
///
/// TMSMSData (which is also a wrapper around nested vector)
/// - vector of TRun
/// - nested vecctor

/// Need to sort species. First species with MS data, then the ones with no data.

/// Class that stores all the mass spec data. This class will iterate over all the species and return the mass spec runs
/// for a species of interest. WE will need to sort the species in order to have first the species WITH mass spec data,
/// then the species without mass spec data. This will avoid have the indices being repeated and store a size_t for each
/// species even it has no data.
class TMSMSData {
private:
	coretools::TNestedVector<TMassSpecRun> _msms_data;
	std::array<double, 256> _binned_likelihoods{};
	void _add_mass_spec_run_for_species(const std::vector<TMassSpecRun> &runs) { _msms_data.push_back(runs); }

public:
	TMSMSData() = default;
	[[nodiscard]] const std::array<double, 256> &get_binned_likelihoods() const { return _binned_likelihoods; }
	[[nodiscard]] double get_likelihood_from_binned_value(uint8_t binned_value) const {
		return _binned_likelihoods[binned_value];
	}
	[[nodiscard]] bool empty() const { return _msms_data.empty(); }
	/// Note : This function assumes that the species were sorted to have first the ones with data
	/// and then the ones without data i.e. the indices of the NestedVector must not contain any repeated indices.
	[[nodiscard]] bool species_has_msms_data(size_t species_idx) const { return species_idx < _msms_data.size(); }
	void add_to_sumlog(coretools::TSumLogProbability &sum_log, uint8_t binned_value) const {
		sum_log.add(get_likelihood_from_binned_value(binned_value));
	}
	void add_to_sumlog(std::array<coretools::TSumLogProbability, 2> &sum_log, uint8_t binned_value) const {
		sum_log[0].add(1.0 - get_likelihood_from_binned_value(binned_value));
		sum_log[1].add(get_likelihood_from_binned_value(binned_value));
	}
};
