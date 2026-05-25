#pragma once

#include "./msms_run.h"
#include "coretools/Containers/TNestedVector.h"
#include "coretools/Containers/TView.h"
#include "coretools/Main/TError.h"
#include "coretools/Math/TSumLog.h"
#include "tree/TTree.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>

/// Class TRun that stores 2 vectors :
/// - vector of 32 bit with value uin8t which is the binned likelihood and molecule index
/// - vector of size_t that tells where the feature starts
/// Use TNestedVector from coretools
///
/// TMSMSData (which is also a wrapper around nested vector)
/// - vector of TRun
/// - nested vecctor

/// Need to sort species. First species with MS data, then the ones with no data.

/// Class that stores all the mass spec data. This class will iterate over all the species and
/// return the mass spec runs for a species of interest. WE will need to sort the species in order
/// to have first the species WITH mass spec data, then the species without mass spec data. This
/// will avoid have the indices being repeated and store a size_t for each species even it has no
/// data.
class TMSMSData {
private:
	coretools::TNestedVector<TMassSpecRun> _msms_data;
	TTree *_molecules_tree;
	TTree *_species_tree;

private:
	void _add_mass_spec_runs_for_species(const std::vector<TMassSpecRun> &runs) {
		_msms_data.push_back(runs);
	}

public:
	explicit TMSMSData(const std::vector<std::unique_ptr<TTree>> &trees);
	~TMSMSData() = default;

	[[nodiscard]] bool empty() const { return _msms_data.empty(); }

	void add_to_sumlog(coretools::TSumLogProbability &sum_log, uint8_t binned_value) const;
	void add_to_sumlog(std::array<coretools::TSumLogProbability, 2> &sum_log,
	                   uint8_t binned_value) const;

	[[nodiscard]] coretools::TConstView<TMassSpecRun>
	get_ms_data_for_species(size_t species_idx) const {
		// with bound check
		return _msms_data.at(species_idx);
	};
	[[nodiscard]] size_t number_of_ms_runs_for_species(size_t species_idx) const {
		// with bound check
		return _msms_data.at(species_idx).size();
	};

	[[nodiscard]] coretools::TConstView<TMassSpecRun>
	get_ms_data_for_species(const std::string &species_name) const {
		size_t species_idx = _species_tree->get_index_within_leaves(species_name);
		return this->get_ms_data_for_species(species_idx);
	};

	[[nodiscard]] size_t get_molecule_index_within_leaves(const std::string &molecule_name) const {
		return _molecules_tree->get_index_within_leaves(molecule_name);
	}
};
