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
	TTree *_molecules_tree = nullptr;
	TTree *_species_tree   = nullptr;
	size_t _molecule_dim   = 0;
	size_t _species_dim    = 0;

	// log P(observation | molecule absent/present), indexed by binned_value (0-255).
	// Initialized to 0.0 (= log(1), no contribution) until populated from file.
	std::array<double, 256> _log_lik_absent{}; // TODO: maybe remove and just add 1-p(score| x=1) ?
	std::array<double, 256> _log_lik_present{};

	void _add_mass_spec_runs_for_species(const std::vector<TMassSpecRun> &runs) {
		_msms_data.push_back(runs);
	}
	void _check_trees_exist() const {
		if (!_species_tree || !_molecules_tree) {
			throw coretools::TDevError("TMSMSData : missing required tree(s) ");
		}
	}

public:
	explicit TMSMSData(const std::vector<std::unique_ptr<TTree>> &trees);
	~TMSMSData() = default;

	[[nodiscard]] bool empty() const { return _msms_data.empty(); }

	// Populate lookup tables from pre-computed log-probability arrays (called during data loading).
	void set_lookup_tables(const std::array<double, 256> &log_lik_absent,
	                       const std::array<double, 256> &log_lik_present) {
		_log_lik_absent  = log_lik_absent;
		_log_lik_present = log_lik_present;
	}

	void add_to_sumlog(coretools::TSumLogProbability &sum_log, uint8_t binned_value) const;
	void add_to_sumlog(std::array<coretools::TSumLogProbability, 2> &sum_log,
	                   uint8_t binned_value) const;

	// Sum log P(MS data for species | molecule = absent/present) into sum_log[0/1].
	// Runs are treated as independent; features within a run as independent (approximation).
	void add_log_likelihood(const std::vector<size_t> &indices_in_leaves,
	                        std::array<coretools::TSumLogProbability, 2> &sum_log) const {
		const auto species_idx  = indices_in_leaves[_species_dim];
		const auto molecule_idx = indices_in_leaves[_molecule_dim];
		const auto runs         = this->get_ms_data_for_species(species_idx);
		for (const auto &run : runs) {
			run.add_log_likelihood_for_molecule(static_cast<uint32_t>(molecule_idx),
			                                    _log_lik_absent, _log_lik_present, sum_log);
		}
	}

	[[nodiscard]] coretools::TConstView<TMassSpecRun>
	get_ms_data_for_species(size_t species_idx) const {
		this->_check_trees_exist();
		return _msms_data.at(species_idx);
	};
	[[nodiscard]] size_t number_of_ms_runs_for_species(size_t species_idx) const {
		this->_check_trees_exist();
		return _msms_data.at(species_idx).size();
	};

	[[nodiscard]] coretools::TConstView<TMassSpecRun>
	get_ms_data_for_species(const std::string &species_name) const {
		size_t species_idx = _species_tree->get_index_within_leaves(species_name);
		return this->get_ms_data_for_species(species_idx);
	};

	[[nodiscard]] size_t get_molecule_index_within_leaves(const std::string &molecule_name) const {
		this->_check_trees_exist();
		return _molecules_tree->get_index_within_leaves(molecule_name);
	}
};
