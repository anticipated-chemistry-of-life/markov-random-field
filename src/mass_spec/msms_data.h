#pragma once

#include "Types.h"
#include "coretools/Containers/TNestedVector.h"
#include "coretools/Containers/TView.h"
#include "coretools/Main/TError.h"
#include "coretools/Math/TSumLog.h"
#include "coretools/algorithms.h"
#include "mass_spec/msms_run.h"
#include "tree/TTree.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

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
class TMSMSData : public stattools::prior::TBaseLikelihoodPrior<TypeMSData, NumDimMSData> {
public:
	// some type aliases, for better readability
	using BoxType = TMSMSData;
	using Base    = stattools::prior::TBaseLikelihoodPrior<TypeMSData, NumDimMSData>;
	using typename Base::Storage;
	using typename Base::UpdatedStorage;

	using TypeParamMassSpecFilter = stattools::TParameter<SpecMassSpecFilter, BoxType>;
	using TypeParamContamination  = stattools::TParameter<SpecContaminationProba, BoxType>;

private:
	coretools::TNestedVector<TMassSpecRun> _msms_data;
	TTree *_molecules_tree = nullptr;
	TTree *_species_tree   = nullptr;
	size_t _molecule_dim   = 0;
	size_t _species_dim    = 0;
	std::array<size_t, 2> _dimensions_for_filters;
	size_t _number_of_filters;

	// log P(observation | molecule absent/present), indexed by binned_value (0-255).
	// Initialized to 0.0 (= log(1), no contribution) until populated from file.
	std::array<double, 256> _log_lik_absent{}; // TODO: maybe remove and just add 1-p(score| x=1) ?
	std::array<double, 256> _log_lik_present{};

	TypeParamContamination *_contamination_probability = nullptr;
	TypeParamMassSpecFilter *_proba_to_pass_filter     = nullptr;
	double _oldLL;
	double _curLL;

	// Markov field parameter (only needed for stattools purposes to build a valid DAG)
	const std::vector<std::unique_ptr<stattools::TParameter<SpecMarkovField, TLotus>>>
	    &_markov_field_stattools_param;

	// Private functions
private:
	void _add_mass_spec_runs_for_species(const std::vector<TMassSpecRun> &runs) {
		_msms_data.push_back(runs);
	}
	void _check_trees_exist() const {
		if (!_species_tree || !_molecules_tree) {
			throw coretools::TDevError("TMSMSData : missing required tree(s) ");
		}
	}

	[[nodiscard]] size_t _get_linear_index_filter_molecule_pair(
	    const std::array<size_t, 2> &filter_molecule_pair_index) const {
		return coretools::getLinearIndex(filter_molecule_pair_index, _dimensions_for_filters);
	}

	[[nodiscard]] std::array<size_t, 2>
	_get_filter_molecule_paire_multidimensional_index_(size_t filter_molecule_pair_index) const {
		return coretools::getSubscriptsAsArray(filter_molecule_pair_index, _dimensions_for_filters);
	}

public:
	explicit TMSMSData(
	    const std::vector<std::unique_ptr<TTree>> &trees, size_t number_of_filters,
	    const std::vector<std::unique_ptr<stattools::TParameter<SpecMarkovField, TLotus>>>
	        &markov_field_stattools_param,
	    TypeParamMassSpecFilter *filter_proba, TypeParamContamination *contamination_proba);
	~TMSMSData() override = default;

	// stattools overrides
	[[nodiscard]] std::string name() const override { return "msmsdata"; }
	void _simulateUnderPrior(Storage *) override {
		throw coretools::TDevError("TMSMSData: simulation under prior is not implemented.");
	}
	void initialize() override;
	void guessInitialValues() override;

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

	[[nodiscard]] double calculateLLRatio(TypeParamContamination *, size_t /*Index*/,
	                                      const Storage &) {
		_curLL = 0.0;
		_oldLL = 0.0;
		throw coretools::TDevError("Function not implemented yet");
	};
	double calculateLLRatio(TypeParamMassSpecFilter *, size_t /*Index*/, const Storage &) {
		_curLL = 0.0;
		_oldLL = 0.0;

		throw coretools::TDevError("Function not implemented yet");
	};
	void updateTempVals(TypeParamContamination *, size_t /*Index*/, bool Accepted) {
		if (!Accepted) {
			_curLL = _oldLL; // reset
		}
	};
};
