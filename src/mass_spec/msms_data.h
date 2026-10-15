#pragma once

#include "Types.h"
#include "constants.h"
#include "coretools/Containers/TNestedVector.h"
#include "coretools/Containers/TView.h"
#include "coretools/Main/TError.h"
#include "coretools/Math/TSumLog.h"
#include "coretools/algorithms.h"
#include "mass_spec/TMassSpecMath.h"
#include "mass_spec/msms_run.h"
#include "stattools/ParametersObservations/TParameter.h"
#include "tree/TTree.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

class TMarkovField;
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
///
/// This class is not a stattools box: its parameters (the filter and contamination
/// probabilities) hang off TDataModel, which owns this object and forwards the MCMC callbacks --
/// the same pattern as TLotus and TSimpleErrorModel (see TDataModel.h). Unlike them it does not
/// satisfy the DataSource concept (data_sources/data_source.h): it scores an assignment problem
/// over a species' runs, not a per-cell state, and it carries its own latent state (the
/// feature<->molecule assignments) that needs its own Metropolis sweep every iteration
/// (update_all_MS_assignments), so TDataModel hand-wires it rather than folding it into
/// ACOL_DATA_SOURCES.
class TMSMSData {
public:
	using TypeParamMassSpecFilter = stattools::TParameter<SpecMassSpecFilter, TDataModel>;
	using TypeParamContamination  = stattools::TParameter<SpecContaminationProba, TDataModel>;

private:
	coretools::TNestedVector<TMassSpecRun> _msms_data;
	TTree *_molecules_tree = nullptr;
	TTree *_species_tree   = nullptr;
	size_t _molecule_dim   = 0;
	size_t _species_dim    = 0;
	std::array<size_t, 2> _dimensions_for_filters{};
	size_t _number_of_filters;
	size_t _number_of_species;
	const TMarkovField &_markov_field;

	TypeParamContamination *_proba_contamination   = nullptr;
	TypeParamMassSpecFilter *_proba_to_pass_filter = nullptr;

	std::vector<uint8_t> _species_has_ms_data;

	/// Private functions
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
	_get_filter_molecule_pair_multidimensional_index_(size_t filter_molecule_pair_index) const {
		return coretools::getSubscriptsAsArray(filter_molecule_pair_index, _dimensions_for_filters);
	}

	/// The probability that a present molecule passes `filter`, in linear space.
	[[nodiscard]] double _filter_probability(size_t filter, uint32_t molecule) const {
		return LINEAR_SPACE_PROBA[_proba_to_pass_filter->value(
		    _get_linear_index_filter_molecule_pair({filter, molecule}))];
	}

	/// The whole mass-spec log-likelihood at `contamination`, given the current field, assignments
	/// and filter probabilities (`mass_spec_math::log_likelihood`). The contamination is passed in
	/// rather than read from the parameter, so the old and the proposed value can both be scored.
	[[nodiscard]] double _log_likelihood(double contamination) const;

public:
	explicit TMSMSData(const std::vector<std::unique_ptr<TTree>> &trees,
	                   const TMarkovField &markov_field, size_t number_of_filters,
	                   TypeParamMassSpecFilter *filter_proba,
	                   TypeParamContamination *contamination_proba);
	~TMSMSData() = default;

	/// Sizes the filter-probability and contamination parameter storages. `box` is the TDataModel
	/// the parameters hang off; it has to be passed in because this class is no longer the box
	/// itself.
	void initialize(TDataModel *box);
	void guess_initial_values();

	[[nodiscard]] bool empty() const { return _msms_data.empty(); }

	// Add the MS-data factor of the field cell (species, molecule) for both field states into
	// sum_log[0] (molecule absent) and sum_log[1] (molecule present), i.e.
	// `mass_spec_math::per_cell_factor`. The feature terms do not read the field and cancel between
	// the two states, so the factor leaves them out.
	//
	// Still dormant: the block update does not read this either (TMarkovField's eight-state block
	// takes the LOTUS and the simple-error term), and adapting a third source is that source's own
	// work.
	void add_log_likelihood(const IndexArray &indices_in_leaves,
	                        std::array<coretools::TSumLogProbability, 2> &sum_log) const {
		const auto species_idx = indices_in_leaves[_species_dim];

		// if the species does not have MS data, skip it
		if (!_species_has_ms_data.at(species_idx)) { return; }
		const auto factor = mass_spec_math::per_cell_factor(
		    this->get_ms_data_for_species(species_idx),
		    static_cast<uint32_t>(indices_in_leaves[_molecule_dim]),
		    [this](size_t filter, uint32_t molecule) {
			    return _filter_probability(filter, molecule);
		    },
		    (double)_proba_contamination->value());
		sum_log[0].add(factor[0]);
		sum_log[1].add(factor[1]);
	}

	[[nodiscard]] coretools::TConstView<TMassSpecRun>
	get_ms_data_for_species(size_t species_idx) const {
		this->_check_trees_exist();
		return _msms_data.at(species_idx);
	};

	[[nodiscard]] coretools::TConstView<TMassSpecRun>
	get_ms_data_for_species(const std::string &species_name) const {
		const size_t species_idx = _species_tree->get_index_within_leaves(species_name);
		return this->get_ms_data_for_species(species_idx);
	};

	[[nodiscard]] coretools::TView<TMassSpecRun> get_ms_data_for_species(size_t species_idx) {
		this->_check_trees_exist();
		return _msms_data.at(species_idx);
	};

	[[nodiscard]] size_t number_of_ms_runs_for_species(size_t species_idx) const {
		this->_check_trees_exist();
		return _msms_data.at(species_idx).size();
	};

	[[nodiscard]] size_t get_molecule_index_within_leaves(const std::string &molecule_name) const {
		this->_check_trees_exist();
		return _molecules_tree->get_index_within_leaves(molecule_name);
	}

	/// Like gamma and the error rate, the contamination probability is a single scalar that enters
	/// (almost) every term of the MS-data likelihood, so we score the whole thing rather than a
	/// localized ratio. The likelihood also depends on Y, the filter probabilities and the
	/// assignments, all of which change in other updates; rather than cache the value and keep it
	/// in sync, we recompute the full likelihood at both the old and the proposed contamination
	/// value here, which is always correct regardless of the order of updates.
	[[nodiscard]] double ll_ratio_after_contamination_move() const {
		return _log_likelihood((double)_proba_contamination->value()) -
		       _log_likelihood((double)_proba_contamination->oldValue());
	}

	[[nodiscard]] double ll_ratio_after_filter_move(size_t index) const;

	/// Log-likelihood ratio (new - old) of a proposed assignment move in a single run. Y and the
	/// mass-spec filter/contamination parameters are held fixed; only the terms that the move
	/// changes are scored:
	///   * the feature term(s): P(feature | assigned molecule), the assigned binned likelihood
	///   mapped
	///     through `_log_lik_present` (for the unknown molecule the per-feature unknown probability
	///     is mapped through the same table);
	///   * the filter/contamination term(s) for the real molecules whose "is assigned" status flips
	///     (none flip for either swap type, so only feature terms change there), via
	///     probability_of_assignment.
	/// Returns 0 for an invalid (no-op) proposal. This is the likelihood ratio only; any
	/// proposal/Hastings ratio must be added by the caller.
	[[nodiscard]] double
	calculate_LL_ratio_for_assignment_move(size_t species_idx, const TMassSpecRun &run,
	                                       const TAssignmentProposal &move) const;

	/// One Metropolis-Hastings update over every feature of every run of every species that has MS
	/// data. See `TMassSpecRun::update_all_assignments`.
	void update_all_MS_assignments();
};
