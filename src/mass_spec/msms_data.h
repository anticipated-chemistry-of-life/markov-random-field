#pragma once

#include "Types.h"
#include "constants.h"
#include "coretools/Containers/TNestedVector.h"
#include "coretools/Containers/TView.h"
#include "coretools/Main/TError.h"
#include "coretools/Math/TSumLog.h"
#include "coretools/algorithms.h"
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
class TMSMSData : public stattools::prior::TBaseLikelihoodPrior<stattools::TObservationBase,
                                                                TypeMSData, NumDimMSData> {
public:
	// some type aliases, for better readability
	using BoxType = TMSMSData;
	using Base    = stattools::prior::TBaseLikelihoodPrior<stattools::TObservationBase, TypeMSData,
	                                                       NumDimMSData>;
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
	std::array<size_t, 2> _dimensions_for_filters{};
	size_t _number_of_filters;
	const TMarkovField &_markov_field;

	TypeParamContamination *_proba_contamination   = nullptr;
	TypeParamMassSpecFilter *_proba_to_pass_filter = nullptr;

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
	_get_filter_molecule_pair_multidimensional_index_(size_t filter_molecule_pair_index) const {
		return coretools::getSubscriptsAsArray(filter_molecule_pair_index, _dimensions_for_filters);
	}

	/// Complete log-likelihood of all the mass-spec data given the current assignments, the current
	/// Y (molecule presence/absence), and the current filter/contamination parameters. Per run it
	/// is
	///   sum_features  log P(feature | assigned molecule)              [feature term, PHRED table]
	/// + sum_molecules log probability_of_assignment(Y, ms, filter, c) [detection term, all
	/// molecules] where `ms` = "molecule is assigned to some feature in this run". The detection
	/// term spans the whole molecule space; the dominant (absent & unassigned -> 1-contamination)
	/// case is added as a /// single bulk term, only present or assigned molecules are scored
	/// individually. The contamination probability is passed in (rather than read from the
	/// parameter) so the same routine can score both the old and the proposed contamination value
	/// in a ratio.
	[[nodiscard]] double _calculate_log_likelihood_of_MSData(double contamination) const;

public:
	explicit TMSMSData(
	    const std::vector<std::unique_ptr<TTree>> &trees, const TMarkovField &markov_field,
	    size_t number_of_filters,
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

	// Add the MS-data log-likelihood contribution of the Y cell (species, molecule) for both Y
	// states into sum_log[0] (molecule absent) and sum_log[1] (molecule present). This is the hook
	// used by the Y Gibbs update, so only the terms that depend on Y[species, molecule] are scored:
	// the per-run detection term (filter / contamination) for this molecule given its current
	// assignment status in each run. The feature term P(feature | assigned molecule) does not
	// depend on Y and would cancel between the two states, so it is intentionally omitted here.
	void add_log_likelihood(const IndexArray &indices_in_leaves,
	                        std::array<coretools::TSumLogProbability, 2> &sum_log) const {
		const auto species_idx  = indices_in_leaves[_species_dim];
		const auto molecule_idx = indices_in_leaves[_molecule_dim];
		const double cont       = (double)_proba_contamination->value();
		for (const auto &run : this->get_ms_data_for_species(species_idx)) {
			if (run.size() == 0) { continue; }
			const bool ms     = run.is_molecule_assigned(molecule_idx);
			const double filt = LINEAR_SPACE_PROBA[_proba_to_pass_filter->value(
			    _get_linear_index_filter_molecule_pair({run.filter_index(), molecule_idx}))];
			sum_log[0].add(TMassSpecRun::probability_of_assignment(false, ms, filt, cont)); // Y = 0
			sum_log[1].add(TMassSpecRun::probability_of_assignment(true, ms, filt, cont));  // Y = 1
		}
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
	[[nodiscard]] double calculateLLRatio(TypeParamContamination *, size_t /*Index*/) {
		const double old_LL =
		    _calculate_log_likelihood_of_MSData((double)_proba_contamination->oldValue());
		const double cur_LL =
		    _calculate_log_likelihood_of_MSData((double)_proba_contamination->value());
		return cur_LL - old_LL;
	};

	double calculateLLRatio(TypeParamMassSpecFilter *, size_t index);

	/// Log-likelihood ratio (new - old) of a proposed assignment move in a single run. Y and the
	/// mass-spec filter/contamination parameters are held fixed; only the terms that the move
	/// changes are scored:
	///   * the feature term(s): P(feature | assigned molecule), the assigned binned likelihood
	///   mapped
	///     through `_log_lik_present` (for the unknown molecule the per-feature unknown probability
	///     is mapped through the same table);
	///   * the filter/contamination term(s) for the real molecules whose "is assigned" status flips
	///     (none flip for a Swap, so only feature terms change there), via
	///     probability_of_assignment.
	/// Returns 0 for an invalid (no-op) proposal. This is the likelihood ratio only; any
	/// proposal/Hastings ratio must be added by the caller.
	[[nodiscard]] double
	calculate_LL_ratio_for_assignment_move(size_t species_idx, const TMassSpecRun &run,
	                                       const TAssignmentProposal &move) const;

	/// TODO: In the current implementation there would be only one update/proposal per MCMC
	/// iteration. We could also propose multiple moves per iteration using an additional loop.
	void update_all_MS_assignments();
};
