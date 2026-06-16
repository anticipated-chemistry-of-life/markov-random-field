#pragma once

/// Class storing a mass spectrometry run. The iterator goes along the
/// molecule dimension. So for each feature, we have a list of molecule indices and binned
/// likelihoods.
#include "./feature_likelihood.h"
#include "Types.h"
#include "coretools/Containers/TNestedVector.h"
#include "coretools/Main/TError.h"
#include "coretools/Math/TSumLog.h"
#include "coretools/devtools.h"
#include <Security/cssmconfig.h>
#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

struct BinarySearchResult {
	bool found{};
	std::optional<uint8_t> binned_likelihood;
	size_t index{};
};

class TMassSpecRun {
private:
	coretools::TNestedVector<TFeatureLikelihood> _features;

	/// Each feature also needs a probability for the unknown molecules. In case a feature has as
	/// true molecule one that is in the unknown set, the probability of that molecule is used.
	/// This vector should be of size "number of features" in the MSMS run.
	std::vector<uint8_t> _probabilities_of_unkowns;

	/// The current assignment of molecules to features. The vector should be of size "number of
	/// features" in the current MSMS run.
	std::vector<TFeatureLikelihood> _current_assignments;

	/// The filter index for the current MSMS run.
	size_t _filter_index;

	// const stattools::TParameter<SpecContaminationProba, TMSMSData> *_proba_of_contamination;
	// const stattools::TParameter<SpecMassSpecFilter, TMSMSData> *_proba_to_pass_filter;

public:
	TMassSpecRun() = default;
	[[nodiscard]] size_t capacity() const { return _features.capacity(); }
	void reserve(size_t n) { _features.reserve(n); }
	[[nodiscard]] bool empty() const { return _features.empty(); }
	void add_empty_feature() { _features.push_back(); }

	[[nodiscard]] coretools::TConstView<TFeatureLikelihood>
	get_likelihoods_for_feature(size_t i) const {
		return _features.at(i);
	}
	[[nodiscard]] size_t size() const { return _features.size(); }
	[[nodiscard]] size_t number_of_features() const { return this->size(); }
	[[nodiscard]] auto begin() const { return _features.begin(); };
	[[nodiscard]] auto end() const { return _features.end(); }
	auto begin() { return _features.begin(); }
	auto end() { return _features.end(); }
	/// Returns the filter index for the current MSMS run.
	[[nodiscard]] size_t filter_index() const { return _filter_index; }
	///
	[[nodiscard]] BinarySearchResult is_molecule_in_feature(size_t feature_idx,
	                                                        uint32_t molecule_index) const {
		const auto &likelihoods = _features.at(feature_idx);
		return is_molecule_in_feature(likelihoods, molecule_index);
	}

	inline static BinarySearchResult
	is_molecule_in_feature(const coretools::TConstView<TFeatureLikelihood> &feature,
	                       uint32_t molecule_index) {
		auto it = std::lower_bound(feature.begin(), feature.end(), molecule_index);
		if (it == feature.end()) { return {false, std::nullopt, feature.size()}; }
		size_t index = std::distance(feature.begin(), it);
		if (it->get_molecule_index() != molecule_index) { return {false, std::nullopt, index}; }
		return {true, it->get_binned_likelihood(), index};
	}

	static double
	get_likelihood_from_binned_value(const std::array<double, 256> &binned_likelihoods,
	                                 uint8_t binned_value) {
		return binned_likelihoods[binned_value];
	}

	// Treats features as independent (approximation: features within a run share molecules).
	void
	add_log_likelihood_for_molecule(uint32_t molecule_idx,
	                                const std::array<double, 256> &log_lik_absent,
	                                const std::array<double, 256> &log_lik_present,
	                                std::array<coretools::TSumLogProbability, 2> &sum_log) const {
		const auto zero = sum_log[0].getSum();
		const auto one  = sum_log[1].getSum();
		OUT(zero, one);
		throw coretools::TDevError(
		    "We still don't know of to add the probability for the same molecule that "
		    "appears in two features of the same MSMS run. Here are the variables of the function "
		    ": ",
		    log_lik_absent, log_lik_present);
		for (size_t f = 0; f < _features.size(); ++f) {
			auto result = is_molecule_in_feature(f, molecule_idx);
			if (result.found) {
				throw coretools::TDevError(
				    "We still don't know of to add the probability for the same molecule that "
				    "appears in two features of the same MSMS run.");
			}
		}
	}

	/// Adds a vector of molecule likelihoods. This represents one feature with all the molecule
	/// likelihoods that are associated to it. The vector will be sorted before being added to
	/// guarantee the order and be able to do binary search.
	void add_likelihood_vector(std::vector<TFeatureLikelihood> &feature_likelihoods) {
		std::sort(feature_likelihoods.begin(), feature_likelihoods.end());
		_features.push_back(feature_likelihoods);
	}

	[[nodiscard]] bool feature_has_unknown_molecule_assigned(size_t feature_idx) const {
		return _current_assignments.at(feature_idx).get_molecule_index() >=
		       TFeatureLikelihood::get_unknown_molecule_index();
	}

	/// Calculates the probability if a molecule is assigned to a feature based on the
	/// observed MS/MS data and the model's probability estimates.
	/// If the current state in Y for a molecule-species pair is 1 it will check if the
	/// molecule is assigned to the feature and return the corresponding probability.
	///
	/// This value will then be multiplied by the likelihood of the molecule being present
	/// in the feature to get the final probability.
	/// TMassSpecRun — pure: caller passes the filter prob (for this molecule) and contamination
	[[nodiscard]] static double probability_of_assignment(bool y, bool ms,
	                                                      double proba_to_pass_filter,
	                                                      double proba_contamination) {
		if (y && ms) { return proba_to_pass_filter; }        // present & detected -> passed filter
		if (y && !ms) { return 1.0 - proba_to_pass_filter; } // present & not detected
		if (!y && ms) { return proba_contamination; }        // absent  & detected -> contamination
		return 1.0 - proba_contamination;                    // absent  & not detected
	}

	[[nodiscard]] bool is_molecule_assigned(uint32_t molecule_idx) const {
		return std::any_of(_current_assignments.begin(), _current_assignments.end(),
		                   [molecule_idx](const TFeatureLikelihood &a) {
			                   return a.get_molecule_index() == molecule_idx;
		                   });
	}
};
