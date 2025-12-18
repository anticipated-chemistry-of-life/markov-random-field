#pragma once

#include "coretools/Containers/TNestedIterator.h"
#include "coretools/Containers/TNestedVector.h"
#include "coretools/Main/TError.h"
#include "coretools/Math/TSumLog.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

struct BinarySearchResult {
	bool found;
	size_t index;
};

/// Class TRun that stores 2 vectors :
/// - vector of 32 bit with value uin8t and molecule index
/// - vetor of size_t that tells where the feature starts
/// Use TNestedVector from coretools
///
/// TMSMSData (which is also a wrapper around nested vector)
/// - vector of TRun
/// - nested vecctor

/// Need to sort species. First species with MS data, then the ones with no data.

/// Class that stores the likelihood of a feature for a defined molecule
class TFeatureLikelihood {
private:
	uint32_t _value = 0;

	/// A mask for the 255 binned likelihoods
	static constexpr uint32_t _one                 = 1;
	static constexpr uint32_t _likelihood_mask     = ~((_one << 24) - 1);
	/// The remaining bits for the linear index
	/// This also defines the maximum number of molecules that can be stored
	static constexpr uint32_t _molecule_index_mask = (_one << 24) - 1;

public:
	TFeatureLikelihood() = default;
	explicit TFeatureLikelihood(const uint32_t linear_index, const uint8_t binned_likelihood) {
		set_molecule_index(linear_index);
		set_binned_likelihood(binned_likelihood);
	}
	uint32_t get_molecule_index() const { return _value & _molecule_index_mask; };
	void set_molecule_index(const uint32_t linear_index) {
		if (linear_index > _molecule_index_mask) {
			throw coretools::TUserError("Molecule index '", linear_index,
			                            "' exceeds the maximum allowed number of molecules : ", _molecule_index_mask,
			                            ".");
		}
		_value = (_value & ~_molecule_index_mask) | (linear_index & _molecule_index_mask);
	}
	uint32_t get_binned_likelihood() const { return (_value & _likelihood_mask) >> 24; }
	void set_binned_likelihood(const uint8_t binned_likelihood) {
		_value = (_value & ~_likelihood_mask) | (static_cast<uint32_t>(binned_likelihood) << 24);
	}

	bool operator<(const TFeatureLikelihood &right) const { return get_molecule_index() < right.get_molecule_index(); }
	bool operator<(const uint32_t right) const { return get_molecule_index() < right; }
	bool operator!=(const uint32_t right) const { return get_molecule_index() != right; }
	bool operator==(const uint32_t right) const { return get_molecule_index() == right; }
};

/// Class storing a mass spectrometry run. The iterator goes along the
/// molecule dimension. So for each feature, we have a list of molecule indices and binned likelihoods.
class TMassSpecRun {
private:
	coretools::TNestedVector<TFeatureLikelihood> _features;

public:
	TMassSpecRun() = default;
	size_t capacity() const { return _features.capacity(); }
	void reserve(size_t n) { _features.reserve(n); }
	bool empty() const { return _features.empty(); }
	void add_empty_feature() { _features.push_back(); }

	/// Adds a molecule likelihood to the last feature (back feature)
	void add_molecule_likelihood_to_last_feature(const TFeatureLikelihood &feature_likelihood) {
		_features.push_back(feature_likelihood);
	}

	/// Adds a vector of molecule likelihoods. This represents one feature with all the molecule likelihoods that are
	/// associated to it.
	void add_likelihood_vector(const std::vector<TFeatureLikelihood> &feature_likelihoods) {
		_features.push_back(feature_likelihoods);
	}

	coretools::TConstView<TFeatureLikelihood> get_likelihoods_for_feature(size_t i) const { return _features.get(i); }
	size_t size() const { return _features.size(); }
	auto begin() const { return _features.begin(); };
	auto end() const { return _features.end(); }
	auto begin() { return _features.begin(); }
	auto end() { return _features.end(); }
	BinarySearchResult is_molecule_in_feature(size_t feature_idx, uint32_t molecule_index) const {
		const auto &likelihoods = _features.get(feature_idx);
		auto it                 = std::lower_bound(likelihoods.begin(), likelihoods.end(), molecule_index);
		if (it == likelihoods.end()) { return {false, likelihoods.size()}; }
		size_t index = std::distance(likelihoods.begin(), it);
		if (it->get_molecule_index() != molecule_index) { return {false, index}; }
		return {true, index};
	}

	const TFeatureLikelihood *get_likelihood_for_molecule_in_feature(size_t feature_idx,
	                                                                 uint32_t molecule_index) const {
		const auto &likelihoods = _features.get(feature_idx);
		auto it                 = std::lower_bound(likelihoods.begin(), likelihoods.end(), molecule_index);
		return it;
	}
};

/// Class that stores all the mass spec data. This class will iterate over all the species and return the mass spec runs
/// for a species of interest. WE will need to sort the species in order to have first the species WITH mass spec data,
/// then the species without mass spec data. This will avoid have the indices being repeated and store a size_t for each
/// species even it has no data.
class TMSMSData {
private:
	coretools::TNestedVector<TMassSpecRun> _msms_data;
	std::array<double, 256> _binned_likelihoods;

public:
	TMSMSData() = default;
	const std::array<double, 256> &get_binned_likelihoods() const { return _binned_likelihoods; }
	double get_likelihood_from_binned_value(uint8_t binned_value) const { return _binned_likelihoods[binned_value]; }
	bool empty() const { return _msms_data.empty(); }
	/// Note : This function assumes that the species were sorted to have first the ones with data
	/// and then the ones without data i.e. the indices of the NestedVector must not contain any repeated indices.
	bool species_has_data(size_t species_idx) const { return species_idx < _msms_data.size(); }
	void add_to_sumlog(coretools::TSumLogProbability &sum_log, uint8_t binned_value) const {
		sum_log.add(get_likelihood_from_binned_value(binned_value));
	}
	void add_to_sumlog(std::array<coretools::TSumLogProbability, 2> &sum_log, uint8_t binned_value) const {
		sum_log[0].add(1.0 - get_likelihood_from_binned_value(binned_value));
		sum_log[1].add(get_likelihood_from_binned_value(binned_value));
	}

	void add_mass_spec_run_for_species(const std::vector<TMassSpecRun> &runs) { _msms_data.push_back(runs); }
};
