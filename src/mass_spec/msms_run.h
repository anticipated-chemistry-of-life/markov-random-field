#pragma once

/// Class storing a mass spectrometry run. The iterator goes along the
/// molecule dimension. So for each feature, we have a list of molecule indices and binned
/// likelihoods.
#include "./feature_likelihood.h"
#include "coretools/Containers/TNestedVector.h"
#include <optional>

struct BinarySearchResult {
	bool found{};
	std::optional<uint8_t> binned_likelihood;
	size_t index{};
};

class TMassSpecRun {
private:
	coretools::TNestedVector<TFeatureLikelihood> _features;

public:
	TMassSpecRun() = default;
	[[nodiscard]] size_t capacity() const { return _features.capacity(); }
	void reserve(size_t n) { _features.reserve(n); }
	[[nodiscard]] bool empty() const { return _features.empty(); }
	void add_empty_feature() { _features.push_back(); }

	[[nodiscard]] coretools::TConstView<TFeatureLikelihood>
	get_likelihoods_for_feature(size_t i) const {
		return _features.get(i);
	}
	[[nodiscard]] size_t size() const { return _features.size(); }
	[[nodiscard]] auto begin() const { return _features.begin(); };
	[[nodiscard]] auto end() const { return _features.end(); }
	auto begin() { return _features.begin(); }
	auto end() { return _features.end(); }
	[[nodiscard]] BinarySearchResult is_molecule_in_feature(size_t feature_idx,
	                                                        uint32_t molecule_index) const {
		const auto &likelihoods = _features.get(feature_idx);
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

	/// Adds a vector of molecule likelihoods. This represents one feature with all the molecule
	/// likelihoods that are associated to it. The vector will be sorted before being added to
	/// guarantee the order and be able to do binary search.
	void add_likelihood_vector(std::vector<TFeatureLikelihood> &feature_likelihoods) {
		std::sort(feature_likelihoods.begin(), feature_likelihoods.end());
		_features.push_back(feature_likelihoods);
	}
};
