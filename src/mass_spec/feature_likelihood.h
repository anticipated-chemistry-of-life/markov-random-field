#pragma once

#include "coretools/Main/TError.h"

#include <cstdint>

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
	void _set_molecule_index(const uint32_t linear_index) {
		if (linear_index > _molecule_index_mask) {
			throw coretools::TUserError("Molecule index '", linear_index,
			                            "' exceeds the maximum allowed number of molecules : ", _molecule_index_mask,
			                            ".");
		}
		_value = (_value & ~_molecule_index_mask) | (linear_index & _molecule_index_mask);
	}
	void _set_binned_likelihood(const uint8_t binned_likelihood) {
		_value = (_value & ~_likelihood_mask) | (static_cast<uint32_t>(binned_likelihood) << 24);
	}

public:
	TFeatureLikelihood() = default;
	explicit TFeatureLikelihood(const uint32_t linear_index, const uint8_t binned_likelihood) {
		this->_set_molecule_index(linear_index);
		this->_set_binned_likelihood(binned_likelihood);
	}
	[[nodiscard]] uint32_t get_molecule_index() const { return _value & _molecule_index_mask; };

	[[nodiscard]] uint32_t get_binned_likelihood() const { return (_value & _likelihood_mask) >> 24; }

	bool operator<(const TFeatureLikelihood &right) const { return get_molecule_index() < right.get_molecule_index(); }
	bool operator<(const uint32_t right) const { return get_molecule_index() < right; }
	bool operator!=(const uint32_t right) const { return get_molecule_index() != right; }
	bool operator==(const uint32_t right) const { return get_molecule_index() == right; }
};
