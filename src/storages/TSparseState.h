#pragma once

#include "constants.h"
#include "storages/storage_concepts.h"
#include <coretools/Main/TError.h>
#include <coretools/algorithms.h>
#include <cstddef>
#include <unordered_map>
#include <vector>

class TSparseBinaryArray {
private:
	IndexArray _dimensions{};
	std::unordered_map<size_t, bool> _states{};

public:
	TSparseBinaryArray() = default;
	using iterator       = std::unordered_map<size_t, bool>::iterator;
	using const_iterator = std::unordered_map<size_t, bool>::const_iterator;
	explicit TSparseBinaryArray(const IndexArray &dimensions) { initialize(dimensions); }
	/// Sizes the array to the container space and puts every cell in state 0.
	void initialize(const IndexArray &dimensions) { _dimensions = dimensions; }
	void initialize(const std::vector<size_t> &dimensions) {
		// the vector must be of the same size as the index array else we have a problem
		if (dimensions.size() != _dimensions.size()) {
			throw std::invalid_argument(
			    "dimensions vector must have the same size as the index array");
		}
		std::copy(dimensions.begin(), dimensions.end(), _dimensions.begin());
	}

	[[nodiscard]] size_t
	get_linear_index_in_container_space(const IndexArray &multidim_index) const {
		return coretools::getLinearIndex(multidim_index, _dimensions);
	}

	[[nodiscard]] IndexArray get_multi_dimensional_index(size_t linear_index) const {
		return coretools::getSubscriptsAsArray(linear_index, _dimensions);
	}

	[[nodiscard]] IsOneResult<const_iterator> is_one(size_t linear_index) const {
		const auto it = _states.find(linear_index);
		return {it != _states.end() && it->second, it != _states.end(), linear_index, it};
	}

	[[nodiscard]] IsOneResult<iterator> is_one(size_t linear_index) {
		auto it = _states.find(linear_index);
		return {it != _states.end() && it->second, it != _states.end(), linear_index, it};
	}

	[[nodiscard]] IsOneResult<const_iterator> is_one(const IndexArray &multidim_index) const {
		return is_one(get_linear_index_in_container_space(multidim_index));
	}

	[[nodiscard]] IsOneResult<iterator> is_one(const IndexArray &multidim_index) {
		return is_one(get_linear_index_in_container_space(multidim_index));
	}

	void set_state(size_t linear_index, bool state) { _states[linear_index] = state; }
	void insert_one(size_t linear_index) {
		DEBUG_ASSERT(_states.find(linear_index) == _states.end());
		_states[linear_index] = true;
	}
	void insert_one(const IndexArray &multidim_index) {
		insert_one(get_linear_index_in_container_space(multidim_index));
	}
	void insert_zero(size_t linear_index) { _states[linear_index] = false; }
	void remove_zeros() {
		for (auto it = _states.begin(); it != _states.end();) {
			if (!it->second) {
				it = _states.erase(it);
			} else {
				++it;
			}
		}
	}

	[[nodiscard]] const IndexArray &dimensions() const { return _dimensions; }

	[[nodiscard]] bool empty() const { return _states.empty(); }
	[[nodiscard]] size_t total_size_of_container_space() const {
		return coretools::containerProduct(_dimensions);
	}

	[[nodiscard]] size_t number_of_ones() const {
		size_t count = 0;
		for (const auto &[linear_index, state] : _states) {
			if (state) { ++count; }
		}
		return count;
	}

	void
	insert_ones_in_container(const std::vector<std::vector<size_t>> &linear_indices_to_insert) {
		for (const auto &linear_indices : linear_indices_to_insert) {
			for (size_t linear_index : linear_indices) { insert_one(linear_index); }
		}
	}
};

static_assert(BinaryFieldStorage<TSparseBinaryArray>,
              "The dense state array must satisfy the binary storage interface.");
static_assert(!FieldStorage<TSparseBinaryArray>,
              "A state array carries no posterior counter, so it is not a field.");
