//
// The state array both dense storages are built on.
//

#pragma once

#include "constants.h"
#include "coretools/Main/TError.h"
#include "coretools/algorithms.h"
#include "storages/storage_concepts.h"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

/// One byte of state per cell of the container space, laid out in the order the linear index
/// already denotes (row-major), so a cell's linear index *is* its position in the array.
///
/// This is the whole of the internal state -- `Z` carries a state and nothing else -- and the
/// state half of the dense field, which puts a counter array of the same length beside it. Holding
/// state and counter in two arrays rather than one packed word is what lets the two dense storages
/// share this class; the sparse pair cannot, because a sparse matrix stores one value per cell and
/// that value has to carry both.
///
/// Every cell of the container space is stored, from the moment the array is sized. That is the
/// point of the dense form -- there is no "is this cell present" case to reason about, which is
/// what makes it the implementation to check the sparse one against -- and it is what three
/// members read differently for:
///   - `remove_zeros` has nothing to remove and does nothing,
///   - `fill_current_state` reports every cell of the clique as existing,
///   - `empty` answers "no cell is one", where sparse answers "nothing is stored". The two part
///     company over a cell inserted as a zero, which the sparse matrix counts as stored: after
///     `insert_zero`, sparse is not empty and this is. The caller asking (TMarkovField, checking
///     that a field it was told to hold fixed was in fact read in) means the former, and it is
///     the only reading dense has to give -- every cell is stored from the moment it is sized.
class TDenseStateArray {
private:
	IndexArray _dimensions{};
	std::vector<uint8_t> _states;

	void _throw_if_outside_container_space(size_t linear_index) const {
		if (linear_index >= _states.size()) {
			throw coretools::TDevError(
			    "You are trying to insert a value at a linear index bigger than the total size of "
			    "the container. The index is: ",
			    linear_index, " and the total size of the container is : ", _states.size());
		}
	}

public:
	TDenseStateArray() = default;
	explicit TDenseStateArray(const IndexArray &dimensions) { initialize_dimensions(dimensions); }

	/// Sizes the array to the container space and puts every cell in state 0.
	void initialize_dimensions(const IndexArray &dimensions) {
		_dimensions = dimensions;
		_states.assign(coretools::containerProduct(dimensions), 0);
	}

	[[nodiscard]] bool is_one(size_t linear_index) const {
		DEBUG_ASSERT(linear_index < _states.size());
		return _states[linear_index] != 0;
	}

	[[nodiscard]] bool is_one(const IndexArray &multidim_index) const {
		const size_t linear_index = get_linear_index_in_container_space(multidim_index);
		DEBUG_ASSERT(linear_index < _states.size());
		return _states[linear_index] != 0;
	}

	void set_state(size_t linear_index, bool state) {
		DEBUG_ASSERT(linear_index < _states.size());
		_states[linear_index] = state ? 1 : 0;
	}

	void insert_one(size_t linear_index) {
		DEBUG_ASSERT(linear_index < _states.size());
		_states[linear_index] = 1;
	}

	void insert_one(const IndexArray &multidim_index) {
		insert_one(get_linear_index_in_container_space(multidim_index));
	}

	void insert_zero(size_t linear_index) {
		DEBUG_ASSERT(linear_index < _states.size());
		_states[linear_index] = 0;
	}

	/// Nothing to remove: the zeros are the array. Sparse drops its zero cells to reclaim the
	/// memory they take, which is invisible through `is_one` -- an absent cell reads as 0 either
	/// way -- so leaving them in place keeps the two implementations in step.
	void remove_zeros() {}

	[[nodiscard]] size_t total_size_of_container_space() const { return _states.size(); }

	/// The size of each dimension of the container space. Two storages are cell-by-cell comparable
	/// (their linear indices denote the same cell) only if these are equal.
	[[nodiscard]] const IndexArray &dimensions() const { return _dimensions; }

	[[nodiscard]] bool empty() const {
		return std::none_of(_states.begin(), _states.end(), [](uint8_t s) { return s != 0; });
	}

	[[nodiscard]] size_t
	get_linear_index_in_container_space(const IndexArray &multidim_index) const {
		return coretools::getLinearIndex(multidim_index, _dimensions);
	}

	[[nodiscard]] IndexArray get_multi_dimensional_index(size_t linear_index) const {
		return coretools::getSubscriptsAsArray(linear_index, _dimensions);
	}
};

static_assert(BinaryFieldStorage<TDenseStateArray>,
              "The dense state array must satisfy the binary storage interface.");
static_assert(!FieldStorage<TDenseStateArray>,
              "A state array carries no posterior counter, so it is not a field.");
