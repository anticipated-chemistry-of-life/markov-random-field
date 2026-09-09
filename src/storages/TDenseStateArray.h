//
// The state array both dense storages are built on.
//

#pragma once

#include "constants.h"
#include "coretools/Main/TError.h"
#include "coretools/algorithms.h"
#include "storages/cell_handle.h"
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
/// what makes it the implementation to check the sparse one against -- and it is what two members
/// read differently for:
///   - `remove_zeros` has nothing to remove and does nothing,
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

	/// One byte per cell, which is what a handle from this array points at.
	using TCell = uint8_t;

	/// Where a cell is, for a caller that is about to write it. A dense array holds every cell of
	/// its container space from the moment it is sized, so the handle is always in the container
	/// and its pointer is never null.
	[[nodiscard]] IsOneResult<TCell> locate(size_t linear_index) {
		DEBUG_ASSERT(linear_index < _states.size());
		TCell *cell = &_states[linear_index];
		return {*cell != 0, true, linear_index, cell};
	}

	[[nodiscard]] IsOneResult<TCell> locate(const IndexArray &multidim_index) {
		return locate(get_linear_index_in_container_space(multidim_index));
	}

	void set_state(size_t linear_index, bool state) {
		DEBUG_ASSERT(linear_index < _states.size());
		_states[linear_index] = state ? 1 : 0;
	}

	void insert_one(size_t linear_index) {
		_throw_if_outside_container_space(linear_index);
		_states[linear_index] = 1;
	}

	void insert_one(const IndexArray &multidim_index) {
		insert_one(get_linear_index_in_container_space(multidim_index));
	}

	void insert_zero(size_t linear_index) {
		_throw_if_outside_container_space(linear_index);
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

	[[nodiscard]] size_t number_of_ones() const {
		return static_cast<size_t>(
		    std::count_if(_states.begin(), _states.end(), [](uint8_t s) { return s != 0; }));
	}

	/// Allocation-free forward walk over the cells that are *one*, in ascending linear-index order.
	/// TSparseBinaryArray::OnesCursor has the same shape, so a merge join reads either unchanged.
	///
	/// The ones, and not the stored cells. Which cells a storage holds is a property of the
	/// backend, and a sum split by that reaches a Metropolis ratio. tests/backend_parity/README.md
	/// records what that cost to find out.
	class OnesCursor {
		const TDenseStateArray *_array = nullptr;
		size_t _index                  = 0;

		void _advance_to_next_one() {
			const size_t total = _array->total_size_of_container_space();
			while (_index < total && !_array->is_one(_index)) { ++_index; }
		}

	public:
		OnesCursor() = default;
		explicit OnesCursor(const TDenseStateArray &array) : _array(&array) {
			_advance_to_next_one();
		}

		[[nodiscard]] bool valid() const {
			return _array != nullptr && _index < _array->total_size_of_container_space();
		}
		[[nodiscard]] size_t linear_index() const { return _index; }
		void advance() {
			++_index;
			_advance_to_next_one();
		}
	};

	[[nodiscard]] OnesCursor ones_cursor() const { return OnesCursor(*this); }
};

static_assert(BinaryStorage<TDenseStateArray>,
              "The dense state array must satisfy the binary storage interface.");
static_assert(LocatableStorage<TDenseStateArray>,
              "The dense state array must point an updater at one of its cells.");
static_assert(!FieldStorage<TDenseStateArray>,
              "A state array carries no posterior counter, so it is not a field.");

// -------------------------------------------------------------------------------------------
// The two bulk paths the dense field and the dense internal state have in common.
//
// Each storage has to answer to a name of its own -- `insert_in_Y` beside `insert_in_Z`,
// `get_full_Y_binary_vector` beside `get_full_Z_binary_vector`, because that is what the sparse
// implementations are called and what the sampler asks for -- so the shape they share lives here
// rather than in a member either could inherit.
// -------------------------------------------------------------------------------------------

/// Writes a one at every index of every batch. The sparse implementations merge the per-thread
/// batches the update hands over and re-sort the matrix once; a dense array has nothing to sort, so
/// this is the writes and nothing else.
///
/// It goes through the storage rather than through its state array because the field's insert does
/// something the array's does not: it starts the cell's counter over, which is what the sparse
/// form's whole-new-entry write also does.
template<BinaryStorage Storage>
void insert_ones_in_batches(Storage &storage, const std::vector<std::vector<size_t>> &batches) {
	for (const auto &batch : batches) {
		for (const size_t linear_index : batch) { storage.insert_one(linear_index); }
	}
}

/// The state of every cell of the container space, in ascending linear-index order.
///
/// The element type is the caller's because the two sparse implementations disagree on it -- the
/// field dumps a vector of bytes and the internal state a vector of words -- and a trace line is
/// written from whatever they return.
template<typename T, BinaryStorage Storage>
[[nodiscard]] std::vector<T> whole_space_states(const Storage &storage) {
	const size_t total = storage.total_size_of_container_space();
	std::vector<T> states;
	states.reserve(total);
	for (size_t i = 0; i < total; ++i) { states.push_back(static_cast<T>(storage.is_one(i))); }
	return states;
}
