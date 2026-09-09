//
// The array both dense storages are built on.
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

/// A cell for every position of the container space, laid out in the order the linear index
/// already denotes (row-major), so a cell's linear index *is* its position in the array.
///
/// The cell is the template argument, because the two dense storages carry different amounts in
/// one: a node state carries a bare state, and the field carries a state and the counter beside it
/// in one packed word. Everything below is written against `cell_is_one` and `write_state`
/// (storages/cell_handle.h), so neither body knows which of the two it is holding.
///
/// Every cell of the container space is stored, from the moment the array is sized. That is the
/// point of the dense form -- there is no "is this cell present" case to reason about, which is
/// what makes it the implementation to check the sparse one against -- and it is what two members
/// read differently for:
///   - `remove_zeros` drops what a cell that is not a one carries rather than the cell itself,
///   - `empty` answers "no cell is one", where sparse answers "nothing is stored". The two part
///     company over a cell inserted as a zero, which the sparse form counts as stored: after
///     `insert_zero`, sparse is not empty and this is. The caller asking (TMarkovField, checking
///     that a field it was told to hold fixed was in fact read in) means the former, and it is
///     the only reading dense has to give -- every cell is stored from the moment it is sized.
template<typename Cell> class TDenseCellArray {
private:
	IndexArray _dimensions{};
	std::vector<Cell> _states;

	void _throw_if_outside_container_space(size_t linear_index) const {
		if (linear_index >= _states.size()) {
			throw coretools::TDevError(
			    "You are trying to insert a value at a linear index bigger than the total size of "
			    "the container. The index is: ",
			    linear_index, " and the total size of the container is : ", _states.size());
		}
	}

protected:
	/// The cells themselves, for a storage built on this one that carries more than a state in
	/// them: the dense field counts every cell that is a one, which is a pass this array has no
	/// reason to know about.
	[[nodiscard]] std::vector<Cell> &cells() { return _states; }
	[[nodiscard]] const std::vector<Cell> &cells() const { return _states; }

public:
	TDenseCellArray() = default;
	explicit TDenseCellArray(const IndexArray &dimensions) { initialize_dimensions(dimensions); }

	/// Sizes the array to the container space and puts every cell in state 0.
	void initialize_dimensions(const IndexArray &dimensions) {
		_dimensions = dimensions;
		_states.assign(coretools::containerProduct(dimensions), Cell());
	}

	[[nodiscard]] bool is_one(size_t linear_index) const {
		DEBUG_ASSERT(linear_index < _states.size());
		return cell_is_one(_states[linear_index]);
	}

	[[nodiscard]] bool is_one(const IndexArray &multidim_index) const {
		return is_one(get_linear_index_in_container_space(multidim_index));
	}

	/// What a handle from this array points at: the whole cell, state and all.
	using TCell = Cell;

	/// Where a cell is, for a caller that is about to write it. A dense array holds every cell of
	/// its container space from the moment it is sized, so the handle is always in the container
	/// and its pointer is never null.
	[[nodiscard]] IsOneResult<TCell> locate(size_t linear_index) {
		DEBUG_ASSERT(linear_index < _states.size());
		TCell *cell = &_states[linear_index];
		return {cell_is_one(*cell), true, linear_index, cell};
	}

	[[nodiscard]] IsOneResult<TCell> locate(const IndexArray &multidim_index) {
		return locate(get_linear_index_in_container_space(multidim_index));
	}

	/// Writes the state of a cell and leaves whatever else it carries as it was.
	void set_state(size_t linear_index, bool state) {
		DEBUG_ASSERT(linear_index < _states.size());
		write_state(_states[linear_index], state);
	}

	/// Writes the cell anew, so a counted cell starts its counter over. That is what an insert
	/// means in the sparse form, where it writes a whole new entry, and the two implementations
	/// have to agree on what a cell holds afterwards.
	void insert_one(size_t linear_index) {
		_throw_if_outside_container_space(linear_index);
		_states[linear_index] = Cell(true);
	}

	void insert_one(const IndexArray &multidim_index) {
		insert_one(get_linear_index_in_container_space(multidim_index));
	}

	void insert_zero(size_t linear_index) {
		_throw_if_outside_container_space(linear_index);
		_states[linear_index] = Cell(false);
	}

	/// Drops what a cell that is not a one carries, and keeps the cell. Sparse erases such a cell
	/// outright to reclaim the memory it takes, which is invisible through `is_one` -- an absent
	/// cell reads as 0 either way -- but not invisible through a counter: a cell the sparse form
	/// erased comes back with no count, so a cell kept here must lose its count too.
	///
	/// A bare state cell carries nothing besides the state, so for a node state this writes a zero
	/// over a zero and changes nothing.
	void remove_zeros() {
		for (auto &cell : _states) {
			if (!cell_is_one(cell)) { cell = Cell(); }
		}
	}

	[[nodiscard]] size_t total_size_of_container_space() const { return _states.size(); }

	/// The size of each dimension of the container space. Two storages are cell-by-cell comparable
	/// (their linear indices denote the same cell) only if these are equal.
	[[nodiscard]] const IndexArray &dimensions() const { return _dimensions; }

	[[nodiscard]] bool empty() const {
		return std::none_of(_states.begin(), _states.end(),
		                    [](const Cell &cell) { return cell_is_one(cell); });
	}

	[[nodiscard]] size_t
	get_linear_index_in_container_space(const IndexArray &multidim_index) const {
		return coretools::getLinearIndex(multidim_index, _dimensions);
	}

	[[nodiscard]] IndexArray get_multi_dimensional_index(size_t linear_index) const {
		return coretools::getSubscriptsAsArray(linear_index, _dimensions);
	}

	[[nodiscard]] size_t number_of_ones() const {
		return static_cast<size_t>(std::count_if(
		    _states.begin(), _states.end(), [](const Cell &cell) { return cell_is_one(cell); }));
	}

	/// Allocation-free forward walk over the cells that are *one*, in ascending linear-index order.
	/// The sparse array's cursor has the same shape, so a merge join reads either unchanged.
	///
	/// The ones, and not the stored cells. Which cells a storage holds is a property of the
	/// backend, and a sum split by that reaches a Metropolis ratio. tests/backend_parity/README.md
	/// records what that cost to find out.
	class OnesCursor {
		const TDenseCellArray *_array = nullptr;
		size_t _index                 = 0;

		void _advance_to_next_one() {
			const size_t total = _array->total_size_of_container_space();
			while (_index < total && !_array->is_one(_index)) { ++_index; }
		}

	public:
		OnesCursor() = default;
		explicit OnesCursor(const TDenseCellArray &array) : _array(&array) {
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

/// The dense array of bare states: one byte per cell of the container space, and nothing else.
/// This is the whole of a node state -- `Z` carries a state and nothing more -- and the storage
/// the simple error model data takes on the dense side.
using TDenseStateArray = TDenseCellArray<uint8_t>;

static_assert(BinaryStorage<TDenseStateArray>,
              "The dense state array must satisfy the binary storage interface.");
static_assert(!FieldStorage<TDenseStateArray>,
              "A state array carries no posterior counter, so it is not a field.");
