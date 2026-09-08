//
// The sparse mirror of the dense state array.
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
#include <unordered_map>
#include <vector>

/// One binary state per *stored* cell, held in a hash map keyed by the cell's linear index. A cell
/// the map does not hold reads as state 0, so `is_one` is total over the container space.
///
/// This is the sparse mirror of TDenseStateArray: the same surface, over the cells a caller put in
/// rather than over all of them. The LOTUS records take it whatever the build selects for the
/// field. It is also the sparse half of the binary-storage alias.
///
/// A stored cell carries a state and not merely a key. A cell that goes to zero is then a write
/// and not an erase, and `remove_zeros` reclaims it between iterations. That is also what lets a
/// sparse node state be built on this later, where an erase inside a parallel region would
/// restructure the map under another thread.
class TSparseBinaryArray {
private:
	IndexArray _dimensions{};
	size_t _total_size = 0;
	std::unordered_map<size_t, uint8_t> _states;

	/// The cells that are one, in ascending linear-index order. This is what `ones_cursor` walks:
	/// a hash map has no order of its own, and the merge joins that read the cursor need one.
	mutable std::vector<size_t> _sorted_ones;
	/// Whether `_sorted_ones` still describes the map. Every mutator sets it. The next cursor
	/// clears it. So the observed data sorts once, a storage that changed sorts again, and a
	/// caller has nothing to remember to call.
	///
	/// `remove_zeros` sets it too, and costs a sort that finds the same ones. Every mutator
	/// setting the flag is one rule. An exemption per mutator is a rule per mutator, and the next
	/// one added gets it wrong.
	mutable bool _ones_are_stale = true;
	/// Whether this array has ever handed out a handle. `locate` sets it, and nothing clears it.
	///
	/// A write through a handle does not pass through this class, so from the first `locate` the
	/// array can no longer tell when its ones changed. The cursor then rebuilds on every call. The
	/// flag above cannot carry this: a cursor taken between the `locate` and the write would clear
	/// it, and the write would leave a cache the array believes in.
	///
	/// An immutable observation locates nothing and still sorts once, which is what the flag above
	/// is for. A storage that is written every iteration rebuilds every iteration either way.
	bool _handed_out_a_handle    = false;

	void _throw_if_outside_container_space(size_t linear_index) const {
		if (linear_index >= _total_size) {
			throw coretools::TDevError(
			    "You are trying to insert a value at a linear index bigger than the total size of "
			    "the container. The index is: ",
			    linear_index, " and the total size of the container is : ", _total_size);
		}
	}

	void _insert(size_t linear_index, bool state) {
		_throw_if_outside_container_space(linear_index);
		_states[linear_index] = state ? 1 : 0;
		_ones_are_stale       = true;
	}

	void _rebuild_sorted_ones() const {
		_sorted_ones.clear();
		_sorted_ones.reserve(_states.size());
		for (const auto &[linear_index, state] : _states) {
			if (state != 0) { _sorted_ones.push_back(linear_index); }
		}
		std::sort(_sorted_ones.begin(), _sorted_ones.end());
		_ones_are_stale = false;
	}

public:
	TSparseBinaryArray() = default;
	explicit TSparseBinaryArray(const IndexArray &dimensions) { initialize_dimensions(dimensions); }

	/// Sizes the container space and drops every stored cell, so every cell reads as state 0.
	void initialize_dimensions(const IndexArray &dimensions) {
		_dimensions = dimensions;
		_total_size = coretools::containerProduct(dimensions);
		_states.clear();
		_sorted_ones.clear();
		_ones_are_stale = true;
	}

	[[nodiscard]] bool is_one(size_t linear_index) const {
		DEBUG_ASSERT(linear_index < _total_size);
		const auto stored = _states.find(linear_index);
		return stored != _states.end() && stored->second != 0;
	}

	/// One byte per stored cell, which is what a handle from this array points at.
	using TCell = uint8_t;

	/// Where a cell is, for a caller that is about to write it. The map holds the cells it was
	/// given, so a cell it does not hold reads as state 0 and comes back with a null pointer. That
	/// is the cell `write_or_defer` defers.
	///
	/// This hands out a write the array cannot see, so from here on the ones cursor rebuilds on
	/// every call. Two threads must not call this on one array at once, because it writes that
	/// decision -- the same reason `ones_cursor` may not be called from two. The observed data,
	/// which is what this array holds, is read one cell at a time and never located.
	[[nodiscard]] IsOneResult<TCell> locate(size_t linear_index) {
		DEBUG_ASSERT(linear_index < _total_size);
		const auto stored = _states.find(linear_index);
		if (stored == _states.end()) { return {false, false, linear_index, nullptr}; }
		_handed_out_a_handle = true;
		return {stored->second != 0, true, linear_index, &stored->second};
	}

	[[nodiscard]] IsOneResult<TCell> locate(const IndexArray &multidim_index) {
		return locate(get_linear_index_in_container_space(multidim_index));
	}

	/// Writes the state of a cell, and stores the cell if the map does not hold it yet.
	void set_state(size_t linear_index, bool state) {
		DEBUG_ASSERT(linear_index < _total_size);
		_states[linear_index] = state ? 1 : 0;
		_ones_are_stale       = true;
	}

	void insert_one(size_t linear_index) { _insert(linear_index, true); }
	void insert_zero(size_t linear_index) { _insert(linear_index, false); }

	/// Drops every stored cell whose state is zero, which is invisible through `is_one` -- an
	/// absent cell reads as zero either way -- and reclaims the memory those cells take.
	void remove_zeros() {
		for (auto it = _states.begin(); it != _states.end();) {
			if (it->second == 0) {
				it = _states.erase(it);
			} else {
				++it;
			}
		}
		_ones_are_stale = true;
	}

	[[nodiscard]] size_t total_size_of_container_space() const { return _total_size; }

	/// The size of each dimension of the container space. Two storages are cell-by-cell comparable
	/// (their linear indices denote the same cell) only if these are equal.
	[[nodiscard]] const IndexArray &dimensions() const { return _dimensions; }

	/// "Nothing is stored", which is what the sparse storages answer. A cell inserted as a zero is
	/// stored, so this is not "no cell is one" -- TDenseStateArray::empty gives that reading.
	[[nodiscard]] bool empty() const { return _states.empty(); }

	[[nodiscard]] size_t number_of_ones() const {
		size_t count = 0;
		for (const auto &[linear_index, state] : _states) {
			if (state != 0) { ++count; }
		}
		return count;
	}

	[[nodiscard]] size_t
	get_linear_index_in_container_space(const IndexArray &multidim_index) const {
		return coretools::getLinearIndex(multidim_index, _dimensions);
	}

	[[nodiscard]] IndexArray get_multi_dimensional_index(size_t linear_index) const {
		return coretools::getSubscriptsAsArray(linear_index, _dimensions);
	}

	/// Allocation-free forward walk over the cells that are *one*, in ascending linear-index order,
	/// with the shape the field's cursor has, so that the merge joins written against a field read
	/// this one unchanged.
	///
	/// The cursor points into the array's sorted cache. It is valid until the next write.
	class OnesCursor {
		const std::vector<size_t> *_ones = nullptr;
		size_t _position                 = 0;

	public:
		OnesCursor() = default;
		explicit OnesCursor(const std::vector<size_t> &ones) : _ones(&ones) {}

		[[nodiscard]] bool valid() const { return _ones != nullptr && _position < _ones->size(); }
		[[nodiscard]] size_t linear_index() const { return (*_ones)[_position]; }
		void advance() { ++_position; }
	};

	/// The ones, in ascending linear-index order.
	///
	/// A write before this call makes it sort again, so it writes the cache it reads. Two threads
	/// must therefore not call it on one array at once. No caller does, because the merge joins
	/// run outside the parallel regions.
	///
	/// An array that has handed out a handle sorts again every call. It cannot see what a handle
	/// wrote, so it trusts nothing it sorted before.
	[[nodiscard]] OnesCursor ones_cursor() const {
		if (_ones_are_stale || _handed_out_a_handle) { _rebuild_sorted_ones(); }
		return OnesCursor(_sorted_ones);
	}
};

static_assert(BinaryStorage<TSparseBinaryArray>,
              "The sparse binary array must satisfy the binary storage interface.");
static_assert(LocatableStorage<TSparseBinaryArray>,
              "The sparse binary array must point an updater at one of its cells.");
static_assert(!FieldStorage<TSparseBinaryArray>,
              "A binary array carries no posterior counter, so it is not a field.");
