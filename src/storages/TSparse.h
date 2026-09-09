//
// The map both sparse storages are built on.
//

#pragma once

#include "constants.h"
#include "coretools/Main/TError.h"
#include "coretools/algorithms.h"
#include "storages/cell_handle.h"
#include "storages/storage_concepts.h"
#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>

/// A flag that starts false, is set once and never cleared, and which several threads may set at
/// the same moment.
///
/// It is atomic because an update writes one storage from every thread, and each of those writes
/// starts by asking the storage where its cell is. A plain `bool` written from two threads is a
/// data race however harmless the value looks. The store is skipped once the flag is set, so after
/// the first one every thread only reads a line nothing writes again, and the relaxed order is
/// enough: what the flag guards is read after the parallel region has joined.
class TStickyFlag {
	std::atomic<bool> _set = false;

public:
	TStickyFlag() = default;
	TStickyFlag(const TStickyFlag &other) : _set(other.is_set()) {}
	TStickyFlag &operator=(const TStickyFlag &other) {
		_set.store(other.is_set(), std::memory_order_relaxed);
		return *this;
	}

	void set() {
		if (!is_set()) { _set.store(true, std::memory_order_relaxed); }
	}
	[[nodiscard]] bool is_set() const { return _set.load(std::memory_order_relaxed); }
};

/// A cell for every position a caller wrote to, held in a hash map keyed by the linear index of
/// that position. A position the map does not hold reads as state 0, so `is_one` is total over the
/// container space and memory tracks the number of ones rather than the size of the space.
///
/// This is the sparse mirror of TDenseCellArray: the same surface, over the cells a caller put in
/// rather than over all of them, and over the same two cells -- a bare state for a node state and
/// the observed data, a state packed with its posterior counter for the field. Everything below is
/// written against `cell_is_one` and `write_state` (storages/cell_handle.h), so neither body knows
/// which of the two it is holding.
///
/// A stored cell carries a state and not merely a key. A cell that goes to zero is then a write
/// and not an erase, which is what lets an update inside a parallel region write it in place: an
/// erase would restructure the map under another thread. `remove_zeros` reclaims such a cell
/// between iterations, where nothing else is running.
///
/// What a stored cell costs is a map node and a bucket slot, which is some tens of bytes for a
/// cell of one or two. The dense form pays the cell alone, for every position of the space. So
/// this form wins on memory well below one one in twenty cells, and not merely below one in two.
/// ADR-0006 argues the choice on fill for that reason.
template<typename Cell> class TSparseStorage {
private:
	IndexArray _dimensions{};
	size_t _total_size = 0;
	std::unordered_map<size_t, Cell> _states;

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
	/// array can no longer tell when its ones changed. The cursor then rebuilds on every call,
	/// which is what the field does: it is written through handles every iteration. The flag above
	/// cannot carry this: a cursor taken between the `locate` and the write would clear it, and
	/// the write would leave a cache the array believes in.
	///
	/// An immutable observation locates nothing and still sorts once, which is what the flag above
	/// is for.
	TStickyFlag _handed_out_a_handle;

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
		_states[linear_index] = Cell(state);
		_ones_are_stale       = true;
	}

	void _rebuild_sorted_ones() const {
		_sorted_ones.clear();
		_sorted_ones.reserve(_states.size());
		for (const auto &[linear_index, cell] : _states) {
			if (cell_is_one(cell)) { _sorted_ones.push_back(linear_index); }
		}
		std::sort(_sorted_ones.begin(), _sorted_ones.end());
		_ones_are_stale = false;
	}

protected:
	/// The cells themselves, for a storage built on this one that carries more than a state in
	/// them: the sparse field counts every cell that is a one, which is a pass this array has no
	/// reason to know about.
	[[nodiscard]] std::unordered_map<size_t, Cell> &cells() { return _states; }
	[[nodiscard]] const std::unordered_map<size_t, Cell> &cells() const { return _states; }

	/// Every stored cell as (linear index, cell), in ascending linear-index order. The map has no
	/// order of its own, and a file written from these has to read the same under either backend.
	[[nodiscard]] std::vector<std::pair<size_t, Cell>> stored_cells_in_order() const {
		std::vector<std::pair<size_t, Cell>> entries(_states.begin(), _states.end());
		std::sort(entries.begin(), entries.end(),
		          [](const auto &left, const auto &right) { return left.first < right.first; });
		return entries;
	}

public:
	TSparseStorage() = default;
	explicit TSparseStorage(const IndexArray &dimensions) { initialize_dimensions(dimensions); }

	/// Sizes the container space and drops every stored cell, so every cell reads as state 0.
	///
	/// It clears the handle flag too. Dropping every cell is what invalidates every handle, so an
	/// array that comes out of here has none, and it starts again as an array whose ones cursor
	/// can trust its cache.
	void initialize_dimensions(const IndexArray &dimensions) {
		_dimensions = dimensions;
		_total_size = coretools::containerProduct(dimensions);
		_states.clear();
		_sorted_ones.clear();
		_ones_are_stale      = true;
		_handed_out_a_handle = TStickyFlag();
	}

	[[nodiscard]] bool is_one(size_t linear_index) const {
		DEBUG_ASSERT(linear_index < _total_size);
		const auto stored = _states.find(linear_index);
		return stored != _states.end() && cell_is_one(stored->second);
	}

	[[nodiscard]] bool is_one(const IndexArray &multidim_index) const {
		return is_one(get_linear_index_in_container_space(multidim_index));
	}

	/// What a handle from this array points at: the whole stored cell, state and all.
	using TCell = Cell;

	/// Where a cell is, for a caller that is about to write it. The map holds the cells it was
	/// given, so a cell it does not hold reads as state 0 and comes back with a null pointer. That
	/// is the cell `write_or_defer` defers, and it is the only cell an update ever has to defer:
	/// every other write lands in place, a one-to-zero transition included.
	///
	/// This hands out a write the array cannot see, so from here on the ones cursor rebuilds on
	/// every call. An update calls this from every thread at once, on cells of its own, which is
	/// why the flag that records it is atomic and the map is only ever read here.
	[[nodiscard]] IsOneResult<TCell> locate(size_t linear_index) {
		DEBUG_ASSERT(linear_index < _total_size);
		const auto stored = _states.find(linear_index);
		if (stored == _states.end()) { return {false, false, linear_index, nullptr}; }
		_handed_out_a_handle.set();
		return {cell_is_one(stored->second), true, linear_index, &stored->second};
	}

	[[nodiscard]] IsOneResult<TCell> locate(const IndexArray &multidim_index) {
		return locate(get_linear_index_in_container_space(multidim_index));
	}

	/// Writes the state of a cell, and stores the cell if the map does not hold it yet. A cell the
	/// map already holds keeps whatever else it carries, so a counted cell keeps its counter.
	void set_state(size_t linear_index, bool state) {
		DEBUG_ASSERT(linear_index < _total_size);
		const auto [stored, inserted] = _states.try_emplace(linear_index, Cell(state));
		if (!inserted) { write_state(stored->second, state); }
		_ones_are_stale = true;
	}

	void insert_one(size_t linear_index) { _insert(linear_index, true); }
	void insert_one(const IndexArray &multidim_index) {
		insert_one(get_linear_index_in_container_space(multidim_index));
	}
	void insert_zero(size_t linear_index) { _insert(linear_index, false); }

	/// Drops every stored cell that is not a one, which is invisible through `is_one` -- an absent
	/// cell reads as zero either way -- and reclaims the memory those cells take. A counted cell
	/// goes with its counter, which is why the dense form clears the counter of a cell it keeps.
	void remove_zeros() {
		for (auto it = _states.begin(); it != _states.end();) {
			if (!cell_is_one(it->second)) {
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
	/// stored, so this is not "no cell is one" -- TDenseCellArray::empty gives that reading.
	[[nodiscard]] bool empty() const { return _states.empty(); }

	/// How many cells the map holds, ones and zeros alike. A property of this backend and not of
	/// what it holds, so nothing that reaches an arithmetic result may read it.
	[[nodiscard]] size_t size() const { return _states.size(); }

	[[nodiscard]] size_t number_of_ones() const {
		size_t count = 0;
		for (const auto &[linear_index, cell] : _states) {
			if (cell_is_one(cell)) { ++count; }
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
	/// with the shape the dense cursor has, so that the merge joins written against one read the
	/// other unchanged.
	///
	/// The cursor points into the array's sorted cache. It is valid until the next write, and
	/// until the next cursor: a second `ones_cursor()` may rebuild the very vector the first one
	/// walks. So one cursor at a time per array. A merge join takes one cursor from each of two
	/// storages, which is what every caller does.
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
		if (_ones_are_stale || _handed_out_a_handle.is_set()) { _rebuild_sorted_ones(); }
		return OnesCursor(_sorted_ones);
	}
};

/// The sparse array of bare states: one byte per stored cell, and nothing else. The LOTUS records
/// take it whatever the build selects for the field, and it is the sparse half of the
/// binary-storage alias.
using TSparseBinary = TSparseStorage<uint8_t>;

static_assert(BinaryStorage<TSparseBinary>,
              "The sparse binary array must satisfy the binary storage interface.");
static_assert(!FieldStorage<TSparseBinary>,
              "A binary array carries no posterior counter, so it is not a field.");
