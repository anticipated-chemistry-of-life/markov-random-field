//
// The dense field.
//

#pragma once

#include "constants.h"
#include "coretools/Main/TError.h"
#include "storages/TDenseStateArray.h"
#include "storages/storage_concepts.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

/// The field, dense: a TDenseStateArray for the state, and a counter of the same length beside it.
/// The two arrays are indexed by the same linear index, so a cell's state and its count are one
/// index computation apart without either having to make room for the other.
///
/// That is what buys the counter its sixteenth bit. The sparse field packs a cell into a single
/// 16-bit word of a sparse matrix, where the state bit takes one of them and the counter is left
/// with 15; here nothing shares the word, so a count runs to 65535 rather than 32767. Twice the
/// range is half the thinning: a chain of n iterations is sampled one iteration in
/// ceil(n / 65535), where the sparse field needs one in ceil(n / 32767).
class TStorageYDense {
private:
	TDenseStateArray _states;
	std::vector<uint16_t> _counts;
	size_t _thinning_factor = 1;
	/// The number of iterations actually counted, and so the largest a cell's counter can be.
	/// Counted rather than derived from the chain length -- see FieldStorage in
	/// storages/storage_concepts.h for why that arithmetic cannot be done up front.
	size_t _total_counts    = 0;

public:
	/// The largest value a counter can hold. A chain thinned by `get_thinning_factor()` cannot
	/// reach past it.
	static constexpr uint16_t MAX_COUNTER = std::numeric_limits<uint16_t>::max();

	TStorageYDense() = default;
	TStorageYDense(size_t n_iterations, const IndexArray &dimensions) {
		initialize(n_iterations, dimensions);
	}

	void initialize(size_t n_iterations, const std::vector<size_t> &dimensions) {
		if (dimensions.size() != NUMBER_OF_TREES) {
			throw coretools::TDevError("dimensions must have size NUMBER_OF_TREES");
		}
		initialize(n_iterations, IndexArray{dimensions[0], dimensions[1]});
	}

	void initialize(size_t n_iterations, const IndexArray &dimensions) {
		// add_to_counter takes the iteration modulo the thinning factor, so it stays at least one
		// even for a chain with no iterations to thin.
		_thinning_factor =
		    std::max<size_t>(1, static_cast<size_t>(std::ceil(static_cast<double>(n_iterations) /
		                                                      static_cast<double>(MAX_COUNTER))));
		_total_counts = 0;
		_states.initialize_dimensions(dimensions);
		_counts.assign(_states.total_size_of_container_space(), 0);
	}

	// -- State. The array answers all of it, except that inserting a cell also starts its counter
	// -- over: in the sparse field an insert writes a whole new entry, counter included, and the
	// -- two implementations have to agree on what a cell holds afterwards.

	/// See TDenseStateArray::is_stored: the dense field holds the whole container space.
	[[nodiscard]] bool is_stored(size_t linear_index) const {
		return _states.is_stored(linear_index);
	}

	[[nodiscard]] bool is_one(size_t linear_index) const { return _states.is_one(linear_index); }
	[[nodiscard]] bool is_one(const IndexArray &multidim_index) const {
		return _states.is_one(_states.get_linear_index_in_container_space(multidim_index));
	}
	void set_state(size_t linear_index, bool state) { _states.set_state(linear_index, state); }

	void insert_one(size_t linear_index) {
		_states.insert_one(linear_index);
		_counts[linear_index] = 0;
	}

	void insert_zero(size_t linear_index) {
		_states.insert_zero(linear_index);
		_counts[linear_index] = 0;
	}

	/// The cells stay -- see TDenseStateArray::remove_zeros -- but their counters do not. Sparse
	/// erases a zero cell outright, counter and all, so a later get_fraction_of_ones reads 0 for
	/// it; zeroing the counter here is what makes the two answer the same.
	void remove_zeros() {
		for (size_t i = 0; i < _counts.size(); ++i) {
			if (!_states.is_one(i)) { _counts[i] = 0; }
		}
	}

	[[nodiscard]] size_t total_size_of_container_space() const {
		return _states.total_size_of_container_space();
	}
	[[nodiscard]] const IndexArray &dimensions() const { return _states.dimensions(); }
	[[nodiscard]] bool empty() const { return _states.empty(); }

	[[nodiscard]] size_t number_of_ones() const {
		size_t count = 0;
		for (size_t i = 0; i < _states.total_size_of_container_space(); ++i) {
			if (_states.is_one(i)) { ++count; }
		}
		return count;
	}

	[[nodiscard]] size_t
	get_linear_index_in_container_space(const IndexArray &multidim_index) const {
		return _states.get_linear_index_in_container_space(multidim_index);
	}
	/// The same conversion under the name the sparse field gives it, because that is the name
	/// production reaches for it by (msms_data.cpp).
	[[nodiscard]] size_t get_linear_index_in_Y_space(const IndexArray &multidim_index) const {
		return get_linear_index_in_container_space(multidim_index);
	}
	[[nodiscard]] IndexArray get_multi_dimensional_index(size_t linear_index) const {
		return _states.get_multi_dimensional_index(linear_index);
	}

	void fill_current_state(const IndexArray &start_index, size_t K, size_t increment,
	                        std::vector<uint8_t> &current_state, std::vector<uint8_t> &exists,
	                        std::vector<size_t> &linear_index) const {
		_states.fill_current_state(start_index, K, increment, current_state, exists, linear_index);
	}

	// -- The posterior counter.

	/// Counts every cell that is currently a one, on one iteration in `get_thinning_factor()`.
	void add_to_counter(size_t iteration) {
		if (iteration % _thinning_factor != 0) { return; }
		for (size_t i = 0; i < _counts.size(); ++i) {
			if (!_states.is_one(i)) { continue; }
			if (_counts[i] == MAX_COUNTER) {
				throw coretools::TDevError("counter exceeds 16-bit maximum (", MAX_COUNTER, ")");
			}
			++_counts[i];
		}
		// after the loop, so a counter that overflows does not leave the denominator counting an
		// iteration the cells never got
		++_total_counts;
	}

	void reset_counts() {
		std::fill(_counts.begin(), _counts.end(), 0);
		_total_counts = 0;
	}

	[[nodiscard]] double get_fraction_of_ones(size_t linear_index) const {
		DEBUG_ASSERT(linear_index < _counts.size());
		// Nothing counted yet: no posterior to report, and nothing to divide by.
		if (_total_counts == 0) { return 0.0; }
		return static_cast<double>(_counts[linear_index]) / static_cast<double>(_total_counts);
	}

	/// How often the cell was counted a one -- the numerator of the fraction above. This is where
	/// the sparse field's `operator[](i).get_counter()` lands once the counter no longer shares a
	/// word with the state, and it is what a test can assert on exactly.
	[[nodiscard]] uint16_t get_counter(size_t linear_index) const {
		DEBUG_ASSERT(linear_index < _counts.size());
		return _counts[linear_index];
	}

	[[nodiscard]] size_t get_total_counts() const { return _total_counts; }
	[[nodiscard]] size_t get_thinning_factor() const { return _thinning_factor; }

	// -- The bulk paths. Deliberately outside the storage concept (see storage_concepts.h): the
	// -- sampler needs them, but they describe what the model does with a field rather than what
	// -- makes a field a field, and both implementations still spell them with a `Y` in the name.

	/// Bulk-insert deferred 0 -> 1 transitions. Mirror of TStorageYMatrix::insert_in_Y.
	///
	/// Every index goes through `insert_one`, counter and all: the sparse form writes a whole new
	/// entry per index, which starts that cell's counter over, and the two have to leave a cell
	/// holding the same thing. The update only defers a cell it found absent, whose counter is 0
	/// either way.
	void insert_in_Y(const std::vector<std::vector<size_t>> &linear_indices_to_insert) {
		insert_ones_in_batches(*this, linear_indices_to_insert);
	}

	/// The state of every cell of the container space, in ascending linear-index order.
	[[nodiscard]] std::vector<uint8_t> get_full_Y_binary_vector() const {
		return whole_space_states<uint8_t>(*this);
	}

	/// One cell as the stored-entry walks below report it. The sparse field hands back its packed
	/// entry type, whose counter is 15 bits; this one is not that type precisely because the dense
	/// counter is 16, and returning a `TStorageY` would mean either dropping the top bit or
	/// throwing on a count the dense field is entitled to hold.
	class TCell {
		bool _state       = false;
		uint16_t _counter = 0;

	public:
		TCell() = default;
		TCell(bool state, uint16_t counter) : _state(state), _counter(counter) {}

		[[nodiscard]] bool is_one() const { return _state; }
		[[nodiscard]] uint16_t get_counter() const { return _counter; }
	};

	/// Every stored cell as (linear index, value), in ascending linear-index order -- which here
	/// is every cell of the container space, since the dense field stores all of them.
	[[nodiscard]] std::vector<std::pair<size_t, TCell>> get_stored_entries() const {
		const size_t total = total_size_of_container_space();
		std::vector<std::pair<size_t, TCell>> entries;
		entries.reserve(total);
		for (size_t i = 0; i < total; ++i) {
			entries.emplace_back(i, TCell(_states.is_one(i), _counts[i]));
		}
		return entries;
	}

	/// Allocation-free forward walk over the cells that are *one*, in ascending linear-index order,
	/// with the shape TStorageYMatrix::OnesCursor has, so that the merge-joins written against the
	/// sparse field (TLotus, the simple error model) read this one unchanged.
	///
	/// The ones and not the stored cells: this field stores the whole container space, so a cursor
	/// over what it stores would visit every cell and hand the caller's closed-form term a
	/// different share of the same sum than the sparse field does -- see the note on OnesCursor.
	class OnesCursor {
		const TStorageYDense *_field = nullptr;
		size_t _index                = 0;

		void _advance_to_next_one() {
			const size_t total = _field->total_size_of_container_space();
			while (_index < total && !_field->is_one(_index)) { ++_index; }
		}

	public:
		OnesCursor() = default;
		explicit OnesCursor(const TStorageYDense &field) : _field(&field) {
			_advance_to_next_one();
		}

		[[nodiscard]] bool valid() const {
			return _field != nullptr && _index < _field->total_size_of_container_space();
		}
		[[nodiscard]] size_t linear_index() const { return _index; }
		void advance() {
			++_index;
			_advance_to_next_one();
		}
	};

	[[nodiscard]] OnesCursor ones_cursor() const { return OnesCursor(*this); }
};

static_assert(FieldStorage<TStorageYDense>,
              "The dense field must satisfy the field storage interface.");
