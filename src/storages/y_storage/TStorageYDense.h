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
	size_t _total_counts    = 0;

public:
	/// The largest value a counter can hold. A chain thinned by `get_thinning_factor()` cannot
	/// reach past it.
	static constexpr uint16_t MAX_COUNTER = std::numeric_limits<uint16_t>::max();

	TStorageYDense() = default;
	TStorageYDense(size_t n_iterations, const IndexArray &dimensions) {
		initialize(n_iterations, dimensions);
	}

	void initialize(size_t n_iterations, const IndexArray &dimensions) {
		// add_to_counter takes the iteration modulo the thinning factor, so it stays at least one
		// even for a chain with no iterations to thin.
		_thinning_factor =
		    std::max<size_t>(1, static_cast<size_t>(std::ceil(static_cast<double>(n_iterations) /
		                                                      static_cast<double>(MAX_COUNTER))));
		_total_counts = n_iterations / _thinning_factor;
		_states.initialize_dimensions(dimensions);
		_counts.assign(_states.total_size_of_container_space(), 0);
	}

	// -- State. The array answers all of it, except that inserting a cell also starts its counter
	// -- over: in the sparse field an insert writes a whole new entry, counter included, and the
	// -- two implementations have to agree on what a cell holds afterwards.

	[[nodiscard]] bool is_one(size_t linear_index) const { return _states.is_one(linear_index); }
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
	[[nodiscard]] bool empty() const { return _states.empty(); }

	[[nodiscard]] size_t
	get_linear_index_in_container_space(const IndexArray &multidim_index) const {
		return _states.get_linear_index_in_container_space(multidim_index);
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
	}

	void reset_counts() { std::fill(_counts.begin(), _counts.end(), 0); }

	[[nodiscard]] double get_fraction_of_ones(size_t linear_index) const {
		DEBUG_ASSERT(linear_index < _counts.size());
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
};

static_assert(FieldStorage<TStorageYDense>,
              "The dense field must satisfy the field storage interface.");
