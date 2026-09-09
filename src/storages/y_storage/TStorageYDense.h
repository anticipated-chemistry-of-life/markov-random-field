//
// The dense field.
//

#pragma once

#include "TStorageY.h"
#include "constants.h"
#include "coretools/Main/TError.h"
#include "storages/TDense.h"
#include "storages/bulk_paths.h"
#include "storages/storage_concepts.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

/// The field, dense: the dense cell array over the packed field cell, and the posterior counter
/// that cell carries beside its state.
///
/// The cell is the sparse field's cell. Both fields therefore hold a 15-bit counter, and a chain
/// of n iterations is sampled one iteration in ceil(n / 32767) whichever backend runs it. That is
/// what makes a posterior field written by one comparable with a posterior field written by the
/// other: the resolution is the cell's, and no longer the backend's.
class TStorageYDense : public TDenseStorage<TStorageY> {
private:
	size_t _thinning_factor = 1;
	/// The number of iterations actually counted, and so the largest a cell's counter can be.
	/// Counted rather than derived from the chain length -- see FieldStorage in
	/// storages/storage_concepts.h for why that arithmetic cannot be done up front.
	size_t _total_counts    = 0;

public:
	/// The largest value a counter can hold. A chain thinned by `get_thinning_factor()` cannot
	/// reach past it.
	static constexpr uint16_t MAX_COUNTER = TStorageY::MAX_COUNTER;

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
		initialize_dimensions(dimensions);
	}

	// -- State. The array answers all of it. An insert writes the cell anew, counter included,
	// -- which is what the sparse field's insert does too.

	/// The same conversion as `get_linear_index_in_container_space`, under the name the sparse
	/// field gives it, because that is the name production reaches for it by (msms_data.cpp).
	[[nodiscard]] size_t get_linear_index_in_Y_space(const IndexArray &multidim_index) const {
		return get_linear_index_in_container_space(multidim_index);
	}

	// -- The posterior counter.

	/// Counts every cell that is currently a one, on one iteration in `get_thinning_factor()`.
	void add_to_counter(size_t iteration) {
		if (iteration % _thinning_factor != 0) { return; }
		for (auto &cell : cells()) { cell.update_counter(); }
		// after the loop, so a counter that overflows does not leave the denominator counting an
		// iteration the cells never got
		++_total_counts;
	}

	void reset_counts() {
		for (auto &cell : cells()) { cell.reset_counter(); }
		_total_counts = 0;
	}

	[[nodiscard]] double get_fraction_of_ones(size_t linear_index) const {
		DEBUG_ASSERT(linear_index < total_size_of_container_space());
		// Nothing counted yet: no posterior to report, and nothing to divide by.
		if (_total_counts == 0) { return 0.0; }
		return static_cast<double>(cells()[linear_index].get_counter()) /
		       static_cast<double>(_total_counts);
	}

	/// How often the cell was counted a one -- the numerator of the fraction above, and what a
	/// test can assert on exactly.
	[[nodiscard]] uint16_t get_counter(size_t linear_index) const {
		DEBUG_ASSERT(linear_index < total_size_of_container_space());
		return cells()[linear_index].get_counter();
	}

	[[nodiscard]] size_t get_total_counts() const { return _total_counts; }
	[[nodiscard]] size_t get_thinning_factor() const { return _thinning_factor; }

	// -- The bulk paths. Deliberately outside the storage concept (see storage_concepts.h): the
	// -- sampler needs them, but they describe what the model does with a field rather than what
	// -- makes a field a field, and both implementations still spell them with a `Y` in the name.

	/// Bulk-insert deferred 0 -> 1 transitions. Mirror of TStorageYSparse::insert_in_Y.
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

	/// Every stored cell as (linear index, cell), in ascending linear-index order -- which here is
	/// every cell of the container space, since the dense field stores all of them.
	[[nodiscard]] std::vector<std::pair<size_t, TStorageY>> get_stored_entries() const {
		const size_t total = total_size_of_container_space();
		std::vector<std::pair<size_t, TStorageY>> entries;
		entries.reserve(total);
		for (size_t i = 0; i < total; ++i) { entries.emplace_back(i, cells()[i]); }
		return entries;
	}

	/// One stored cell, by linear index. The sparse field spells the same lookup, so a test that
	/// reads a counter reads either field through it.
	[[nodiscard]] TStorageY operator[](size_t linear_index) const {
		DEBUG_ASSERT(linear_index < total_size_of_container_space());
		return cells()[linear_index];
	}
};

static_assert(FieldStorage<TStorageYDense>,
              "The dense field must satisfy the field storage interface.");
