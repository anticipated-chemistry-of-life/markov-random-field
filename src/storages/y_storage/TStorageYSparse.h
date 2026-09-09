//
// The sparse field.
//

#pragma once

#include "TStorageY.h"
#include "constants.h"
#include "coretools/Main/TError.h"
#include "storages/TSparse.h"
#include "storages/bulk_paths.h"
#include "storages/storage_concepts.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

/// The field, sparse: the sparse cell array over the packed field cell, and the posterior counter
/// that cell carries beside its state.
///
/// Memory tracks the number of ones rather than the size of the leaf space, which is what this
/// backend exists for (ADR-0006). A cell the map does not hold reads as state 0 and carries no
/// count, so a run over a field that is mostly zeros pays for the ones alone.
///
/// The cell is the dense field's cell, so both fields hold a 15-bit counter and thin a chain
/// identically. Nothing about a posterior field's resolution follows from which backend wrote it.
class TStorageYSparse : public TSparseStorage<TStorageY> {
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

	TStorageYSparse() = default;
	TStorageYSparse(size_t n_iterations, const IndexArray &dimensions) {
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

	// -- State. The array answers all of it. An insert writes a whole new entry, counter included,
	// -- which is what the dense field's insert does too.

	/// The same conversion as `get_linear_index_in_container_space`, under the name production
	/// reaches for it by (msms_data.cpp).
	[[nodiscard]] size_t get_linear_index_in_Y_space(const IndexArray &multidim_index) const {
		return get_linear_index_in_container_space(multidim_index);
	}

	// -- The posterior counter.

	/// Counts every cell that is currently a one, on one iteration in `get_thinning_factor()`.
	void add_to_counter(size_t iteration) {
		if (iteration % _thinning_factor != 0) { return; }
		for (auto &[linear_index, cell] : cells()) { cell.update_counter(); }
		// after the loop, so a counter that overflows does not leave the denominator counting an
		// iteration the cells never got
		++_total_counts;
	}

	void reset_counts() {
		for (auto &[linear_index, cell] : cells()) { cell.reset_counter(); }
		_total_counts = 0;
	}

	[[nodiscard]] double get_fraction_of_ones(size_t linear_index) const {
		// Nothing counted yet: no posterior to report, and nothing to divide by.
		if (_total_counts == 0) { return 0.0; }
		return static_cast<double>((*this)[linear_index].get_counter()) /
		       static_cast<double>(_total_counts);
	}

	[[nodiscard]] size_t get_total_counts() const { return _total_counts; }
	[[nodiscard]] size_t get_thinning_factor() const { return _thinning_factor; }

	/// One cell, by linear index. A cell the map does not hold comes back default-constructed --
	/// state 0 and no count -- which is what its absence says.
	[[nodiscard]] TStorageY operator[](size_t linear_index) const {
		const auto stored = cells().find(linear_index);
		return stored == cells().end() ? TStorageY() : stored->second;
	}

	[[nodiscard]] TStorageY operator[](const IndexArray &multidim_index) const {
		return (*this)[get_linear_index_in_container_space(multidim_index)];
	}

	// -- The bulk paths. Deliberately outside the storage concept (see storage_concepts.h): the
	// -- sampler needs them, but they describe what the model does with a field rather than what
	// -- makes a field a field, and both implementations still spell them with a `Y` in the name.

	/// Bulk-insert deferred 0 -> 1 transitions. Mirror of TStorageYDense::insert_in_Y.
	void insert_in_Y(const std::vector<std::vector<size_t>> &linear_indices_to_insert) {
		insert_ones_in_batches(*this, linear_indices_to_insert);
	}

	/// The state of every cell of the container space, in ascending linear-index order.
	[[nodiscard]] std::vector<uint8_t> get_full_Y_binary_vector() const {
		return whole_space_states<uint8_t>(*this);
	}

	/// Every stored cell as (linear index, cell), in ascending linear-index order. Which cells
	/// those are is a property of this backend: the map holds the cells it was given, ones and
	/// zeros alike, where the dense field holds the whole container space. Only the posterior
	/// field is written from this, and it drops a cell that is neither a one nor ever counted
	/// (TMarkovField::_write_only_values_in_Y_vector), which is what makes the two files agree.
	[[nodiscard]] std::vector<std::pair<size_t, TStorageY>> get_stored_entries() const {
		return stored_cells_in_order();
	}
};

static_assert(FieldStorage<TStorageYSparse>,
              "The sparse field must satisfy the field storage interface.");
