//
// The dense internal state.
//

#pragma once

#include "TStorageZ.h"
#include "storages/TDenseStateArray.h"

#include <cstddef>
#include <utility>
#include <vector>

/// The dense internal state: the dense state array, plus the three bulk paths `Z` is asked for by
/// name. The state array is inherited rather than wrapped because `Z` carries one binary state per
/// cell and nothing else -- there is no member to add, only the whole-space dump, the bulk insert
/// and the stored-entry walk, none of which the storage concept covers (see storage_concepts.h) and
/// all of which the sparse implementation spells with a `Z` in the name too.
///
/// Every cell of the container space is stored here, so "the stored entries" is "every cell". That
/// is the one difference from the sparse implementation that is visible in output:
/// `write_Z_to_file` asked for only the stored cells writes the whole container space under this
/// backend. Production only ever asks it for the whole space anyway.
class TStorageZDense : public TDenseStateArray {
public:
	TStorageZDense() = default;
	explicit TStorageZDense(const IndexArray &dimensions) : TDenseStateArray(dimensions) {}

	/// Bulk-insert deferred 0 -> 1 transitions. Mirror of TStorageZMatrix::insert_in_Z.
	void insert_in_Z(const std::vector<std::vector<size_t>> &linear_indices_to_insert) {
		insert_ones_in_batches(*this, linear_indices_to_insert);
	}

	/// The state of every cell of the container space, in ascending linear-index order.
	[[nodiscard]] std::vector<size_t> get_full_Z_binary_vector() const {
		return whole_space_states<size_t>(*this);
	}

	/// Every stored cell as (linear index in Z space, value), in ascending linear-index order --
	/// which here is every cell of the container space.
	[[nodiscard]] std::vector<std::pair<size_t, TStorageZ>> get_stored_entries() const {
		const size_t total = total_size_of_container_space();
		std::vector<std::pair<size_t, TStorageZ>> entries;
		entries.reserve(total);
		for (size_t i = 0; i < total; ++i) { entries.emplace_back(i, TStorageZ(is_one(i))); }
		return entries;
	}
};

static_assert(BinaryFieldStorage<TStorageZDense>,
              "The dense internal state must satisfy the binary storage interface.");
static_assert(!FieldStorage<TStorageZDense>,
              "The internal state carries no posterior counter, so it is not a field.");
