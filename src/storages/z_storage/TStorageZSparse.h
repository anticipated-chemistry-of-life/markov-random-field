//
// The sparse node state.
//

#pragma once

#include "TStorageZ.h"
#include "constants.h"
#include "storages/TSparse.h"
#include "storages/bulk_paths.h"
#include "storages/storage_concepts.h"

#include <cstddef>
#include <utility>
#include <vector>

/// The sparse node state: the sparse array of bare states, plus the three bulk paths `Z` is asked
/// for by name. The array is inherited rather than wrapped because `Z` carries one binary state
/// per cell and nothing else -- there is no member to add, only the whole-space dump, the bulk
/// insert and the stored-entry walk, none of which the storage concept covers (see
/// storage_concepts.h) and all of which the dense implementation spells with a `Z` in the name too.
///
/// `locate` comes with the array too, so the node state points an updater at one of its cells
/// without a member of its own.
///
/// There is one node state per dimension `d`. The dimensions in Z space are the number of leaves
/// in each dimension except for `d`, which spans every node of its tree -- leaves included
/// (ADR-0005). For example, with leaf counts [2, 3] and 19 nodes in dimension `d == 0`, the
/// dimensions in Z space are [19, 3], where the first two rows are the leaves the field also
/// holds. `node_state_dimensions` in tree/node_state_shape.h is where that rule lives.
///
/// Only the cells a caller put in are stored, so "the stored entries" is a subset of the container
/// space where the dense implementation reports all of it. That is the one difference between the
/// two that is visible in output: `write_Z_to_file` asked for only the stored cells writes fewer
/// rows under this backend. Production only ever asks it for the whole space anyway.
class TStorageZSparse : public TSparseBinary {
public:
	TStorageZSparse() = default;
	explicit TStorageZSparse(const IndexArray &dimensions) : TSparseBinary(dimensions) {}

	/// Bulk-insert deferred 0 -> 1 transitions. Mirror of TStorageZDense::insert_in_Z.
	void insert_in_Z(const std::vector<std::vector<size_t>> &linear_indices_to_insert) {
		insert_ones_in_batches(*this, linear_indices_to_insert);
	}

	/// The state of every cell of the container space, in ascending linear-index order.
	[[nodiscard]] std::vector<size_t> get_full_Z_binary_vector() const {
		return whole_space_states<size_t>(*this);
	}

	/// Every stored cell as (linear index in Z space, value), in ascending linear-index order.
	[[nodiscard]] std::vector<std::pair<size_t, TStorageZ>> get_stored_entries() const {
		std::vector<std::pair<size_t, TStorageZ>> entries;
		entries.reserve(size());
		for (const auto &[linear_index, state] : stored_cells_in_order()) {
			entries.emplace_back(linear_index, TStorageZ(state != 0));
		}
		return entries;
	}
};

static_assert(BinaryStorage<TStorageZSparse>,
              "The sparse node state must satisfy the binary storage interface.");
static_assert(!FieldStorage<TStorageZSparse>,
              "The node state carries no posterior counter, so it is not a field.");
