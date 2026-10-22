//
// The dense internal state.
//

#pragma once

#include "storages/TDense.h"
#include "storages/bulk_paths.h"

#include <cstddef>
#include <vector>

/// The dense internal state: the dense state array, plus the two bulk paths `Z` is asked for by
/// name. The state array is inherited rather than wrapped because `Z` carries one binary state per
/// cell and nothing else -- there is no member to add, only the whole-space dump and the bulk
/// insert, neither of which the storage concept covers (see storage_concepts.h) and both of which
/// the sparse implementation spells with a `Z` in the name too.
///
/// `locate` comes with the array too, so the node state points an updater at one of its cells
/// without a member of its own.
///
/// The array is laid out by whoever sizes it. A tree sizes its node state with its own dimension
/// contiguous, so each of its cliques is one run of adjacent cells (ADR-0013).
class TStorageZDense : public TDenseBinary {
public:
	TStorageZDense() = default;
	explicit TStorageZDense(const IndexArray &dimensions,
	                        size_t contiguous_dimension = NUMBER_OF_TREES - 1)
	    : TDenseBinary(dimensions, contiguous_dimension) {}

	/// Bulk-insert deferred 0 -> 1 transitions. Mirror of TStorageZSparse::insert_in_Z.
	void insert_in_Z(const std::vector<std::vector<size_t>> &linear_indices_to_insert) {
		insert_ones_in_batches(*this, linear_indices_to_insert);
	}

	/// The state of every cell of the container space, in ascending linear-index order.
	[[nodiscard]] std::vector<size_t> get_full_Z_binary_vector() const {
		return whole_space_states<size_t>(*this);
	}
};

static_assert(BinaryStorage<TStorageZDense>,
              "The dense node state must satisfy the binary storage interface.");
static_assert(!FieldStorage<TStorageZDense>,
              "The node state carries no posterior counter, so it is not a field.");
