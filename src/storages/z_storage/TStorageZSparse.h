//
// The sparse node state.
//

#pragma once

#include "constants.h"
#include "storages/TSparse.h"
#include "storages/bulk_paths.h"
#include "storages/storage_concepts.h"

#include <cstddef>
#include <vector>

/// The sparse node state: the sparse array of bare states, plus the two bulk paths `Z` is asked
/// for by name. The array is inherited rather than wrapped because `Z` carries one binary state
/// per cell and nothing else -- there is no member to add, only the whole-space dump and the bulk
/// insert, neither of which the storage concept covers (see storage_concepts.h) and both of which
/// the dense implementation spells with a `Z` in the name too.
///
/// `locate` comes with the array too, so the node state points an updater at one of its cells
/// without a member of its own.
///
/// There is one node state per dimension `d`. The dimensions in Z space are the number of leaves
/// in each dimension except for `d`, which spans every node of its tree -- leaves included
/// (ADR-0005). For example, with leaf counts [2, 3] and 19 nodes in dimension `d == 0`, the
/// dimensions in Z space are [19, 3], where the first two rows are the leaves the field also
/// holds. `node_state_dimensions` in tree/node_state_shape.h is where that rule lives.
class TStorageZSparse : public TSparseBinary {
public:
	TStorageZSparse() = default;
	explicit TStorageZSparse(const IndexArray &dimensions,
	                         size_t contiguous_dimension = NUMBER_OF_TREES - 1)
	    : TSparseBinary(dimensions, contiguous_dimension) {}

	/// Bulk-insert deferred 0 -> 1 transitions. Mirror of TStorageZDense::insert_in_Z.
	void insert_in_Z(const std::vector<std::vector<size_t>> &linear_indices_to_insert) {
		insert_ones_in_batches(*this, linear_indices_to_insert);
	}

	/// The state of every cell of the container space, in ascending linear-index order.
	[[nodiscard]] std::vector<size_t> get_full_Z_binary_vector() const {
		return whole_space_states<size_t>(*this);
	}
};

static_assert(BinaryStorage<TStorageZSparse>,
              "The sparse node state must satisfy the binary storage interface.");
static_assert(!FieldStorage<TStorageZSparse>,
              "The node state carries no posterior counter, so it is not a field.");
