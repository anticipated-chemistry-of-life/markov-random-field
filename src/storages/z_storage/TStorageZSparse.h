//
// The sparse internal state.
//

#pragma once

#include "storages/TSparseState.h"

#include <cstddef>

/// The dense internal state: the dense state array, plus the three bulk paths `Z` is asked for by
/// name. The state array is inherited rather than wrapped because `Z` carries one binary state per
/// cell and nothing else -- there is no member to add, only the whole-space dump, the bulk insert
/// and the stored-entry walk, none of which the storage concept covers (see storage_concepts.h) and
/// all of which the sparse implementation spells with a `Z` in the name too.
///
/// Every cell of the container space is stored here, so "the stored entries" is "every cell". That
/// is the same reading `TDenseStateArray::fill_current_state` gives `exists`, and it is the one
/// difference from the sparse implementation that is visible in output: `write_Z_to_file` asked for
/// only the stored cells writes the whole container space under this backend. Production only ever
/// asks it for the whole space anyway.
class TStorageZSparse : public TSparseBinaryArray {
public:
	TStorageZSparse() = default;
	explicit TStorageZSparse(const IndexArray &dimensions) : TSparseBinaryArray(dimensions) {}
};

static_assert(BinaryFieldStorage<TStorageZSparse>,
              "The dense internal state must satisfy the binary storage interface.");
static_assert(!FieldStorage<TStorageZSparse>,
              "The internal state carries no posterior counter, so it is not a field.");
