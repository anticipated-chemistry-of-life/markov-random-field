//
// The one write an updater makes to a cell it may not own.
//

#pragma once

#include "storages/cell_handle.h"
#include "storages/storage_concepts.h"
#include <cstddef>
#include <vector>

/// Writes one cell of `storage`, or records it for a later insert.
///
/// This is the whole write of an update that runs inside a parallel region. The caller names a
/// storage and a linear index, and says nothing about which backend it holds.
///
/// A storage that answers `locate` hands out a handle, and the write is the one branch in
/// cell_handle.h. A sorted-vector matrix has no cell to point at, because it keeps every cell
/// twice, so it makes the write itself and says whether it held the cell. Both answer the same two
/// questions: the write landed, or the cell has to wait.
///
/// A cell the storage holds is written in place. An absent cell that turns into a one waits in
/// `deferred_inserts`, and the caller commits the list once the region ends. An absent cell
/// written to zero is left out, because a cell the storage does not hold already reads as zero.
/// ADR-0006 gives the argument.
///
/// It answers whether the write landed in the storage, which is what a caller that reads the cell
/// back before the commit has to know. Until then the storage still reads the old state, so a
/// walk that reads what it has just written keeps that answer itself
/// (tree/clique/TCliqueView.h). Every other caller ignores it.
template<BinaryStorage Storage>
bool write_or_defer(Storage &storage, size_t linear_index, bool state,
                    std::vector<size_t> &deferred_inserts) {
	if constexpr (LocatableStorage<Storage>) {
		return write_or_defer(storage.locate(linear_index), state, deferred_inserts);
	} else {
		const bool held = storage.write_state_if_held(linear_index, state);
		if (!held && state) { deferred_inserts.push_back(linear_index); }
		return held;
	}
}
