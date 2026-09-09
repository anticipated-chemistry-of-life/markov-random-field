//
// The handle a storage hands an updater, and the one write that acts on it.
//

#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

/// Everything an updater needs to write one cell: its state, whether the storage holds it, its
/// linear index, and where it is.
///
/// The name is `is_one`'s, because this is what `is_one` answers when the caller is about to
/// write. A read-only caller keeps asking `is_one` and keeps getting a `bool`.
///
/// `cell` is a raw pointer and not an iterator. A hash map iterates over key-value pairs, and a
/// caller that dereferenced one would have to know which half of the pair the state is in. A
/// pointer says "the cell is here" under every storage.
///
/// `in_container` is true exactly when `cell` is not null. A dense storage holds every cell of its
/// container space, so it always answers true. A sparse storage answers true for the cells it was
/// given, and false for the rest -- those are what `write_or_defer` below defers.
///
/// The pointer is valid until the next write that can restructure the storage. An insert or a
/// `remove_zeros` between the `locate` and the write invalidates it.
template<typename Cell> struct IsOneResult {
	bool is_one         = false;
	bool in_container   = false;
	size_t linear_index = 0;
	Cell *cell          = nullptr;
};

/// Whether a cell holds state 1.
///
/// A cell is one of two things. A bare byte carries a state and nothing else, and is its own
/// answer. A cell that carries more than a state -- the packed field cell, whose posterior counter
/// shares its word with the state bit -- answers for itself. Every storage is built over one of
/// the two, so these two overloads are what let one array body serve a state and a counted cell.
inline bool cell_is_one(uint8_t cell) { return cell != 0; }
template<typename Cell> [[nodiscard]] bool cell_is_one(const Cell &cell) { return cell.is_one(); }

/// Writes the state of the cell a handle points at, and leaves everything else the cell carries
/// as it was. A counted cell keeps its counter, which is the rule every state write follows: only
/// an insert starts a counter over.
inline void write_state(uint8_t &cell, bool state) { cell = state ? 1 : 0; }
template<typename Cell> void write_state(Cell &cell, bool state) { cell.set_state(state); }

/// Writes a cell the storage holds, or records the cell for a later insert.
///
/// This is the whole difference between the two backends, and it is one branch. A dense storage
/// holds every cell, so the write always lands here. A sparse storage cannot insert a cell inside
/// a parallel region -- the insert restructures the container while other threads read it -- so an
/// absent cell that turns into a one waits in `deferred_inserts`, and the caller commits the list
/// once the region ends.
///
/// An absent cell written to zero is left out. A cell the storage does not hold already reads as
/// zero, so inserting it would store a cell to say what its absence says.
///
/// The helper keeps no record of the cell. A caller that writes one absent cell twice before it
/// commits defers it twice, and a one followed by a zero still inserts a one. An update visits
/// each cell once, which is what makes that a non-question.
///
/// The helper takes a handle and not a storage. One body therefore serves every storage, and it is
/// a template only over the cell the handle points at.
///
/// It answers whether the write landed in the storage. A caller that reads the cell back before
/// the deferred list is committed needs that answer, because the storage still reads the old
/// state; every other caller ignores it.
template<typename Cell>
bool write_or_defer(const IsOneResult<Cell> &handle, bool state,
                    std::vector<size_t> &deferred_inserts) {
	if (handle.in_container) {
		write_state(*handle.cell, state);
	} else if (state) {
		deferred_inserts.push_back(handle.linear_index);
	}
	return handle.in_container;
}
