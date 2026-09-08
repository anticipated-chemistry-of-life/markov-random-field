//
// The interface the sampler needs from a binary storage, written as C++20 concepts.
//

#pragma once

#include "constants.h"
#include "storages/cell_handle.h"
#include <concepts>
#include <cstddef>
#include <vector>

/// The strided view a storage opens over itself: a start, a count and a stride.
///
/// The sampler reads and writes a run of cells through a window rather than through a cache it
/// owns. Each storage brings the traversal that suits it, and the two are free to differ: the
/// dense window indexes the state vector, and the sparse window materialises its line once. See
/// ADR-0006.
///
/// A window shows its own write to a later read on the same window. The sparse window has to
/// buffer a write to a cell it does not hold, because the insert would reallocate a row inside the
/// parallel region. A read that returned the old state would send the two backends down different
/// chains inside a single update.
///
/// A window ends once, in one of two ways. `close` commits what the window buffered, and the
/// destructor runs it again for a window the caller lets go of. `take_buffered_inserts` hands the
/// buffer to the caller as linear indices and writes nothing, which is the only exit a window
/// inside a parallel region may take. The dense window hands out an empty list, so one loop body
/// serves both backends. ADR-0006 gives the argument.
///
/// A cell is addressed by its position in the window; `linear_index` is what turns that back into
/// the storage's own index.
template<typename T>
concept StorageWindow = requires(T &window, const T &const_window, size_t k, bool state) {
	{ const_window.size() } -> std::same_as<size_t>;
	{ const_window.is_one(k) } -> std::same_as<bool>;
	{ const_window.linear_index(k) } -> std::same_as<size_t>;
	{ window.set_state(k, state) } -> std::same_as<void>;
	{ window.take_buffered_inserts() } -> std::same_as<std::vector<size_t>>;
	{ window.close() } -> std::same_as<void>;
};

/// The surface every binary storage shares: a binary state per cell, the size of the space that
/// state lives in, and the conversion between a linear index and a multi-dimensional one.
///
/// Three kinds of thing satisfy it. The field carries a posterior counter on top (`FieldStorage`
/// below). A node state carries a bare state. The observed data -- the LOTUS records and the
/// simple error model data -- carries a bare state and never changes after it is read in.
///
/// A cell that is not stored reads as state 0, so `is_one` is total over the container space and
/// there is no "does this cell exist" question on this interface. Which cells a storage holds is
/// its own business. `insert_one` and `insert_zero` are the bounds-checked way in for a cell that
/// may be absent; `set_state` is the in-place write for one that is already known to be
/// addressable.
///
/// Three things stay deliberately outside the concept. The bulk-insert and whole-space dump paths
/// still spell `Y` and `Z` in their names. The ones cursor is offered by the field and the
/// observed data alone. The reporting accessors belong to the field. An implementation still has
/// to provide whatever production code calls of those, and the compiler says so, but none of them
/// is what makes a storage a storage.
template<typename T>
concept BinaryStorage = requires(T &storage, const T &const_storage, size_t linear_index,
                                 bool state, const IndexArray &multidim_index) {
	// State.
	{ const_storage.is_one(linear_index) } -> std::same_as<bool>;
	{ storage.set_state(linear_index, state) } -> std::same_as<void>;
	{ storage.insert_one(linear_index) } -> std::same_as<void>;
	{ storage.insert_zero(linear_index) } -> std::same_as<void>;
	{ storage.remove_zeros() } -> std::same_as<void>;

	// Dimensions.
	{ const_storage.total_size_of_container_space() } -> std::same_as<size_t>;
	{ const_storage.empty() } -> std::same_as<bool>;

	// Index conversion.
	{ const_storage.get_linear_index_in_container_space(multidim_index) } -> std::same_as<size_t>;
	{ const_storage.get_multi_dimensional_index(linear_index) } -> std::same_as<IndexArray>;
};

/// A storage that can point an updater at one of its cells.
///
/// `locate` answers the four questions a write asks at once -- the state, whether the storage
/// holds the cell, its linear index and where it is -- so that an in-place write and a deferred
/// insert are one branch apart. `write_or_defer` in cell_handle.h is that branch.
///
/// `is_one` is untouched and still returns a `bool`. The many read-only call sites ask a question
/// with one answer, and they say so.
///
/// The two sorted-vector matrix storages do not satisfy this. A sparse matrix keeps every cell
/// twice, once in its row and once in its column, so it has no single cell to point at. They
/// answer `locate` once they own their cells.
template<typename T>
concept LocatableStorage = BinaryStorage<T> && requires(T &storage, size_t linear_index,
                                                        const IndexArray &multidim_index) {
	typename T::TCell;
	{ storage.locate(linear_index) } -> std::same_as<IsOneResult<typename T::TCell>>;
	{ storage.locate(multidim_index) } -> std::same_as<IsOneResult<typename T::TCell>>;
};

/// A storage the sampler reads and writes through a window: the field and the node states.
///
/// `open_window` is the traversal such a storage brings with it, described at `StorageWindow`
/// above. It is the only way in to a run of cells. The observed data does not satisfy this and
/// does not need to: a data source reads one cell at a time.
template<typename T>
concept WindowedStorage = BinaryStorage<T> && requires(T &storage, const IndexArray &multidim_index,
                                                       size_t n_cells, size_t stride) {
	typename T::TWindow;
	requires StorageWindow<typename T::TWindow>;
	{ storage.open_window(multidim_index, n_cells, stride) } -> std::same_as<typename T::TWindow>;
};

/// The field on top of that: every cell also carries how often it was a one, which is what the
/// posterior fraction of ones is read off at the end of a chain.
///
/// The counter is thinned rather than incremented every iteration, because an implementation may
/// hold it in fewer bits than there are iterations. `add_to_counter` is therefore called with the
/// iteration and decides for itself whether this one counts.
///
/// `get_total_counts` reports how many iterations it has actually counted -- a running total, and
/// deliberately not the arithmetic it looks like. `n_iterations / thinning factor` is wrong twice
/// over: it floors where `add_to_counter` rounds up, which is how a cell that was a one throughout
/// came to report a fraction just above 1; and it assumes the counted iterations start at zero,
/// which they never do, since the caller's index climbs through burn-in and `reset_counts` clears
/// the counters without clearing it. How many multiples of the thinning factor fall in a chain of
/// a given length depends on where that chain started, which is not knowable when the field is
/// sized. Counting removes both questions: the denominator is the numerator's ceiling by
/// construction, so the fraction stays a probability.
template<typename T>
concept FieldStorage = BinaryStorage<T> && requires(T &field, const T &const_field,
                                                    size_t iteration, size_t linear_index) {
	{ field.add_to_counter(iteration) } -> std::same_as<void>;
	{ field.reset_counts() } -> std::same_as<void>;
	{ const_field.get_fraction_of_ones(linear_index) } -> std::same_as<double>;
	{ const_field.get_total_counts() } -> std::same_as<size_t>;
	{ const_field.get_thinning_factor() } -> std::same_as<size_t>;
};
