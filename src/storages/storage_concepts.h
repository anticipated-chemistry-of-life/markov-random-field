//
// The interface the sampler needs from a binary storage, written as C++20 concepts.
//

#pragma once

#include "constants.h"
#include <concepts>
#include <cstddef>
#include <cstdint>
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

/// The surface the field (`Y`) and the internal state (`Z`) share: a binary state per cell,
/// the size of the space that state lives in, and the conversion between a linear index and a
/// multi-dimensional one.
///
/// "Field" here is the Markov-random-field sense -- both `Y` and `Z` are binary fields over their
/// own space -- and not CONTEXT.md's *Field*, which is `Y` alone. That narrower one is what
/// `FieldStorage` below names, and `Z` deliberately does not satisfy it.
///
/// A cell that is not stored reads as state 0, so `is_one` is total over the container space and
/// there is no "does this cell exist" question on the point-lookup path. `insert_one` and
/// `insert_zero` are the bounds-checked way in for a cell that may be absent; `set_state` is the
/// in-place write for one that is already known to be addressable.
///
/// `open_window` is the traversal the storage brings with it: a strided view the sampler reads and
/// writes through, described at `StorageWindow` above.
///
/// `fill_current_state` writes, for the `K` cells starting at `start_index` and running `increment`
/// apart, the state, whether the cell is stored, and the linear index of each. No update calls it
/// any more -- the field update reads its cells through a window (#54) -- so what is left is the
/// tests and the benchmark, and #55 takes it off this interface.
///
/// Deliberately outside the concept: the bulk-insert and whole-space dump paths, which the two
/// implementations still spell with `Y` and `Z` in their names, and the field-only reporting
/// accessors. An implementation still has to provide whatever production code calls of those --
/// the compiler says so -- but they are not part of what makes a storage a storage.
template<typename T>
concept BinaryFieldStorage =
    requires(T &storage, const T &const_storage, size_t linear_index, bool state,
             const IndexArray &multidim_index, size_t n_cells, size_t increment,
             std::vector<uint8_t> &states, std::vector<uint8_t> &exists,
             std::vector<size_t> &linear_indices) {
	    // State.
	    { const_storage.is_one(linear_index) } -> std::same_as<bool>;
	    // Whether the cell is held, as opposed to reading as zero because it is absent. Dense
	    // always says yes; sparse searches. No update asks any more: a window defers the write it
	    // cannot make in place, and answers a later read from its own line. The one caller that
	    // had to know which storage it was talking to has gone. What is left is the tests, and #55
	    // takes it off this interface.
	    { const_storage.is_stored(linear_index) } -> std::same_as<bool>;
	    { storage.set_state(linear_index, state) } -> std::same_as<void>;
	    { storage.insert_one(linear_index) } -> std::same_as<void>;
	    { storage.insert_zero(linear_index) } -> std::same_as<void>;
	    { storage.remove_zeros() } -> std::same_as<void>;

	    // Dimensions.
	    { const_storage.total_size_of_container_space() } -> std::same_as<size_t>;
	    { const_storage.empty() } -> std::same_as<bool>;

	    // Index conversion.
	    {
		    const_storage.get_linear_index_in_container_space(multidim_index)
	    } -> std::same_as<size_t>;
	    { const_storage.get_multi_dimensional_index(linear_index) } -> std::same_as<IndexArray>;

	    // One clique's current state, in one call.
	    {
		    const_storage.fill_current_state(multidim_index, n_cells, increment, states, exists,
		                                     linear_indices)
	    } -> std::same_as<void>;

	    // The strided window the sampler reads and writes through. The field update walks a row of
	    // it, and the node-state walk a column of it.
	    typename T::TWindow;
	    requires StorageWindow<typename T::TWindow>;
	    {
		    storage.open_window(multidim_index, n_cells, increment)
	    } -> std::same_as<typename T::TWindow>;
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
concept FieldStorage = BinaryFieldStorage<T> && requires(T &field, const T &const_field,
                                                         size_t iteration, size_t linear_index) {
	{ field.add_to_counter(iteration) } -> std::same_as<void>;
	{ field.reset_counts() } -> std::same_as<void>;
	{ const_field.get_fraction_of_ones(linear_index) } -> std::same_as<double>;
	{ const_field.get_total_counts() } -> std::same_as<size_t>;
	{ const_field.get_thinning_factor() } -> std::same_as<size_t>;
};
