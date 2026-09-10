//
// The two bulk paths every storage has in common, under the names the sampler asks for them by.
//

#pragma once

#include "storages/storage_concepts.h"
#include <cstddef>
#include <vector>

// Each storage has to answer to a name of its own -- `insert_in_Y` beside `insert_in_Z`,
// `get_full_Y_binary_vector` beside `get_full_Z_binary_vector`, because that is what the sampler
// asks for -- so the shape the four share lives here rather than in a member they could inherit.
// Neither is part of the storage concept (see storage_concepts.h): they describe what the model
// does with a storage rather than what makes a storage a storage.

/// Writes a one at every index of every batch, which is how an update's deferred inserts land once
/// its parallel region has ended.
///
/// It goes through the storage rather than through its cells because an insert does something a
/// state write does not: it writes the cell anew, so a counted cell starts its counter over. Both
/// backends have to leave a cell holding the same thing, and an update only ever defers a cell it
/// found absent, whose counter is 0 either way.
template<BinaryStorage Storage>
void insert_ones_in_batches(Storage &storage, const std::vector<std::vector<size_t>> &batches) {
	for (const auto &batch : batches) {
		for (const size_t linear_index : batch) { storage.insert_one(linear_index); }
	}
}

/// The state of every cell of the container space, in ascending linear-index order.
///
/// The element type is the caller's because the two fields and the two node states disagree on it
/// -- a field dumps a vector of bytes and a node state a vector of words -- and a trace line is
/// written from whatever they return.
template<typename T, BinaryStorage Storage>
[[nodiscard]] std::vector<T> whole_space_states(const Storage &storage) {
	const size_t total = storage.total_size_of_container_space();
	std::vector<T> states;
	states.reserve(total);
	for (size_t i = 0; i < total; ++i) { states.push_back(static_cast<T>(storage.is_one(i))); }
	return states;
}
