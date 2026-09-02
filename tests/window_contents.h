//
// A window's cells as vectors, for the tests that compare a whole line in one assertion.
//

#pragma once

#include "storages/storage_concepts.h"
#include <cstddef>
#include <cstdint>
#include <vector>

/// The state of every cell of a window, as a vector a test can compare in one assertion.
///
/// A window answers one cell at a time. A loop of EXPECT_EQ over it prints one failure line per
/// cell. Comparing two vectors prints the whole line instead. That is what says *where* a backend
/// went wrong, and not only that it did.
template<StorageWindow Window> std::vector<uint8_t> states_of(const Window &window) {
	std::vector<uint8_t> states;
	states.reserve(window.size());
	for (size_t k = 0; k < window.size(); ++k) {
		states.push_back(static_cast<uint8_t>(window.is_one(k)));
	}
	return states;
}

/// The linear index of every cell of a window, in the same shape.
template<StorageWindow Window> std::vector<size_t> linear_indices_of(const Window &window) {
	std::vector<size_t> linear_indices;
	linear_indices.reserve(window.size());
	for (size_t k = 0; k < window.size(); ++k) { linear_indices.push_back(window.linear_index(k)); }
	return linear_indices;
}
