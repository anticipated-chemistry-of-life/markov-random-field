//
// Created by VISANI Marco on 18.10.2024.
//

#ifndef TSTORAGEZVECTOR_H
#define TSTORAGEZVECTOR_H
#include "TStorageZ.h"
#include <cassert>
#include <vector>

class TStorageZVector {
private:
	std::vector<TStorageZ> _vec;
	[[nodiscard]] std::pair<bool, size_t> _binary_search(uint32_t coordinate) const {

		// lower_bound return the first element that is not less than the value
		auto it = std::lower_bound(_vec.begin(), _vec.end(), coordinate);

		// if our coordinate is bigger than the biggest element in the vector
		// we say that we haven't found our element and that if we want to
		// insert it, we should insert it at the end of the vector
		if (it == _vec.end()) { return {false, _vec.size()}; }

		// else our coordinate is in the range of the coordinates in the vector
		// meaning that if we haven't found it, we will insert it at that position
		// to keep the vector sorted
		size_t index = std::distance(_vec.begin(), it);
		if (*it != coordinate) { return {false, index}; }

		// if we found the coordinate we return the index and true
		return {true, index};
	};

public:
	TStorageZVector() = default;
	[[nodiscard]] bool is_one(const int32_t coordinate) const {
		const int32_t abs_coordinate = std::abs(coordinate);
		auto [found, index]          = _binary_search(abs_coordinate);
		if (found) { return _vec[index].is_one(); }
		return false;
	};

	void set_to_one(int32_t coordinate) {
		assert(0 < coordinate && "Coordinate must be positive");
		const int32_t abs_coordinate = std::abs(coordinate);
		auto [found, index]          = _binary_search(abs_coordinate);
		if (found) {
			_vec[index].set_state(true);
		} else {
			_vec.insert(_vec.begin() + index, TStorageZ(coordinate));
		}
	}

	void set_to_zero(int32_t coordinate) {
		const int32_t abs_coordinate = std::abs(coordinate);
		auto [found, index]          = _binary_search(abs_coordinate);
		if (found) { _vec[index].set_state(false); }
	}

	void remove_zeros() {
		_vec.erase(std::remove_if(_vec.begin(), _vec.end(), [](const TStorageZ &storage) { return !storage.is_one(); }),
		           _vec.end());
	}
	size_t size() const { return _vec.size(); }
};

#endif // TSTORAGEZVECTOR_H
