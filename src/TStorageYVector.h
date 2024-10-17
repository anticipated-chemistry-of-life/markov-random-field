//
// Created by VISANI Marco on 17.10.2024.
//

#ifndef TSTORAGEYVECTOR_H
#define TSTORAGEYVECTOR_H
#include "TStorageY.h"

#include <vector>

class TStorageYVector {
private:
	size_t _thinning_factor;
	size_t _total_counts;
	std::vector<TStorageY> _vec;
	std::pair<bool, size_t> _binary_search(uint64_t coordinate) const {

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
	TStorageYVector(size_t n_iterations) {
		constexpr uint16_t max_value = std::numeric_limits<uint16_t>::max();
		_thinning_factor             = std::ceil(static_cast<double>(n_iterations) / static_cast<double>(max_value));
		_total_counts                = n_iterations / _thinning_factor;
	};

	// we want to check if coordinate exists in the vector or not.
	// the coordinate is the position of the element in the Y dimension
	// so if an element is in the vector that means it is currently a one or
	// it has been in the past iteration.
	bool is_one(uint64_t coordinate) const {
		auto [found, index] = _binary_search(coordinate);
		if (found) { return _vec[index].is_one(); }
		return false;
	};

	/* set_to_one will set the element at the coordinate to 1
	 * if the element is already in the vector, we just set it to 1
	 * if the element is not in the vector, we insert it in the vector
	 * and set it to 1
	 */
	void set_to_one(uint64_t coordinate) {
		auto [found, index] = _binary_search(coordinate);
		if (found) {
			_vec[index].set_state(true);
		} else {
			_vec.insert(_vec.begin() + index, TStorageY(coordinate));
		}
	}

	void set_to_zero(uint64_t coordinate) {
		auto [found, index] = _binary_search(coordinate);
		if (found) { _vec[index].set_state(false); }
	}

	void add_to_counter(size_t iteration) {
		if (iteration % _thinning_factor == 0) {
			for (auto &elem : _vec) { elem.update_counter(); }
		}
	}

	double get_fraction_of_ones(uint64_t coordinate) const {
		auto [found, index] = _binary_search(coordinate);
		if (!found) { return 0.0; }
		return static_cast<double>(_vec[index].get_counter()) / static_cast<double>(_total_counts);
	}

	size_t get_total_counts() const { return _total_counts; }

	size_t length_of_vector() const { return _vec.size(); }
};

#endif // TSTORAGEYVECTOR_H
