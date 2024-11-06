//
// Created by VISANI Marco on 17.10.2024.
//

#ifndef TSTORAGEYVECTOR_H
#define TSTORAGEYVECTOR_H
#include "TStorageY.h"
#include "coretools/algorithms.h"

#include <cstddef>
#include <cstdint>
#include <vector>

class TStorageYVector {
private:
	size_t _thinning_factor;
	size_t _total_counts;
	std::vector<TStorageY> _vec;

	/// Gives the number of dimensions of our model. This should usually be 2 (species and molecules) but as we have
	/// defined in our model, this can be more.
	std::vector<size_t> _dimensions;

public:
	using value_type     = uint64_t;
	using const_iterator = typename std::vector<TStorageY>::const_iterator;
	explicit TStorageYVector(const size_t n_iterations, const std::vector<size_t> &dimensions) {
		// TODO : NOTE that dimensions correspond to the number of leaf nodes in each dimension !!!
		constexpr uint16_t max_value = std::numeric_limits<uint16_t>::max();
		_thinning_factor             = std::ceil(static_cast<double>(n_iterations) / static_cast<double>(max_value));
		_total_counts                = n_iterations / _thinning_factor;
		_dimensions                  = dimensions;
	};

	/// we want to check if coordinate exists in the vector or not.
	/// the coordinate is the position of the element in the Y dimension
	/// so if an element is in the vector that means it is currently a one or
	/// it has been in the past iteration.
	/// @param coordinate the position of the element in the Y vector.
	/// @return true if the element is one, false otherwise.
	[[nodiscard]] bool is_one(const uint64_t index_in_Y) const { return _vec[index_in_Y].is_one(); };

	/** set_to_one will set the element at the coordinate to 1
	 * if the element is already in the vector, we just set it to 1
	 * if the element is not in the vector, we insert it in the vector
	 * and set it to 1
	 * @param coordinate the position of the element in the Y vector
	 */
	void set_to_one(uint64_t coordinate) {
		auto [found, index] = binary_search(coordinate);
		if (found) {
			_vec[index].set_state(true);
		} else {
			_vec.insert(_vec.begin() + index, TStorageY(coordinate));
		}
	}

	/// Does the same as set_to_one but sets the element to zero
	/// @param coordinate the position of the element in the Y vector
	void set_to_zero(uint64_t coordinate) {
		auto [found, index] = binary_search(coordinate);
		if (found) { _vec[index].set_state(false); }
	}

	void add_to_counter(size_t iteration) {
		if (iteration % _thinning_factor == 0) {
			for (auto &elem : _vec) { elem.update_counter(); }
		}
	}

	double get_fraction_of_ones(uint64_t coordinate) const {
		auto [found, index] = binary_search(coordinate);
		if (!found) { return 0.0; }
		return static_cast<double>(_vec[index].get_counter()) / static_cast<double>(_total_counts);
	}

	size_t get_total_counts() const { return _total_counts; }

	size_t size() const { return _vec.size(); }
	auto begin() const { return _vec.begin(); }
	auto end() const { return _vec.end(); }
	const std::vector<TStorageY> &get_vector() const { return _vec; }
	const TStorageY &operator[](size_t index) const { return _vec[index]; }

	void reset_counts() {
		for (auto &elem : _vec) { elem.reset_counter(); }
	}

	/// Remove all the elements that have the state to zero.
	/// @return void
	void remove_zeros() {
		_vec.erase(std::remove_if(_vec.begin(), _vec.end(), [](const TStorageY &elem) { return !elem.is_one(); }),
		           _vec.end());
	}

	/// Given a multi-dimensional index, we want to get its linear index.
	/// @param multi_dim_index the multi-dimensional index
	/// @return the linear index
	uint64_t get_linear_coordinate(const std::vector<size_t> &multi_dim_index) const {
		return coretools::getLinearIndex(multi_dim_index, _dimensions);
	}

	/// @brief Binary search to find the coordinate in the vector
	[[nodiscard]] std::pair<bool, size_t> binary_search(uint64_t coordinate) const {

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

	[[nodiscard]] std::pair<bool, size_t> binary_search(const std::vector<size_t> &multi_dim_index) const {
		return binary_search(get_coordinate(multi_dim_index));
	}
};

#endif // TSTORAGEYVECTOR_H
