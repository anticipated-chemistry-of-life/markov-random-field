//
// Created by VISANI Marco on 17.10.2024.
//

#ifndef TSTORAGEYVECTOR_H
#define TSTORAGEYVECTOR_H
#include "TStorageY.h"
#include "Types.h"
#include "coretools/algorithms.h"

#include <cstddef>
#include <cstdint>
#include <vector>

class TStorageYVector {
private:
	size_t _thinning_factor;
	size_t _total_counts;
	std::vector<TStorageY> _vec;

	/// _dimensions_Y_space is the number of leaf nodes in each dimension
	std::vector<size_t> _dimensions_Y_space;

public:
	using value_type     = uint64_t;
	using const_iterator = typename std::vector<TStorageY>::const_iterator;
	TStorageYVector()    = default;
	~TStorageYVector()   = default;
	TStorageYVector(const size_t n_iterations, const std::vector<size_t> &dimensions_Y_space) {
		initialize(n_iterations, dimensions_Y_space);
	};

	void initialize(const size_t n_iterations, const std::vector<size_t> &dimensions_Y_space) {
		constexpr uint16_t max_value = std::numeric_limits<uint16_t>::max();
		_thinning_factor             = std::ceil(static_cast<double>(n_iterations) / static_cast<double>(max_value));
		_total_counts                = n_iterations / _thinning_factor;
		_dimensions_Y_space          = dimensions_Y_space;
	}

	/// We want to check if element at position index_in_TStorageYVector is one.
	/// This index *must* be previously obtained by a binary search.
	/// @param index_in_TStorageYVector the position of the element in the Y vector.
	/// @return true if the element is one, false otherwise.
	[[nodiscard]] bool is_one(const size_t index_in_TStorageYVector) const {
		if (index_in_TStorageYVector >= _vec.size()) {
			UERROR("Index '", index_in_TStorageYVector,
			       "' is out of range. The length of the vector is : ", _vec.size(), ".");
		}
		return _vec[index_in_TStorageYVector].is_one();
	};

	/** set_to_one will set the element at the index_in_TStorageYVector to 1.
	 * @param index_in_TStorageYVector the position of the element in the Y vector
	 */
	void set_to_one(size_t index_in_TStorageYVector) {
		if (index_in_TStorageYVector >= _vec.size()) {
			UERROR("Index '", index_in_TStorageYVector,
			       "' is out of range. The length of the vector is : ", _vec.size());
		}
		_vec[index_in_TStorageYVector].set_state(true);
	}
	void set_to_zero(size_t index_in_TStorageYVector) {
		if (index_in_TStorageYVector >= _vec.size()) {
			UERROR("Index '", index_in_TStorageYVector,
			       "' is out of range. The length of the vector is : ", _vec.size());
		}
		_vec[index_in_TStorageYVector].set_state(false);
	}

	void insert_one(uint64_t linear_index_in_Y_space) {
		auto [found, index] = binary_search(linear_index_in_Y_space);
		if (found) {
			_vec[index].set_state(true);
		} else {
			_vec.insert(_vec.begin() + index, TStorageY(linear_index_in_Y_space));
		}
	}

	/// Does the same as set_to_one but sets the element to zero
	void insert_zero(uint64_t linear_index_in_Y_space) {
		auto [found, index] = binary_search(linear_index_in_Y_space);
		if (found) {
			_vec[index].set_state(false);
		} else {
			_vec.insert(_vec.begin() + index, TStorageY(linear_index_in_Y_space));
			_vec[index].set_state(false);
		}
	}

	void add_to_counter(size_t iteration) {
		if (iteration % _thinning_factor == 0) {
// TODO: can we do this in parallel ?
#pragma omp parallel for num_threads(NUMBER_OF_THREADS)
			for (auto &elem : _vec) { elem.update_counter(); }
		}
	}

	double get_fraction_of_ones(uint64_t linear_index_in_Y_space) const {
		auto [found, index] = binary_search(linear_index_in_Y_space);
		if (!found) { return 0.0; }
		return static_cast<double>(_vec[index].get_counter()) / static_cast<double>(_total_counts);
	}

	size_t get_total_counts() const { return _total_counts; }

	size_t size() const { return _vec.size(); }
	auto begin() const { return _vec.begin(); }
	auto end() const { return _vec.end(); }
	const std::vector<TStorageY> &get_vector() const { return _vec; }
	const TStorageY &operator[](size_t index_in_TStorageYVector) const { return _vec[index_in_TStorageYVector]; }

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
	uint64_t get_linear_index_in_Y_space(const std::vector<size_t> &multidim_index_in_Y_space) const {
		return coretools::getLinearIndex(multidim_index_in_Y_space, _dimensions_Y_space);
	}
	uint64_t get_linear_index_in_container_space(const std::vector<size_t> &multidim_index_in_Y_space) const {
		return get_linear_index_in_Y_space(multidim_index_in_Y_space);
	}

	/// Given a linear index, we want to get the multi-dimensional index.
	[[nodiscard]] std::vector<size_t> get_multi_dimensional_index(uint64_t linear_index_in_Y_space) const {
		auto tmp = static_cast<size_t>(linear_index_in_Y_space);
		return coretools::getSubscripts(tmp, _dimensions_Y_space);
	}

	/// Returns the product of the dimensions in the container. This is the
	/// maximum number of ones that can be stored in the vector given
	/// the dimensions of the container. For example, if the container
	/// has dimension sizes [2, 3, 4], the maximum number of ones that can be
	/// stored in the vector is 24.
	[[nodiscard]]
	size_t total_size_of_container_space() const {
		return coretools::containerProduct(_dimensions_Y_space);
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

	[[nodiscard]] std::pair<bool, size_t> binary_search(const std::vector<size_t> &multidim_index_in_Y_space) const {
		return binary_search(get_linear_index_in_Y_space(multidim_index_in_Y_space));
	}

	void insert_in_Y(const std::vector<std::vector<TStorageY>> &linear_indices_in_Y_space_to_insert) {
		auto size_to_insert =
		    std::accumulate(linear_indices_in_Y_space_to_insert.begin(), linear_indices_in_Y_space_to_insert.end(), 0,
		                    [](size_t sum, const std::vector<TStorageY> &i) { return sum + i.size(); });

		this->_vec.reserve(this->_vec.size() + size_to_insert);

		for (const auto &vec : linear_indices_in_Y_space_to_insert) {
			this->_vec.insert(_vec.end(), vec.begin(), vec.end());
		}
		std::sort(_vec.begin(), _vec.end());
	}

	std::vector<size_t> get_full_Y_binary_vector() const {
		std::vector<size_t> Y_as_vector;
		for (size_t i = 0; i < total_size_of_container_space(); ++i) {
			auto [found, idx] = binary_search(i);
			if (found) {
				Y_as_vector.push_back(_vec[idx].is_one());
			} else {
				Y_as_vector.push_back(false);
			}
		}
		return Y_as_vector;
	}

	bool empty() const { return _vec.empty(); }
};

#endif // TSTORAGEYVECTOR_H
