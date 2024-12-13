//
// Created by VISANI Marco on 18.10.2024.
//

#ifndef TSTORAGEZVECTOR_H
#define TSTORAGEZVECTOR_H
#include "TStorageZ.h"
#include "coretools/Main/TError.h"
#include "coretools/algorithms.h"
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <numeric>
#include <vector>

/// There is one TStorageZVector per dimension `d`. The dimensions in Z space correspond
/// to the number of leaves in each dimension except for dimension `d` where the
/// dimension is given by the number of internal nodes.
/// For example dimension of interest `d` where d = 2 and number of leaves
/// in each dimension is [2, 3, 4, 5] and the number of internal nodes in dimension `d` is 17.
/// then the dimensions in Z space are [2, 3, 17, 5].
class TStorageZVector {
private:
	std::vector<size_t> _dimensions_in_Z_space;
	std::vector<TStorageZ> _vec;

public:
	using value_type     = uint32_t;
	using const_iterator = typename std::vector<TStorageZ>::const_iterator;
	TStorageZVector()    = default;
	explicit TStorageZVector(const std::vector<size_t> &dimensions_in_Z_space) {
		initialize_dimensions(dimensions_in_Z_space);
	}

	void initialize_dimensions(const std::vector<size_t> &dimensions_in_Z_space) {
		_dimensions_in_Z_space = dimensions_in_Z_space;
	}
	[[nodiscard]] bool is_one(const size_t index_in_TStorageZVector) const {
		if (index_in_TStorageZVector >= _vec.size()) {
			UERROR("Index '", index_in_TStorageZVector,
			       "' is out of range. The length of the vector is : ", _vec.size());
		}
		return _vec[index_in_TStorageZVector].is_one();
	};

	void set_to_one(size_t index_in_TStorageZVector) {
		if (index_in_TStorageZVector >= _vec.size()) {
			UERROR("Index '", index_in_TStorageZVector,
			       "' is out of range. The length of the vector is : ", _vec.size());
		}
		_vec[index_in_TStorageZVector].set_state(true);
	}

	void set_to_zero(size_t index_in_TStorageZVector) {
		if (index_in_TStorageZVector >= _vec.size()) {
			UERROR("Index '", index_in_TStorageZVector,
			       "' is out of range. The length of the vector is : ", _vec.size());
		}
		_vec[index_in_TStorageZVector].set_state(false);
	}

	void insert_one(uint32_t linear_index_in_Z_space) {
		auto [found, index] = binary_search(linear_index_in_Z_space);
		if (found) {
			_vec[index].set_state(true);
		} else {
			_vec.insert(_vec.begin() + index, TStorageZ(linear_index_in_Z_space));
		}
	}

	void insert_zero(uint32_t linear_index_in_Z_space) {
		auto [found, index] = binary_search(linear_index_in_Z_space);
		if (found) {
			_vec[index].set_state(false);
		} else {
			_vec.insert(_vec.begin() + index, TStorageZ(linear_index_in_Z_space));
			_vec[index].set_state(false);
		}
	}

	void remove_zeros() {
		_vec.erase(std::remove_if(_vec.begin(), _vec.end(), [](const TStorageZ &storage) { return !storage.is_one(); }),
		           _vec.end());
	}
	size_t size() const { return _vec.size(); }
	auto begin() const { return _vec.begin(); }
	auto cbegin() const { return _vec.cbegin(); }
	auto cend() const { return _vec.cend(); }
	auto end() const { return _vec.end(); }
	TStorageZ operator[](size_t index_in_TStorageZVector) const { return _vec[index_in_TStorageZVector]; }

	[[nodiscard]] uint64_t get_linear_index_in_Z_space(const std::vector<size_t> &multidim_index_in_Z_space) const {
		return coretools::getLinearIndex(multidim_index_in_Z_space, _dimensions_in_Z_space);
	}

	uint64_t get_linear_index_in_container_space(const std::vector<size_t> &multidim_index_in_Z_space) const {
		return get_linear_index_in_Z_space(multidim_index_in_Z_space);
	}

	[[nodiscard]] std::pair<bool, size_t> binary_search(uint32_t linear_index_in_Z_space) const {

		// lower_bound return the first element that is not less than the value
		auto it = std::lower_bound(_vec.begin(), _vec.end(), linear_index_in_Z_space);

		// if our coordinate is bigger than the biggest element in the vector
		// we say that we haven't found our element and that if we want to
		// insert it, we should insert it at the end of the vector
		if (it == _vec.end()) { return {false, _vec.size()}; }

		// else our coordinate is in the range of the coordinates in the vector
		// meaning that if we haven't found it, we will insert it at that position
		// to keep the vector sorted
		size_t index = std::distance(_vec.begin(), it);
		if (*it != linear_index_in_Z_space) { return {false, index}; }

		// if we found the coordinate we return the index and true
		return {true, index};
	};

	[[nodiscard]] std::pair<bool, size_t> binary_search(const std::vector<size_t> &multidim_index_in_Z_space) const {
		return binary_search(get_linear_index_in_Z_space(multidim_index_in_Z_space));
	}

	void insert_in_Z(const std::vector<std::vector<TStorageZ>> &linear_indices_in_Z_space_to_insert) {
		auto size_to_insert =
		    std::accumulate(linear_indices_in_Z_space_to_insert.begin(), linear_indices_in_Z_space_to_insert.end(), 0,
		                    [](size_t sum, const std::vector<TStorageZ> &i) { return sum + i.size(); });

		this->_vec.reserve(this->_vec.size() + size_to_insert);

		for (const auto &vec : linear_indices_in_Z_space_to_insert) {
			this->_vec.insert(_vec.end(), vec.begin(), vec.end());
		}
		std::sort(_vec.begin(), _vec.end());
	}
};

#endif // TSTORAGEZVECTOR_H
