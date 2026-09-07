#pragma once

#include "constants.h"
#include "coretools/algorithms.h"
#include "storages/storage_concepts.h"
#include "storages/y_storage/TStorageY.h"
#include <coretools/Main/TError.h>
#include <cstddef>

class TStorageYDense {
private:
	size_t _thinning_factor = 1;
	size_t _total_counts    = 0;
	std::vector<TStorageY> _vec;
	/// _dimensions_Y_space is the number of leaf nodes in each dimension
	IndexArray _dimensions_Y_space{};

public:
	using iterator       = typename std::vector<TStorageY>::iterator;
	using const_iterator = typename std::vector<TStorageY>::const_iterator;

	void initialize(size_t n_iterations, const IndexArray &dimensions) {
		_thinning_factor = std::max<size_t>(
		    1, static_cast<size_t>(std::ceil(static_cast<double>(n_iterations) /
		                                     static_cast<double>(TStorageY::MAX_COUNTER))));

		_total_counts       = 0;
		_dimensions_Y_space = dimensions;
		_vec.assign(total_size_of_container_space(), TStorageY(false));
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

	[[nodiscard]] bool empty() const { return _vec.empty(); }
	[[nodiscard]] IsOneResult<const_iterator> is_one(size_t linear_index) const {
		DEBUG_ASSERT(linear_index < _vec.size());
		return {_vec[linear_index].is_one(), linear_index < _vec.size(), linear_index,
		        _vec.cbegin() + linear_index};
	}

	[[nodiscard]] IsOneResult<const_iterator> is_one(const IndexArray &index) const {
		const size_t linear_index = get_linear_index_in_container_space(index);
		return {_vec[linear_index].is_one(), linear_index < _vec.size(), linear_index,
		        _vec.cbegin() + linear_index};
	}

	[[nodiscard]] IsOneResult<iterator> is_one(size_t linear_index) {
		DEBUG_ASSERT(linear_index < _vec.size());
		return {_vec[linear_index].is_one(), linear_index < _vec.size(), linear_index,
		        _vec.begin() + linear_index};
	}

	[[nodiscard]] IsOneResult<iterator> is_one(const IndexArray &index) {
		const size_t linear_index = get_linear_index_in_container_space(index);
		return {_vec[linear_index].is_one(), linear_index < _vec.size(), linear_index,
		        _vec.begin() + linear_index};
	}

	void set_state(size_t linear_index, bool state) {
		DEBUG_ASSERT(linear_index < _vec.size());
		_vec[linear_index].set_state(state);
	}

	void insert_one(size_t linear_index) {
		DEBUG_ASSERT(linear_index < _vec.size());
		_vec[linear_index].set_state(true);
	}

	void insert_zero(size_t linear_index) {
		DEBUG_ASSERT(linear_index < _vec.size());
		_vec[linear_index].set_state(false);
	}

	[[nodiscard]] size_t get_linear_index_in_container_space(const IndexArray &index) const {
		DEBUG_ASSERT(index.size() == _dimensions_Y_space.size());
		return coretools::getLinearIndex(index, _dimensions_Y_space);
	}

	[[nodiscard]] IndexArray get_multi_dimensional_index(size_t linear_index_in_Y_space) const {
		DEBUG_ASSERT(linear_index_in_Y_space < _vec.size());
		return coretools::getSubscriptsAsArray(linear_index_in_Y_space, _dimensions_Y_space);
	}

	void remove_zeros() {
		for (auto &i : _vec) {
			if (!i.is_one()) { i.set_counter(0); }
		}
	}

	void reset_counts() {
		for (auto &i : _vec) { i.set_counter(0); }
		_total_counts = 0;
	}

	[[nodiscard]] double get_fraction_of_ones(size_t linear_index) const {
		DEBUG_ASSERT(linear_index < _vec.size());
		if (_total_counts == 0) { return 0.0; }
		return static_cast<double>(_vec[linear_index].get_counter()) /
		       static_cast<double>(_total_counts);
	}

	[[nodiscard]] size_t get_total_counts() const { return _total_counts; }
	[[nodiscard]] size_t get_thinning_factor() const { return _thinning_factor; }
	void add_to_counter(size_t iteration) {
		if (iteration % _thinning_factor != 0) { return; }
		for (auto &i : _vec) {
			if (!i.is_one()) { continue; }
			i.update_counter();
		}
		++_total_counts;
	}

	[[nodiscard]] const IndexArray &dimensions() const { return _dimensions_Y_space; }

	void
	insert_ones_in_container(const std::vector<std::vector<size_t>> &linear_indices_to_insert) {
		for (const auto &indices : linear_indices_to_insert) {
			for (size_t linear_index : indices) {
				DEBUG_ASSERT(linear_index < _vec.size());
				DEBUG_ASSERT(_vec[linear_index].get_counter() == 0);
				_vec[linear_index].set_state(true);
			}
		}
	}
};
static_assert(FieldStorage<TStorageYDense>);
