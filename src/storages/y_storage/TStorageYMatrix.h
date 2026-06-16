//
// Created by VISANI Marco on 17.10.2024.
//

#pragma once

#include "TStorageY.h"
#include "cli.h"
#include "constants.h"
#include "coretools/Main/TError.h"
#include "coretools/Math/TSparseMatrix.h"
#include "coretools/algorithms.h"
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

class TStorageYMatrix {
private:
	size_t _thinning_factor;
	size_t _total_counts;
	coretools::TSparseMatrix<TStorageY> _mat;

	/// _dimensions_Y_space is the number of leaf nodes in each dimension
	std::array<size_t, NUMBER_OF_TREES> _dimensions_Y_space;

	void _insert(uint64_t linear_index_in_Y_space, bool state) {
		if (linear_index_in_Y_space >= this->total_size_of_container_space()) {
			throw coretools::TDevError(
			    "You are trying to insert a value at an linear index bigger than the total "
			    "size of the container. ",
			    "The index is: ", linear_index_in_Y_space,
			    " and the total size of the container is : ",
			    this->total_size_of_container_space());
		}
		const auto multidim_index = this->get_multi_dimensional_index(linear_index_in_Y_space);
		_mat.set(multidim_index[0], multidim_index[1], TStorageY(state));
	};

public:
	TStorageYMatrix()  = default;
	~TStorageYMatrix() = default;
	TStorageYMatrix(const size_t n_iterations,
	                const std::array<size_t, NUMBER_OF_TREES> &dimensions_Y_space) {
		initialize(n_iterations, dimensions_Y_space);
	};

	void initialize(const size_t n_iterations,
	                const std::array<size_t, NUMBER_OF_TREES> &dimensions_Y_space) {
		constexpr int16_t max_value = std::numeric_limits<int16_t>::max();
		_thinning_factor =
		    std::ceil(static_cast<double>(n_iterations) / static_cast<double>(max_value));
		_total_counts       = n_iterations / _thinning_factor;
		_dimensions_Y_space = dimensions_Y_space;
		_mat.resize(dimensions_Y_space[0], dimensions_Y_space[1]);
	}

	void initialize(const size_t n_iterations, const std::vector<size_t> &dimensions_Y_space) {
		if (dimensions_Y_space.size() != NUMBER_OF_TREES) {
			throw coretools::TDevError("dimensions_Y_space must have size NUMBER_OF_TREES");
		}
		std::array<size_t, NUMBER_OF_TREES> dims = {dimensions_Y_space[0], dimensions_Y_space[1]};
		initialize(n_iterations, dims);
	}

	/// We want to check if element at position index_in_TStorageYMatrix is one.
	/// This index *must* be previously obtained by a binary search.
	/// @param index_in_TStorageYMatrix the position of the element in the Y vector.
	/// @return true if the element is one, false otherwise.
	[[nodiscard]] inline bool is_one(const size_t index_in_TStorageYMatrix) const {
		const auto multidim_index = get_multi_dimensional_index(index_in_TStorageYMatrix);
		return _mat.get(multidim_index[0], multidim_index[1]).is_one();
	};

	/** set_to_one will set the element at the index_in_TStorageYMatrix to 1.
	 * @param index_in_TStorageYMatrix the position of the element in the Y vector
	 */
	void set_to_one(size_t index_in_TStorageYMatrix) {
		const auto multidim_index = get_multi_dimensional_index(index_in_TStorageYMatrix);
		auto s                    = _mat.get(multidim_index[0], multidim_index[1]);
		s.set_state(true);
		_mat.set(multidim_index[0], multidim_index[1], s);
	}
	void set_to_zero(size_t index_in_TStorageYMatrix) {
		const auto multidim_index = get_multi_dimensional_index(index_in_TStorageYMatrix);
		auto s                    = _mat.get(multidim_index[0], multidim_index[1]);
		s.set_state(false);
		_mat.set(multidim_index[0], multidim_index[1], s);
	}

	void insert_one(uint64_t linear_index_in_Y_space) { _insert(linear_index_in_Y_space, true); }

	/// Does the same as set_to_one but sets the element to zero
	void insert_zero(uint64_t linear_index_in_Y_space) { _insert(linear_index_in_Y_space, false); }

	void add_to_counter(size_t iteration) {
		if (iteration % _thinning_factor == 0) {
			_mat.updateValues([](TStorageY &elem) { elem.update_counter(); });
		}
	}

	[[nodiscard]] double get_fraction_of_ones(uint64_t linear_index_in_Y_space) const {
		const auto multidim_index = get_multi_dimensional_index(linear_index_in_Y_space);
		// A missing cell returns a default-constructed TStorageY (counter == 0),
		// so the result is 0.0 for absent entries without an explicit empty check.
		const auto storage        = _mat.get(multidim_index[0], multidim_index[1]);
		return static_cast<double>(storage.get_counter()) / static_cast<double>(_total_counts);
	}

	size_t get_total_counts() const { return _total_counts; }

	TStorageY operator[](size_t index_in_TStorageYMatrix) const {
		const auto multidim_index = get_multi_dimensional_index(index_in_TStorageYMatrix);
		return _mat.get(multidim_index[0], multidim_index[1]);
	}

	TStorageY
	operator[](const std::array<size_t, NUMBER_OF_TREES> &index_in_TStorageYMatrix) const {
		return _mat.get(index_in_TStorageYMatrix[0], index_in_TStorageYMatrix[1]);
	}

	void reset_counts() {
		_mat.updateValues([](TStorageY &val) { val.reset_counter(); });
	}

	/// Remove all the elements that have the state to zero.
	/// @return void
	void remove_zeros() {
		_mat.erase_if([](const TStorageY &elem) { return !elem.is_one(); });
	}

	/// Given a multi-dimensional index, we want to get its linear index.
	/// @param multi_dim_index the multi-dimensional index
	/// @return the linear index
	uint64_t
	get_linear_index_in_Y_space(const std::vector<size_t> &multidim_index_in_Y_space) const {
		return coretools::getLinearIndex(multidim_index_in_Y_space, _dimensions_Y_space);
	}
	uint64_t get_linear_index_in_container_space(
	    const std::vector<size_t> &multidim_index_in_Y_space) const {
		return get_linear_index_in_Y_space(multidim_index_in_Y_space);
	}

	/// Given a linear index, we want to get the multi-dimensional index.
	[[nodiscard]] std::array<size_t, NUMBER_OF_TREES>
	get_multi_dimensional_index(uint64_t linear_index_in_Y_space) const {
		const auto tmp = static_cast<size_t>(linear_index_in_Y_space);
		return coretools::getSubscriptsAsArray(tmp, _dimensions_Y_space);
	}

	/// Given a linear index, we want to get the multi-dimensional index.
	[[nodiscard]] std::array<size_t, NUMBER_OF_TREES>
	get_multi_dimensional_index(size_t linear_index_in_Y_space) const {
		return coretools::getSubscriptsAsArray(linear_index_in_Y_space, _dimensions_Y_space);
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

	std::vector<uint8_t> get_full_Y_binary_vector() const {
		std::vector<uint8_t> Y_as_vector;
		size_t total_size = total_size_of_container_space();
		Y_as_vector.reserve(total_size);
		for (size_t i = 0; i < _mat.nRows(); ++i) {
			const auto row = _mat.getRow(i);
			for (const auto &val : row) { Y_as_vector.push_back(val.is_one()); }
		}
		return Y_as_vector;
	}

	void insert_in_Y(const std::vector<std::vector<size_t>> &linear_indices_in_Y_space_to_insert) {
		// Step 1: Merge all vectors into a single sorted vector
		std::vector<size_t> merged_vec;
		for (const auto &vec : linear_indices_in_Y_space_to_insert) {
			merged_vec.insert(merged_vec.end(), vec.begin(), vec.end());
		}

		for (const auto &it : merged_vec) {
			const auto multi_dim_index = get_multi_dimensional_index(it);
			_mat.setRaw(multi_dim_index[0], multi_dim_index[1], TStorageY(true));
		}
		_mat.cleanUp();
	}

	[[nodiscard]] bool empty() const { return _mat.nNonZero() == 0; }
	[[nodiscard]] size_t number_of_dimensions() const { return _dimensions_Y_space.size(); }
	[[nodiscard]] size_t get_thinning_factor() const { return _thinning_factor; }
	[[nodiscard]] size_t number_of_ones() const {
		size_t count = 0;
		for (size_t row = 0; row < _mat.nRows(); ++row) {
			for (auto it = _mat.begin_row(row); it != _mat.end_row(row); ++it) {
				if (it->val.is_one()) { ++count; }
			}
		}
		return count;
	}
};
