//
// Created by VISANI Marco on 17.10.2024.
//

#pragma once

#include "TStorageY.h"
#include "constants.h"
#include "coretools/Main/TError.h"
#include "coretools/Math/TSparseMatrix.h"
#include "coretools/algorithms.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

class TStorageYMatrix {
private:
	size_t _thinning_factor;
	size_t _total_counts;
	coretools::TSparseMatrix<TStorageY> _mat;

	/// _dimensions_Y_space is the number of leaf nodes in each dimension
	IndexArray _dimensions_Y_space;

	/// Allocation-free (row, col) of a linear index, for the hot internal paths. The public
	/// get_multi_dimensional_index returns a std::vector and is meant for external callers.
	[[nodiscard]] IndexArray _row_col(size_t linear_index_in_Y_space) const {
		return coretools::getSubscriptsAsArray(linear_index_in_Y_space, _dimensions_Y_space);
	}

	void _insert(size_t linear_index_in_Y_space, bool state) {
		if (linear_index_in_Y_space >= this->total_size_of_container_space()) {
			throw coretools::TDevError(
			    "You are trying to insert a value at an linear index bigger than the total "
			    "size of the container. ",
			    "The index is: ", linear_index_in_Y_space,
			    " and the total size of the container is : ",
			    this->total_size_of_container_space());
		}
		const auto multidim_index = _row_col(linear_index_in_Y_space);
		_mat.set(multidim_index[0], multidim_index[1], TStorageY(state));
	};

public:
	TStorageYMatrix()  = default;
	~TStorageYMatrix() = default;
	TStorageYMatrix(const size_t n_iterations, const IndexArray &dimensions_Y_space) {
		initialize(n_iterations, dimensions_Y_space);
	};

	void initialize(const size_t n_iterations, const IndexArray &dimensions_Y_space) {
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
		IndexArray dims = {dimensions_Y_space[0], dimensions_Y_space[1]};
		initialize(n_iterations, dims);
	}

	/// We want to check if element at position index_in_TStorageYMatrix is one.
	/// This index *must* be previously obtained by a binary search.
	/// @param index_in_TStorageYMatrix the position of the element in the Y vector.
	/// @return true if the element is one, false otherwise.
	[[nodiscard]] inline bool is_one(const size_t index_in_TStorageYMatrix) const {
		const auto multidim_index = _row_col(index_in_TStorageYMatrix);
		return _mat.get(multidim_index[0], multidim_index[1]).is_one();
	};

	[[nodiscard]] inline bool is_one(const IndexArray &index_in_Y_space) const {
		return _mat.get(index_in_Y_space[0], index_in_Y_space[1]).is_one();
	};

	/** set_to_one will set the element at the index_in_TStorageYMatrix to 1.
	 * @param index_in_TStorageYMatrix the position of the element in the Y vector
	 */
	void set_to_one(size_t index_in_TStorageYMatrix) {
		const auto multidim_index = _row_col(index_in_TStorageYMatrix);
		auto s                    = _mat.get(multidim_index[0], multidim_index[1]);
		s.set_state(true);
		_mat.set(multidim_index[0], multidim_index[1], s);
	}
	void set_to_zero(size_t index_in_TStorageYMatrix) {
		const auto multidim_index = _row_col(index_in_TStorageYMatrix);
		auto s                    = _mat.get(multidim_index[0], multidim_index[1]);
		s.set_state(false);
		_mat.set(multidim_index[0], multidim_index[1], s);
	}

	/// Set the state of the cell at `linear_index_in_Y_space` to `state`, preserving its counter.
	/// If the cell does not exist yet, TSparseMatrix::set inserts it (counter starts at 0), so this
	/// unifies the old set_to_one / set_to_zero / insert-later branches into a single call.
	void set_state(size_t linear_index_in_Y_space, bool state) {
		const auto multidim_index = _row_col(linear_index_in_Y_space);
		auto s                    = _mat.get(multidim_index[0], multidim_index[1]);
		s.set_state(state);
		_mat.set(multidim_index[0], multidim_index[1], s);
	}

	void insert_one(size_t linear_index_in_Y_space) { _insert(linear_index_in_Y_space, true); }

	/// Does the same as set_to_one but sets the element to zero
	void insert_zero(size_t linear_index_in_Y_space) { _insert(linear_index_in_Y_space, false); }

	void add_to_counter(size_t iteration) {
		if (iteration % _thinning_factor == 0) {
			_mat.updateValues([](TStorageY &elem) { elem.update_counter(); });
		}
	}

	[[nodiscard]] double get_fraction_of_ones(size_t linear_index_in_Y_space) const {
		const auto multidim_index = _row_col(linear_index_in_Y_space);
		// A missing cell returns a default-constructed TStorageY (counter == 0),
		// so the result is 0.0 for absent entries without an explicit empty check.
		const auto storage        = _mat.get(multidim_index[0], multidim_index[1]);
		return static_cast<double>(storage.get_counter()) / static_cast<double>(_total_counts);
	}

	[[nodiscard]] size_t get_total_counts() const { return _total_counts; }

	TStorageY operator[](size_t index_in_TStorageYMatrix) const {
		const auto multidim_index = _row_col(index_in_TStorageYMatrix);
		return _mat.get(multidim_index[0], multidim_index[1]);
	}

	TStorageY operator[](const IndexArray &index_in_TStorageYMatrix) const {
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
	[[nodiscard]] size_t
	get_linear_index_in_Y_space(const IndexArray &multidim_index_in_Y_space) const {
		return coretools::getLinearIndex(multidim_index_in_Y_space, _dimensions_Y_space);
	}
	[[nodiscard]] size_t
	get_linear_index_in_container_space(const IndexArray &multidim_index_in_Y_space) const {
		return get_linear_index_in_Y_space(multidim_index_in_Y_space);
	}

	/// Given a linear index, we want to get the multi-dimensional index (in Y / leaves space).
	/// Returns a std::vector for callers that need it (get_clique, x_is_one, ...); the hot internal
	/// paths use the allocation-free _row_col instead.
	[[nodiscard]] IndexArray get_multi_dimensional_index(size_t linear_index_in_Y_space) const {
		return _row_col(linear_index_in_Y_space);
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

	/// Fast current-state fill for a clique of `K` nodes running along one dimension, starting at
	/// `start_index` (a multi-dimensional index in leaves space). The clique's nodes occupy the
	/// contiguous coordinates [start, start + K) along the variable dimension, all other
	/// coordinates fixed. Outputs, for every k in [0, K):
	///   - current_state[k] : whether that cell is currently a one,
	///   - exists[k]        : whether that cell is stored in the matrix,
	///   - linear_index[k]  : the linear index in Y space of that cell.
	///
	/// This replaces the binary-search machinery of fill_current_state for Y: because the matrix
	/// keeps rows and columns sorted, we walk the relevant row/column once (O(nnz in that line)):
	///   - increment == 1   : the variable dimension is the last one -> walk the matrix row.
	///   - increment  > 1   : the variable dimension is the first one -> walk the matrix column.
	void fill_current_state(const IndexArray &start_index, size_t K, size_t increment,
	                        std::vector<uint8_t> &current_state, std::vector<uint8_t> &exists,
	                        std::vector<size_t> &linear_index) const {
		current_state.assign(K, 0);
		exists.assign(K, 0);
		linear_index.assign(K, 0);

		const size_t start_linear = coretools::getLinearIndex(start_index, _dimensions_Y_space);
		for (size_t k = 0; k < K; ++k) { linear_index[k] = start_linear + k * increment; }

		if (increment == 1) {
			// variable dimension is the last one -> a single matrix row, entries sorted by column.
			const size_t row       = start_index[0];
			const size_t start_col = start_index[1];
			for (auto it = _mat.begin_row(row); it != _mat.end_row(row); ++it) {
				if (it->index < start_col) { continue; }
				const size_t k = it->index - start_col;
				if (k >= K) { break; } // sorted -> no later entry can fall in range
				current_state[k] = it->val.is_one();
				exists[k]        = 1;
			}
		} else {
			// variable dimension is the first one -> a single matrix column, entries sorted by row.
			const size_t col       = start_index[1];
			const size_t start_row = start_index[0];
			for (auto it = _mat.begin_col(col); it != _mat.end_col(col); ++it) {
				if (it->index < start_row) { continue; }
				const size_t k = it->index - start_row;
				if (k >= K) { break; }
				current_state[k] = it->val.is_one();
				exists[k]        = 1;
			}
		}
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
			const auto multi_dim_index = _row_col(it);
			_mat.setRaw(multi_dim_index[0], multi_dim_index[1], TStorageY(true));
		}
		_mat.cleanUp();
	}

	[[nodiscard]] bool empty() const { return _mat.nNonZero() == 0; }
	[[nodiscard]] double sparsity() const {
		return static_cast<double>(_mat.nNonZero()) /
		       static_cast<double>(total_size_of_container_space());
	}
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

	/// Returns every stored cell as (linear index in Y space, value), in ascending linear-index
	/// order. Iterating rows then columns is naturally row-major, which matches the linear-index
	/// layout (linear = row * nCols + col). Used to write only the non-default entries of Y.
	[[nodiscard]] std::vector<std::pair<size_t, TStorageY>> get_stored_entries() const {
		std::vector<std::pair<size_t, TStorageY>> entries;
		entries.reserve(_mat.nNonZero());
		const size_t n_cols = _dimensions_Y_space[1];
		for (size_t row = 0; row < _mat.nRows(); ++row) {
			for (auto it = _mat.begin_row(row); it != _mat.end_row(row); ++it) {
				entries.emplace_back(row * n_cols + it->index, it->val);
			}
		}
		return entries;
	}

	/// Allocation-free forward walk over the stored cells in ascending linear-index order
	/// (row-major: linear = row * nCols + col) -- the streaming form of get_stored_entries().
	/// Within a row the sparse entries are already column-sorted and rows are visited in order,
	/// so the linear index increases monotonically. Lets callers merge-join two matrices without
	/// materializing a vector for either.
	class StoredCursor {
		const coretools::TSparseMatrix<TStorageY> *_mat = nullptr;
		size_t _n_cols                                  = 0;
		size_t _row                                     = 0;
		decltype(std::declval<const coretools::TSparseMatrix<TStorageY>>().begin_row(0)) _it{};

		void _skip_empty_rows() {
			while (_row < _mat->nRows() && _it == _mat->end_row(_row)) {
				++_row;
				if (_row < _mat->nRows()) { _it = _mat->begin_row(_row); }
			}
		}

	public:
		StoredCursor() = default;
		StoredCursor(const coretools::TSparseMatrix<TStorageY> &mat, size_t n_cols)
		    : _mat(&mat), _n_cols(n_cols) {
			if (_mat->nRows() > 0) {
				_it = _mat->begin_row(0);
				_skip_empty_rows();
			}
		}

		[[nodiscard]] bool valid() const { return _mat != nullptr && _row < _mat->nRows(); }
		[[nodiscard]] size_t linear_index() const { return _row * _n_cols + _it->index; }
		[[nodiscard]] bool is_one() const { return _it->val.is_one(); }
		[[nodiscard]] const TStorageY &value() const { return _it->val; }
		void advance() {
			++_it;
			_skip_empty_rows();
		}
	};

	[[nodiscard]] StoredCursor stored_cursor() const {
		return StoredCursor(_mat, _dimensions_Y_space[1]);
	}
};
