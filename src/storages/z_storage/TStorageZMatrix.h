//
// Created by VISANI Marco on 18.10.2024.
//

#ifndef TStorageZMatrix_H
#define TStorageZMatrix_H

#include "TStorageZ.h"
#include "constants.h"
#include "coretools/Main/TError.h"
#include "coretools/Math/TSparseMatrix.h"
#include "coretools/algorithms.h"
#include "storages/TSparseWindow.h"
#include "storages/storage_concepts.h"
#include <cstddef>
#include <utility>
#include <vector>

/// There is one TStorageZMatrix per dimension `d`. The dimensions in Z space are the number of
/// leaves in each dimension except for `d`, which spans every node of its tree -- leaves included
/// (ADR-0005). For example, with leaf counts [2, 3] and 19 nodes in dimension `d == 0`, the
/// dimensions in Z space are [19, 3], where the first two rows are the leaves the field also
/// holds. `node_state_dimensions` in tree/node_state_shape.h is where that rule lives.
///
/// Backed by a coretools::TSparseMatrix<TStorageZ>, mirroring TStorageYMatrix: rows and
/// columns are kept sorted, so a clique's nodes can be range-walked in O(nnz in that
/// line) and the (row, col) position encodes the linear index in Z space.
class TStorageZMatrix {
private:
	IndexArray _dimensions_in_Z_space{};
	coretools::TSparseMatrix<TStorageZ> _mat;

	/// Allocation-free (row, col) of a linear index, for the hot internal paths.
	[[nodiscard]] IndexArray _row_col(size_t linear_index_in_Z_space) const {
		return coretools::getSubscriptsAsArray(linear_index_in_Z_space, _dimensions_in_Z_space);
	}

	void _insert(size_t linear_index_in_Z_space, bool state) {
		if (linear_index_in_Z_space >= this->total_size_of_container_space()) {
			throw coretools::TDevError(
			    "You are trying to insert a value at an linear index bigger than the total "
			    "size of the container. ",
			    "The index is: ", linear_index_in_Z_space,
			    " and the total size of the container is : ",
			    this->total_size_of_container_space());
		}
		const auto md = _row_col(linear_index_in_Z_space);
		_mat.set(md[0], md[1], TStorageZ(state));
	}

public:
	TStorageZMatrix()  = default;
	~TStorageZMatrix() = default;
	explicit TStorageZMatrix(const IndexArray &dimensions_in_Z_space) {
		initialize_dimensions(dimensions_in_Z_space);
	}

	void initialize_dimensions(const IndexArray &dimensions_in_Z_space) {
		_dimensions_in_Z_space = dimensions_in_Z_space;
		_mat.resize(dimensions_in_Z_space[0], dimensions_in_Z_space[1]);
	}

	/// Point lookup by linear index in Z space. A missing cell reads as false.
	[[nodiscard]] inline bool is_one(size_t linear_index_in_Z_space) const {
		const auto md = _row_col(linear_index_in_Z_space);
		return _mat.get(md[0], md[1]).is_one();
	}

	/// Set the state of the cell at `linear_index_in_Z_space` in place. If the cell does not exist
	/// yet, TSparseMatrix::set inserts it.
	void set_state(size_t linear_index_in_Z_space, bool state) {
		const auto md = _row_col(linear_index_in_Z_space);
		auto s        = _mat.get(md[0], md[1]);
		s.set_state(state);
		_mat.set(md[0], md[1], s);
	}

	void insert_one(size_t linear_index_in_Z_space) { _insert(linear_index_in_Z_space, true); }
	void insert_one(const IndexArray &multi_dim_index_in_Z_space) {
		insert_one(get_linear_index_in_Z_space(multi_dim_index_in_Z_space));
	}
	void insert_zero(size_t linear_index_in_Z_space) { _insert(linear_index_in_Z_space, false); }

	/// Remove all the elements whose state is zero.
	void remove_zeros() {
		_mat.erase_if([](const TStorageZ &elem) { return !elem.is_one(); });
	}

	[[nodiscard]] size_t
	get_linear_index_in_Z_space(const IndexArray &multidim_index_in_Z_space) const {
		return coretools::getLinearIndex(multidim_index_in_Z_space, _dimensions_in_Z_space);
	}
	[[nodiscard]] size_t
	get_linear_index_in_container_space(const IndexArray &multidim_index_in_Z_space) const {
		return get_linear_index_in_Z_space(multidim_index_in_Z_space);
	}

	[[nodiscard]] IndexArray get_multi_dimensional_index(size_t linear_index_in_Z_space) const {
		return _row_col(linear_index_in_Z_space);
	}

	[[nodiscard]] size_t total_size_of_container_space() const {
		return coretools::containerProduct(_dimensions_in_Z_space);
	}

	using TWindow = TSparseWindow<TStorageZ>;

	/// A strided window over `n_cells` cells, the first at `start_index` and the rest `stride`
	/// apart. The window walks the line once on open, and buffers an insert until it closes.
	[[nodiscard]] TWindow open_window(const IndexArray &start_index, size_t n_cells,
	                                  size_t stride) {
		const size_t start_linear = coretools::getLinearIndex(start_index, _dimensions_in_Z_space);
		DEBUG_ASSERT(n_cells == 0 ||
		             start_linear + (n_cells - 1) * stride < total_size_of_container_space());
		return {_mat, _dimensions_in_Z_space[1], start_linear, n_cells, stride};
	}

	/// Bulk-insert deferred 0 -> 1 transitions (linear indices in Z space), then re-sort once.
	/// Mirror of TStorageYMatrix::insert_in_Y.
	void insert_in_Z(const std::vector<std::vector<size_t>> &linear_indices_in_Z_space_to_insert) {
		std::vector<size_t> merged_vec;
		for (const auto &vec : linear_indices_in_Z_space_to_insert) {
			merged_vec.insert(merged_vec.end(), vec.begin(), vec.end());
		}
		for (const auto &it : merged_vec) {
			const auto md = _row_col(it);
			_mat.setRaw(md[0], md[1], TStorageZ(true));
		}
		_mat.cleanUp();
	}

	[[nodiscard]] std::vector<size_t> get_full_Z_binary_vector() const {
		std::vector<size_t> Z_as_vector;
		Z_as_vector.reserve(total_size_of_container_space());
		for (size_t i = 0; i < _mat.nRows(); ++i) {
			const auto row = _mat.getRow(i);
			for (const auto &val : row) { Z_as_vector.push_back(val.is_one()); }
		}
		return Z_as_vector;
	}

	/// Returns every stored cell as (linear index in Z space, value), in ascending linear-index
	/// order (rows then columns is row-major, matching linear = row * nCols + col).
	[[nodiscard]] std::vector<std::pair<size_t, TStorageZ>> get_stored_entries() const {
		std::vector<std::pair<size_t, TStorageZ>> entries;
		entries.reserve(_mat.nNonZero());
		const size_t n_cols = _dimensions_in_Z_space[1];
		for (size_t row = 0; row < _mat.nRows(); ++row) {
			for (auto it = _mat.begin_row(row); it != _mat.end_row(row); ++it) {
				entries.emplace_back(row * n_cols + it->index, it->val);
			}
		}
		return entries;
	}

	[[nodiscard]] bool empty() const { return _mat.nNonZero() == 0; }
	[[nodiscard]] size_t size() const { return _mat.nNonZero(); }
};

static_assert(BinaryFieldStorage<TStorageZMatrix>,
              "The sparse node state must satisfy the binary storage interface.");

#endif // TStorageZMatrix_H
