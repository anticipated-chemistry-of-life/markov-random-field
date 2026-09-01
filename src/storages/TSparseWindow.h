//
// The window the sparse storages open over themselves.
//

#pragma once

#include "coretools/Main/TError.h"
#include "coretools/Math/TSparseMatrix.h"
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

/// A strided window over a sparse matrix of binary cells: the matrix, the linear index the window
/// starts at, how many cells it holds, and how far apart they are.
///
/// The window materialises its line once, on open, exactly as the clique fill does. A point lookup
/// in this matrix costs a search, so the sparse path pays that search once per line. See ADR-0006.
///
/// A write to a cell the matrix already holds goes in place. A write to one it does not hold waits
/// in a buffer, because the insert would reallocate a row and a column while other threads read
/// them.
///
/// A window ends once, in one of two ways. `close` commits the buffer into the matrix.
/// `take_buffered_inserts` hands the buffer to the caller and writes nothing, which is the only
/// exit a window inside a parallel region may take. See ADR-0006.
///
/// A window shows its own buffered write to a later read on the same window. That is what lets a
/// post-order walk read the states it has just assigned, and what keeps the sparse and the dense
/// backend on one chain inside a single update.
///
/// `Cell` is the value the matrix holds. It answers `is_one`, takes `set_state`, and is built from
/// a state -- which the field's cell and the node state's cell both do.
template<typename Cell> class TSparseWindow {
private:
	/// Where a cell of this window is. `held` is in the matrix and is written in place. `buffered`
	/// is a write the matrix cannot take yet, and waits for `close`.
	enum class THolding : uint8_t { absent, held, buffered };

	coretools::TSparseMatrix<Cell> *_mat = nullptr;
	size_t _n_cols                       = 0;
	size_t _start_linear                 = 0;
	size_t _stride                       = 1;

	/// The state of every cell of the window: what the matrix held when the window opened, plus
	/// every write since. A read answers from here, which is what makes a buffered write visible.
	std::vector<uint8_t> _states;
	/// Where each cell is.
	std::vector<THolding> _holding;
	/// The cells whose insert waits for `close`.
	std::vector<size_t> _to_insert;
	bool _open = false;

	/// The (row, column) of the window's `k`th cell.
	[[nodiscard]] std::pair<size_t, size_t> _row_col(size_t k) const {
		const size_t linear = linear_index(k);
		return {linear / _n_cols, linear % _n_cols};
	}

	/// Reads the window's cells out of the matrix, in one walk of one line.
	void _materialise() {
		const size_t n_cells = _states.size();
		if (n_cells == 0) { return; }

		const size_t start_row = _start_linear / _n_cols;
		const size_t start_col = _start_linear % _n_cols;

		// The column count is part of the test because a single-column container gives a window
		// along either dimension a stride of one, and only the first dimension can be longer than
		// one cell there. A container is one column wide whenever the other tree has a single leaf.
		if (_stride == 1 && _n_cols > 1) {
			// the variable dimension is the last one -> a single matrix row, sorted by column
			for (auto it = _mat->begin_row(start_row); it != _mat->end_row(start_row); ++it) {
				if (it->index < start_col) { continue; }
				const size_t k = it->index - start_col;
				if (k >= n_cells) { break; } // sorted -> no later entry can fall in range
				_states[k]  = it->val.is_one();
				_holding[k] = THolding::held;
			}
		} else {
			// the variable dimension is the first one -> a single matrix column, sorted by row
			for (auto it = _mat->begin_col(start_col); it != _mat->end_col(start_col); ++it) {
				if (it->index < start_row) { continue; }
				const size_t k = it->index - start_row;
				if (k >= n_cells) { break; }
				_states[k]  = it->val.is_one();
				_holding[k] = THolding::held;
			}
		}
	}

	/// Writes the state of a cell the matrix holds, and keeps whatever else that cell carries.
	/// The field's cell carries a posterior counter, which a state write must not disturb.
	void _write_in_place(size_t k, bool state) {
		const auto [row, col] = _row_col(k);
		Cell cell             = _mat->get(row, col);
		cell.set_state(state);
		_mat->set(row, col, cell);
	}

public:
	TSparseWindow(coretools::TSparseMatrix<Cell> &mat, size_t n_cols, size_t start_linear,
	              size_t n_cells, size_t stride)
	    : _mat(&mat), _n_cols(n_cols), _start_linear(start_linear), _stride(stride),
	      _states(n_cells, 0), _holding(n_cells, THolding::absent), _open(true) {
		// A window runs along one dimension of a two-dimensional container, so it steps by one
		// cell or by one row. The line walk above reads the stride that way.
		DEBUG_ASSERT(stride == 1 || stride == n_cols);
		_materialise();
	}

	/// Commits the buffer, for a window the caller lets go of rather than closing.
	~TSparseWindow() { close(); }

	// A window owns writes that are not in the matrix yet, so it is neither copied nor moved.
	TSparseWindow(const TSparseWindow &)            = delete;
	TSparseWindow &operator=(const TSparseWindow &) = delete;
	TSparseWindow(TSparseWindow &&)                 = delete;
	TSparseWindow &operator=(TSparseWindow &&)      = delete;

	[[nodiscard]] size_t size() const { return _states.size(); }

	/// The linear index in container space of the window's `k`th cell.
	[[nodiscard]] size_t linear_index(size_t k) const {
		DEBUG_ASSERT(k < _states.size());
		return _start_linear + k * _stride;
	}

	[[nodiscard]] bool is_one(size_t k) const {
		DEBUG_ASSERT(k < _states.size());
		return _states[k] != 0;
	}

	void set_state(size_t k, bool state) {
		DEBUG_ASSERT(k < _states.size());
		DEBUG_ASSERT(_open);
		_states[k] = state ? 1 : 0;
		switch (_holding[k]) {
		case THolding::held: _write_in_place(k, state); break;
		case THolding::absent:
			_holding[k] = THolding::buffered;
			_to_insert.push_back(k);
			break;
		case THolding::buffered: break; // already buffered; the state above is what it commits
		}
	}

	/// Hands out the cells this window could not write in place, as linear indices in container
	/// space, and closes the window. The caller commits them, which is what lets a window end
	/// inside a parallel region. An insert writes one row and one column of the matrix, and
	/// threads share one of the two however the work is split. See ADR-0006.
	///
	/// A cell that ended the window as a zero is left out, because a cell the matrix does not hold
	/// already reads as zero.
	[[nodiscard]] std::vector<size_t> take_buffered_inserts() {
		std::vector<size_t> linear_indices;
		linear_indices.reserve(_to_insert.size());
		for (const size_t k : _to_insert) {
			if (_states[k] == 0) { continue; }
			linear_indices.push_back(linear_index(k));
		}
		_to_insert.clear();
		_open = false;
		return linear_indices;
	}

	/// Inserts every buffered cell that ended the window as a one. A cell that ended as a zero is
	/// dropped, because a cell the matrix does not hold already reads as zero.
	///
	/// This inserts, so it reallocates rows and columns. It must run outside a parallel region. A
	/// window that ran inside one hands its buffer out instead.
	void close() {
		if (!_open) { return; }
		_open = false;
		for (const size_t k : _to_insert) {
			if (_states[k] == 0) { continue; }
			const auto [row, col] = _row_col(k);
			_mat->set(row, col, Cell(true));
			_holding[k] = THolding::held;
		}
		_to_insert.clear();
	}
};
