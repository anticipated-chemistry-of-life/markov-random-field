//
// The simple error model's pass through one block update (data_sources/data_source_pass.h).
//
// What the pass says about a cell is P(D_cell | Y) at the current error rate. What it tallies is
// the number of cells D contradicts, the one statistic the epsilon move reads
// (TSimpleErrorModelMath.h). The pass is built from the data and the error rate alone, so a test
// builds one without a tree or stattools; `TSimpleErrorModel` hands one out and takes it back
// (`begin_block_update`, `end_block_update`), and nothing else installs a disagreement count.
//
// Not guarded by USE_SIMPLE_ERROR_MODEL, like TSimpleErrorModelMath.h, so it compiles and is unit
// tested in every build configuration.
//

#pragma once

#include "coretools/Main/TError.h"
#include "data_sources/data_source_pass.h"
#include "simple_error_model/TSimpleErrorModelMath.h"
#include "storages/storage_backend.h"
#include "storages/storage_concepts.h"
#include <array>
#include <atomic>
#include <cstddef>

namespace simple_error_model {

class TBlockPass {
private:
	/// The observed data, with the field's dimensions. The pass only reads it.
	const TBinaryStorage *_D = nullptr;
	/// The error rate the update draws under. It does not move during the update.
	double _eps              = 0.0;
	/// The cells D contradicts, over every row that has ended so far.
	size_t _n_disagree       = 0;

public:
	/// What one row carries from cell to cell.
	struct TRow {
		/// The linear index of the row's cell at molecule leaf 0. D has the field's dimensions and
		/// the field is row-major, so a cell's index in D is this plus its molecule leaf.
		size_t first_field_cell = 0;
		/// The state D holds for the cell `cell_factor` was last asked about, kept so that `record`
		/// does not look it up a second time.
		bool observed           = false;
		/// The cells of this row D contradicts so far.
		size_t n_disagree       = 0;
	};

	TBlockPass(const TBinaryStorage &D, double eps) : _D(&D), _eps(eps) {}

	/// Nothing to ready: D is not written during the update.
	void begin_pass() {}

	[[nodiscard]] TRow begin_row([[maybe_unused]] size_t species_leaf,
	                             size_t first_field_cell) const {
		DEBUG_ASSERT(_D->get_linear_index_in_container_space(IndexArray{species_leaf, 0}) ==
		             first_field_cell);
		return TRow{.first_field_cell = first_field_cell, .observed = false, .n_disagree = 0};
	}

	/// {P(D_cell | Y = 0), P(D_cell | Y = 1)}.
	[[nodiscard]] std::array<double, 2> cell_factor(TRow &row, size_t molecule_leaf) const {
		row.observed = _D->is_one(row.first_field_cell + molecule_leaf);
		std::array<double, 2> prob{};
		probabilities_for_both_Y_states(row.observed, _eps, prob);
		return prob;
	}

	/// The cell disagrees when D holds the other state. The fold has asserted that this is the cell
	/// `cell_factor` was last asked about (field/TBlockData.h).
	void record(TRow &row, size_t /*molecule_leaf*/, bool drawn_y) const {
		row.n_disagree += static_cast<size_t>(row.observed != drawn_y);
	}

	/// Rows end on many threads at once, and each adds its count atomically. The count is an
	/// integer, so the sum is the same whichever rows got there first, and so at any thread count.
	/// That is this pass's rule to keep, not the block update's (ADR-0017).
	void end_row(const TRow &row) {
		std::atomic_ref<size_t>(_n_disagree).fetch_add(row.n_disagree, std::memory_order_relaxed);
	}

	/// The update visited every cell once, so the count is the field's. A debug build recounts it
	/// from the field the update left; the recount walks the field, so a release build does not.
	template<FieldStorage Field> void end_pass([[maybe_unused]] const Field &Y) const {
#ifndef NDEBUG
		DEBUG_ASSERT(_n_disagree == count_disagreements(Y, *_D));
#endif
	}

	/// The cells D contradicts in the field the update left, once the pass has ended.
	[[nodiscard]] size_t n_disagree() const noexcept { return _n_disagree; }
};

static_assert(DataSourcePass<TBlockPass>,
              "The simple error model's pass must answer what the block update asks a source.");

} // namespace simple_error_model
