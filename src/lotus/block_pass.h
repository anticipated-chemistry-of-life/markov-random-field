//
// LOTUS's pass through one block update (data_sources/data_source_pass.h).
//
// The pass says P(L_cell | Y) about a cell, under the current reporting model. It tallies the
// field's present cells by kind. The gamma and error-rate moves read that tally (TLotusMath.h,
// ADR-0015). A test builds the pass from the records, the paper-count classes and the reporting
// model alone, without a tree or stattools. `TLotus` hands a pass out and takes it back
// (`begin_block_update`, `end_block_update`). Nothing else installs present cells.
//
// Not guarded by USE_LOTUS, like TLotusMath.h, so it compiles and is unit tested in every build
// configuration.
//

#pragma once

#include "coretools/Main/TError.h"
#include "data_sources/data_source_pass.h"
#include "lotus/TLotusMath.h"
#include "storages/TSparse.h"
#include "storages/storage_concepts.h"
#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace lotus_math {

class TBlockPass {
private:
	/// The LOTUS records, with the field's dimensions. The pass reads them, and refreshes their
	/// sorted-ones cache.
	const TSparseBinary *_L            = nullptr;
	/// The paper-count class of every leaf: which kind a present cell is counted as.
	const TPaperCountClasses *_classes = nullptr;
	/// The reporting model the update draws under. It does not move during the update.
	const TReportingModel *_reporting  = nullptr;
	/// The present cells, over every row that has ended so far.
	TPresentCells _present_cells;

public:
	/// What one row carries from cell to cell.
	struct TRow {
		/// A forward cursor into the records' sorted-ones cache, seeded at the row's first cell.
		/// The row's cells arrive in ascending order, which is all the cursor asks.
		TSparseBinary::OnesCursor records;
		/// The linear index of the row's cell at molecule leaf 0. L has the field's dimensions and
		/// the field is row-major, so a cell's index in L is this plus its molecule leaf.
		size_t first_field_cell = 0;
		/// The species leaf the row is, which research effort reads.
		size_t species_leaf     = 0;
		/// The paper-count class of the row's species leaf: the row of kinds it counts into.
		uint32_t species_class  = 0;
		/// Whether LOTUS reports the cell `cell_factor` was last asked about, kept so that `record`
		/// does not ask the cursor a second time.
		bool reported           = false;
		/// The row's present cells, `[reported * n_molecule_classes + molecule_class]`. A row has
		/// at most a molecule leaf's worth, so a 32-bit count holds it.
		std::vector<uint32_t> present;
	};

	TBlockPass(const TSparseBinary &L, const TPaperCountClasses &classes,
	           const TReportingModel &reporting)
	    : _L(&L), _classes(&classes), _reporting(&reporting),
	      _present_cells(L.total_size_of_container_space(), L.number_of_ones(), classes) {}

	/// Makes the records' sorted-ones cache fresh, so that every row may seed a cursor of its own
	/// from inside the parallel region. The sort is what is not safe to do concurrently.
	void begin_pass() { _L->refresh_ones(); }

	[[nodiscard]] TRow begin_row(size_t species_leaf, size_t first_field_cell) const {
		DEBUG_ASSERT(_L->get_linear_index_in_container_space(IndexArray{species_leaf, 0}) ==
		             first_field_cell);
		return TRow{.records          = _L->fresh_ones_cursor_from(first_field_cell),
		            .first_field_cell = first_field_cell,
		            .species_leaf     = species_leaf,
		            .species_class    = _classes->class_of(0, species_leaf),
		            .reported         = false,
		            .present = std::vector<uint32_t>(2 * _present_cells.n_molecule_classes)};
	}

	/// {P(L_cell | Y = 0), P(L_cell | Y = 1)}. The one place the cursor is asked about a cell.
	[[nodiscard]] std::array<double, 2> cell_factor(TRow &row, size_t molecule_leaf) const {
		row.reported = row.records.advance_to_and_check(row.first_field_cell + molecule_leaf);
		const IndexArray cell{row.species_leaf, molecule_leaf};
		return {_reporting->probability(false, row.reported, cell),
		        _reporting->probability(true, row.reported, cell)};
	}

	/// Counts a present cell as its kind: its molecule class, and whether LOTUS reports it. An
	/// absent cell's term does not depend on where it is, so it is not counted. The fold has
	/// asserted that this is the cell `cell_factor` was last asked about (field/TBlockData.h).
	void record(TRow &row, size_t molecule_leaf, bool drawn_y) const {
		if (!drawn_y) { return; }
		const size_t reported = static_cast<size_t>(row.reported);
		++row.present[reported * _present_cells.n_molecule_classes +
		              _classes->class_of(1, molecule_leaf)];
	}

	/// Rows end on many threads at once, and each adds its row of counts entry by entry,
	/// atomically. The counts are integers, so they are the same whichever rows got there first,
	/// and so at any thread count (ADR-0015, ADR-0017).
	void end_row(const TRow &row) {
		const size_t n_molecule_classes = _present_cells.n_molecule_classes;
		for (size_t reported = 0; reported < 2; ++reported) {
			for (uint32_t m = 0; m < n_molecule_classes; ++m) {
				const uint32_t n = row.present[reported * n_molecule_classes + m];
				if (n == 0) { continue; }
				std::atomic_ref<uint64_t>(
				    _present_cells.ones[reported][_present_cells.slot(row.species_class, m)])
				    .fetch_add(n, std::memory_order_relaxed);
			}
		}
	}

	/// The update visited every cell once, so the counts are the field's. A debug build recounts
	/// them from the field the update left; the recount walks the field, so a release build does
	/// not.
	template<FieldStorage Field> void end_pass([[maybe_unused]] const Field &Y) const {
#ifndef NDEBUG
		DEBUG_ASSERT(_present_cells == count_present_cells(Y, *_L, *_classes));
#endif
	}

	/// The present cells of the field the update left, once the pass has ended.
	[[nodiscard]] const TPresentCells &present_cells() const & noexcept { return _present_cells; }
	/// The same, moved out of a pass that is done with.
	[[nodiscard]] TPresentCells present_cells() && noexcept { return std::move(_present_cells); }
};

static_assert(DataSourcePass<TBlockPass>,
              "LOTUS's pass must answer what the block update asks a source.");

} // namespace lotus_math
