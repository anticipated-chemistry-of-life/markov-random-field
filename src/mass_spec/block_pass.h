//
// Mass spectrometry's pass through one block update (data_sources/data_source_pass.h).
//
// The pass says {P(runs | Y = 0), P(runs | Y = 1)} about a cell, given the current assignments,
// filter probabilities and contamination (`mass_spec_math::per_cell_factor`). It tallies nothing:
// the mass-spec likelihood reads the field itself (TMassSpecMath.h), so the row and end hooks are
// empty. A test builds the pass from the runs and the parameter values alone, without a tree or
// stattools. `TMSMSData` hands a pass out and takes it back (`begin_block_update`,
// `end_block_update`).
//
// Not guarded by USE_MS_DATA, like TMassSpecMath.h, so it compiles and is unit tested in every
// build configuration.
//

#pragma once

#include "coretools/Containers/TNestedVector.h"
#include "coretools/Containers/TView.h"
#include "coretools/Main/TError.h"
#include "data_sources/data_source_pass.h"
#include "mass_spec/TMassSpecMath.h"
#include "mass_spec/msms_run.h"
#include "storages/storage_concepts.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace mass_spec_math {

class TBlockPass {
private:
	/// The runs of every species leaf, at their current assignments. The pass only reads them.
	const coretools::TNestedVector<TMassSpecRun> *_runs = nullptr;
	/// The probability that a present molecule passes a filter, in linear space,
	/// `[filter * _n_molecules + molecule]`. Read once, when the pass is built: the parameter does
	/// not move during the update, so these are the doubles a read per cell would give.
	std::vector<double> _filter_probabilities;
	size_t _n_molecules   = 0;
	/// The contamination probability, read once for the same reason.
	double _contamination = 0.0;

public:
	/// What one row carries from cell to cell: the runs of the row's species leaf. Empty for a
	/// species without runs, whose every cell factor is then `{1, 1}`.
	struct TRow {
		coretools::TConstView<TMassSpecRun> runs;
	};

	/// `filter_probabilities` is laid out `[filter * n_molecules + molecule]`, over every filter
	/// a run of `runs` names.
	TBlockPass(const coretools::TNestedVector<TMassSpecRun> &runs,
	           std::vector<double> filter_probabilities, size_t n_molecules, double contamination)
	    : _runs(&runs), _filter_probabilities(std::move(filter_probabilities)),
	      _n_molecules(n_molecules), _contamination(contamination) {
		DEBUG_ASSERT(_n_molecules > 0 && _filter_probabilities.size() % _n_molecules == 0);
	}

	/// Nothing to ready: neither the runs nor the parameters are written during the update.
	void begin_pass() {}

	[[nodiscard]] TRow begin_row(size_t species_leaf, size_t /*first_field_cell*/) const {
		return TRow{.runs = _runs->at(species_leaf)};
	}

	/// {P(runs | Y = 0), P(runs | Y = 1)}: `{1, 1}` for a species without runs.
	[[nodiscard]] std::array<double, 2> cell_factor(TRow &row, size_t molecule_leaf) const {
		DEBUG_ASSERT(molecule_leaf < _n_molecules);
		return per_cell_factor(
		    row.runs, static_cast<uint32_t>(molecule_leaf),
		    [this](size_t filter, uint32_t molecule) {
			    return _filter_probabilities[filter * _n_molecules + molecule];
		    },
		    _contamination);
	}

	/// Nothing to tally.
	void record(TRow & /*row*/, size_t /*molecule_leaf*/, bool /*drawn_y*/) const {}
	void end_row(const TRow & /*row*/) {}
	template<FieldStorage Field> void end_pass(const Field & /*Y*/) {}

	/// The filter probabilities the pass read, moved out of a pass that is done with, so that the
	/// data source can refill the same buffer for the next update.
	[[nodiscard]] std::vector<double> filter_probabilities() && noexcept {
		return std::move(_filter_probabilities);
	}
};

static_assert(DataSourcePass<TBlockPass>,
              "Mass spectrometry's pass must answer what the block update asks a source.");

} // namespace mass_spec_math
