//
// What the field update asks about a cell: everything the data says about it.
//
// The pass in field/field_update.h reads three states from its storages -- the field cell and the
// two tree field cells at that leaf pair -- and takes the link from field/TFieldMath.h. It asks
// this for the rest. Each compiled-in data source says what it makes of the field cell, under both
// of its states, and is told afterwards which one the cell was given.
//
// This is the whole of what the pass must not know about. So the pass keeps no data source, and a
// test can run its loop against a model of its own.
//
// The two trees are not here. Each draws its own tree field as the leaf block of its node state
// (ADR-0005), so what a tree's process says about a leaf is that tree's business and reaches the
// walk through its own seam (field/tree_field_link.h).
//
// The per-source likelihood bookkeeping lives here too. It is the data's business, and the pass
// only hands each cell back once.
//

#pragma once

#include "TDataModel.h"
#include "constants.h"
#include "coretools/Math/TSumLog.h"
#include "coretools/Types/probability.h"
#include "field/field_update.h"
#include "omp.h"
#include <array>
#include <cstddef>
#include <vector>

//-----------------------------------
// Data-source bookkeeping
//-----------------------------------

/// Per-cell outcome of one field update. Each data source keeps its own likelihood bookkeeping
/// (they are independent terms), so the results are handed back separately instead of merged.
struct TCellOutcome {
#ifdef USE_LOTUS
	/// P(L_cell | Y = the drawn state). Neutral value 1.0 (log 0).
	double prob_lotus_new_state = 1.0;
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
	/// Whether the observed D cell contradicts the state the field was just given.
	bool simple_model_disagrees = false;
#endif
};

/// Per-thread accumulators for one full field update, committed to the data sources at the end.
///
/// The accumulators are bundled into a single object on purpose: `#ifdef` cannot appear inside a
/// `#pragma omp` line, and a `default(none) shared(...)` clause has to name every variable it
/// touches. One object keeps that clause identical in every build configuration.
class TDataUpdateAccumulator {
private:
#ifdef USE_LOTUS
	std::vector<coretools::TSumLogProbability> _lotus_LL;
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
	std::vector<size_t> _n_disagree;
#endif

public:
	/// Sizing happens in the body rather than in a member-initializer list, so that adding or
	/// removing a source does not require rebalancing the commas of a #ifdef'd init list.
	explicit TDataUpdateAccumulator([[maybe_unused]] size_t n_threads) {
#ifdef USE_LOTUS
		_lotus_LL.resize(n_threads);
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
		_n_disagree.assign(n_threads, 0);
#endif
	}

	/// Hot path: called once per updated cell, from inside the parallel region. Only ever touches
	/// the slot of the calling thread.
	void add([[maybe_unused]] size_t thread, [[maybe_unused]] const TCellOutcome &outcome) {
#ifdef USE_LOTUS
		_lotus_LL[thread].add(outcome.prob_lotus_new_state);
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
		_n_disagree[thread] += static_cast<size_t>(outcome.simple_model_disagrees);
#endif
	}

	/// Sums the per-thread slots and installs the results in the data sources. Called once, after
	/// the parallel region.
	void commit([[maybe_unused]] TDataModel &data_model) {
#ifdef USE_LOTUS
		double sum_new_LL = 0.0;
		for (auto &i : _lotus_LL) { sum_new_LL += i.getSum(); }
		data_model.get_lotus().update_cur_LL(sum_new_LL);
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
		size_t total_disagree = 0;
		for (const auto &i : _n_disagree) { total_disagree += i; }
		// The update visits every cell exactly once, so this is the complete disagreement count.
		data_model.get_simple_error_model().set_n_disagree(total_disagree);
#endif
	}
};

//-----------------------------------
// The model
//-----------------------------------

/// Every compiled-in data source, as the field update wants them.
///
/// An inferred chain is the only one that has a model to ask. A simulated one draws its whole
/// configuration forward and runs no update (TMarkovField::simulate).
///
/// Every thread of the update asks this one object. `factors` therefore reads and writes nothing
/// another thread also touches, and `record` writes the accumulator slot of the calling thread
/// alone.
class TFieldModel {
private:
	TDataModel &_data_model;
	TDataUpdateAccumulator &_accumulator;

public:
	TFieldModel(TDataModel &data_model, TDataUpdateAccumulator &accumulator)
	    : _data_model(data_model), _accumulator(accumulator) {}

	[[nodiscard]] field_update::TCellFactors factors([[maybe_unused]] const IndexArray &cell) const {
		// 1.0 is the neutral value, adding log(1) = 0. A build that left a source out keeps it.
		field_update::TCellFactors factors;
#ifdef USE_LOTUS
		const TLotus &lotus = _data_model.get_lotus();
		std::array<double, 2> prob_lotus{1.0, 1.0};
		lotus.calculate_LL_update_Y(cell, lotus.holds_a_record(cell), prob_lotus);
		factors.lotus = {coretools::P(prob_lotus[0]), coretools::P(prob_lotus[1])};
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
		const TSimpleErrorModel &simple = _data_model.get_simple_error_model();
		std::array<double, 2> prob_simple{};
		simple.probabilities_for_Y_update(simple.observed_state_of(cell), prob_simple);
		factors.simple_error = {coretools::P(prob_simple[0]), coretools::P(prob_simple[1])};
#endif
		return factors;
	}

	/// Each source scored both field states before the draw, and only now knows which one to keep.
	/// This writes the accumulator, so it is not const.
	void record([[maybe_unused]] const IndexArray &cell,
	            [[maybe_unused]] const field_update::TCellFactors &factors,
	            [[maybe_unused]] bool drawn) {
		TCellOutcome outcome;
#ifdef USE_LOTUS
		outcome.prob_lotus_new_state = factors.lotus[static_cast<size_t>(drawn)].get();
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
		const TSimpleErrorModel &simple = _data_model.get_simple_error_model();
		outcome.simple_model_disagrees  = simple.observed_state_of(cell) != drawn;
#endif
		_accumulator.add(static_cast<size_t>(omp_get_thread_num()), outcome);
	}
};

static_assert(field_update::FieldModel<TFieldModel>,
              "The field model must answer what the field update asks a model.");
