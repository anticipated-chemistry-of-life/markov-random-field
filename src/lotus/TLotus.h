//
// The LOTUS data source.
//
// L is a sparse binary matrix of reported occurrences, indexed on every tree and therefore the
// same shape as the field. What distinguishes it from the other data sources is that the
// probability of reporting a present metabolite is not a flat error rate but a per-cell research
// effort derived from paper counts and the inferred rate gamma.
//
// This class is not a stattools box: its parameters (gamma, epsilon) hang off TDataModel, which
// owns this object and forwards the MCMC callbacks. See TDataModel.h.
//

#pragma once

#include "Types.h"

#ifdef USE_LOTUS

#include "cli.h"
#include "constants.h"
#include "data_sources/data_source.h"
#include "lotus/TLotusMath.h"
#include "ntfy/TNtfyNotifier.h"
#include "omp.h"
#include "stattools/ParametersObservations/TParameter.h"
#include "storages/TSparse.h"
#include "storages/storage_backend.h"
#include "tree/TTree.h"
#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <vector>

class TLotus {
public:
	using TypeParamGamma     = stattools::TParameter<SpecGamma, TDataModel>;
	using TypeParamErrorRate = stattools::TParameter<SpecErrorRate, TDataModel>;

private:
	// trees should be a const ref because we don't want to change the trees and don't want to copy
	// them
	const std::vector<std::unique_ptr<TTree>> &_trees;

	// data
	TSparseBinary _L;

	/// Raw publication counts per (tree, leaf). Constant data; the log transform and the
	/// detection rates are applied by the reporting model.
	std::vector<std::vector<size_t>> _paper_counts;

	/// The per-cell emission, memoized against the current gamma and error rate. Rebuilt as a whole
	/// value whenever either parameter moves, so there is nothing to refresh in place and nothing
	/// to revert.
	std::optional<lotus_math::TReportingModel> _reporting_model;

	// parameters gamma
	TypeParamGamma *_gamma = nullptr;

	// Error rate of lotus
	TypeParamErrorRate *_error_rate = nullptr;

	// temporary values
	double _oldLL = 0.0;
	double _curLL = 0.0;

	// private functions
	/// Gather the raw paper counts of every tree. Both the inference path and the simulation path
	/// need this.
	void _gather_paper_counts();
	/// Build a reporting model from the current gamma and error rate.
	[[nodiscard]] lotus_math::TReportingModel _build_reporting_model() const;
	[[nodiscard]] const lotus_math::TReportingModel &_reporting() const {
		return _reporting_model.value();
	}

public:
	TLotus(const std::vector<std::unique_ptr<TTree>> &trees, TypeParamGamma *gamma,
	       TypeParamErrorRate *error_rate);
	~TLotus() = default;

	/// Reads the LOTUS file (inference only) and sizes the parameter storages. `box` is the
	/// TDataModel the parameters hang off; it has to be passed in because this class is no longer
	/// the box itself.
	void initialize(TDataModel *box, bool simulate);
	void load_from_file(const std::string &filename);
	void guess_initial_values(const TFieldStorage &Y);

	[[nodiscard]] double calculate_log_likelihood_of_L(const TFieldStorage &Y) const;
	/// Satisfies DataSource.
	[[nodiscard]] double log_likelihood() const { return _curLL; }

	// --- hooks used by the field update (see TMarkovField::_update_Y) ---

	/// The cursor a row of the block update walks `_L`'s sorted-ones cache with. One per field
	/// row, seeded at that row's first cell, and owned by the row rather than by the thread that
	/// happens to draw it.
	using RecordCursor = TSparseBinary::OnesCursor;

	/// Makes `_L`'s sorted-ones cache fresh, so that every row of the block update may seed a
	/// cursor of its own from inside the parallel region. Has to run once, single-threaded,
	/// before that region starts, because the sort is what is not safe to do concurrently.
	///
	/// It is told nothing about how the traversal shares its rows out. The cursor used to be one
	/// per thread, which meant this had to be handed the exact chunk size
	/// `schedule(static, ...)` would use -- and so the block update had to derive its schedule
	/// from `ProgramOptions::NUMBER_OF_THREADS` rather than from the team it actually got.
	void prepare_for_block_update() const { _L.refresh_ones(); }

	/// A cursor seeded at the first cell of one field row. `prepare_for_block_update` has to have
	/// run first. Safe to call from every row at once: it reads the cache and never rebuilds it.
	[[nodiscard]] RecordCursor cursor_for_row(size_t first_cell_of_row) const {
		return _L.fresh_ones_cursor_from(first_cell_of_row);
	}

	/// Whether LOTUS holds a record for one cell of the field. L has the field's dimensions, so
	/// the field's index is already L's. The update asks this once per leaf pair, in ascending
	/// linear-index order down the row -- a row walks its columns in order -- which is what lets
	/// this walk the sorted-ones cache forward with the row's own cursor instead of hashing.
	/// Nothing else calls this out of order: it is `TBlockModel::factors`'s alone.
	[[nodiscard]] bool holds_a_record(RecordCursor &cursor,
	                                  const IndexArray &index_in_leaves_space) const {
		const size_t linear_index = _L.get_linear_index_in_container_space(index_in_leaves_space);
		return cursor.advance_to_and_check(linear_index);
	}

	/// prob[0] = P(L_cell | Y = 0), prob[1] = P(L_cell | Y = 1). `reports_the_cell` is whether
	/// LOTUS holds a record for it, which the caller reads through the accessor above.
	void calculate_LL_update_Y(const IndexArray &index_in_leaves_space, bool reports_the_cell,
	                           std::array<double, 2> &prob) const;
	/// The field update accumulates the new likelihood as it goes and installs it here at the end.
	void update_cur_LL(double cur_LL) { _curLL = cur_LL; }

	// --- MCMC moves on gamma / epsilon ---
	// Both recompute the full LOTUS likelihood: gamma and epsilon enter every cell. Each builds a
	// candidate reporting model and installs it; on rejection the candidate is simply dropped, so
	// the only thing to restore is the cached likelihood.

	[[nodiscard]] double ll_ratio_after_parameter_move(const TFieldStorage &Y);
	void revert_parameter_move();

	// --- simulation ---

	/// Sizes L and sets gamma / epsilon to the values the data should be simulated under. The
	/// parameter storages are sized in initialize(), on both paths. Must run before
	/// simulate_from_Y.
	void prepare_for_simulation();
	/// Draws every cell of L given the simulated Y. Satisfies DataSource.
	void simulate_from_Y(const TFieldStorage &Y);
	/// Writes the simulated L as <prefix>_simulated_lotus.tsv: a header naming every tree, then
	/// one row of leaf node ids per cell whose state is 1. Satisfies DataSource.
	void write_simulated(const std::string &prefix) const;

	// --- accessors ---

	[[nodiscard]] const TSparseBinary &get_L() const { return _L; }
	[[nodiscard]] std::vector<std::string> tree_names() const;
	[[nodiscard]] std::vector<TNtfyNotifier::ParamStats> gamma_stats() const;
	[[nodiscard]] TNtfyNotifier::ParamStats error_rate_stats() const;

	/// Satisfies DataSource.
	void contribute_stats(TNotifierStats &stats) const;

	static std::string get_filename_lotus() { return ProgramOptions::LOTUS_FILENAME; }
};

static_assert(DataSource<TLotus>,
              "TLotus must satisfy the DataSource concept (data_sources/data_source.h)");

#endif // USE_LOTUS
