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
#include "lotus/block_pass.h"
#include "ntfy/TNtfyNotifier.h"
#include "omp.h"
#include "stattools/ParametersObservations/TParameter.h"
#include "storages/TSparse.h"
#include "storages/storage_backend.h"
#include "tree/TTree.h"
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
	/// The same counts, as the class of each leaf: what the present-cell counts are indexed by.
	lotus_math::TPaperCountClasses _classes;

	/// The per-cell emission, memoized against the current gamma and error rate. Rebuilt as a whole
	/// value whenever either parameter moves, so there is nothing to refresh in place and nothing
	/// to revert.
	std::optional<lotus_math::TReportingModel> _reporting_model;

	/// The field, as much of it as the likelihood of L reads: its present cells by kind. Counted
	/// from the field when the chain starts, and by every block update after (ADR-0015).
	lotus_math::TPresentCells _present_cells;

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
	[[nodiscard]] std::vector<double> _gammas() const;
	/// The likelihood of L at the current present cells, gamma and error rate.
	[[nodiscard]] double _log_likelihood_of_present_cells() const;
	/// Replaces the present cells, and the likelihood read from them.
	void _install_present_cells(lotus_math::TPresentCells cells);

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

	/// Satisfies DataSource.
	[[nodiscard]] double log_likelihood() const { return _curLL; }

	// --- the block update ---

	/// This source's part in one block update (lotus/block_pass.h).
	using TBlockPass = lotus_math::TBlockPass;

	/// A pass over L under the current reporting model. L, the paper counts, gamma and the error
	/// rate do not move until the pass is handed back.
	[[nodiscard]] TBlockPass begin_block_update() const {
		return TBlockPass(_L, _classes, _reporting());
	}

	/// Installs the present cells the pass counted over the whole field, and the likelihood read
	/// from them.
	void end_block_update(TBlockPass &&pass);

	// --- MCMC moves on gamma / epsilon ---
	// Both enter every present cell, and neither moves the field: each builds a candidate reporting
	// model and reads the likelihood off the present cells the field left. On rejection the
	// candidate is simply dropped, so the only thing to restore is the cached likelihood.

	[[nodiscard]] double ll_ratio_after_parameter_move();
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
