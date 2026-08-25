//
// The LOTUS data source.
//
// L is a sparse binary matrix of reported occurrences. Unlike the other data sources it may live in
// fewer dimensions than Y (the remaining ones are collapsed away by TCollapser), and the
// probability of reporting a present metabolite is not a flat error rate but a per-cell research
// effort derived from paper counts and the inferred rate gamma.
//
// This class is not a stattools box: its parameters (gamma, epsilon) hang off TDataModel, which
// owns this object and forwards the MCMC callbacks. See TDataModel.h.
//

#pragma once

#include "Types.h"

#ifdef USE_LOTUS

#include "TCollapser.h"
#include "TCurrentState.h"
#include "cli.h"
#include "constants.h"
#include "lotus/TLotusMath.h"
#include "ntfy/TNtfyNotifier.h"
#include "stattools/ParametersObservations/TParameter.h"
#include "storages/y_storage/TStorageYMatrix.h"
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
	TStorageYMatrix _L;

	/// Raw publication counts per (kept dimension, leaf). Constant data; the log transform and the
	/// detection rates are applied by the reporting model.
	std::vector<std::vector<size_t>> _paper_counts;

	/// The per-cell emission, memoized against the current gamma and error rate. Rebuilt as a whole
	/// value whenever either parameter moves, so there is nothing to refresh in place and nothing
	/// to revert.
	std::optional<lotus_math::TReportingModel> _reporting_model;

	// how to collapse
	TCollapser _collapser;

	// parameters gamma
	TypeParamGamma *_gamma = nullptr;

	// Error rate of lotus
	TypeParamErrorRate *_error_rate = nullptr;

	// temporary values
	double _oldLL = 0.0;
	double _curLL = 0.0;
	TCurrentState _tmp_state_along_last_dim;

	// private functions
	/// Gather the raw paper counts of every kept dimension. Both the inference path and the
	/// simulation path need this, and neither can do it before the collapser is initialized.
	void _gather_paper_counts();
	/// Build a reporting model from the current gamma and error rate.
	[[nodiscard]] lotus_math::TReportingModel _build_reporting_model() const;
	[[nodiscard]] const lotus_math::TReportingModel &_reporting() const {
		return _reporting_model.value();
	}
	[[nodiscard]] double
	_calculate_log_likelihood_of_L_no_collapsing(const TStorageYMatrix &Y) const;

public:
	TLotus(const std::vector<std::unique_ptr<TTree>> &trees, TypeParamGamma *gamma,
	       TypeParamErrorRate *error_rate);
	~TLotus() = default;

	/// Reads the LOTUS file (inference only) and sizes the parameter storages. `box` is the
	/// TDataModel the parameters hang off; it has to be passed in because this class is no longer
	/// the box itself.
	void initialize(TDataModel *box, bool simulate);
	void load_from_file(const std::string &filename);
	void guess_initial_values(const TStorageYMatrix &Y);

	[[nodiscard]] double calculate_log_likelihood_of_L(const TStorageYMatrix &Y) const;
	[[nodiscard]] double cur_LL() const { return _curLL; }

	// --- hooks used by the Y sweep (see TMarkovField::_update_Y) ---

	void fill_tmp_state_along_last_dim(const IndexArray &start_index_clique_along_last_dim,
	                                   size_t K);
	void calculate_LL_update_Y(const IndexArray &index_in_leaves_space, size_t index_for_tmp_state,
	                           bool old_state, std::array<double, 2> &prob) const;
	/// The Y sweep accumulates the new likelihood as it goes and installs it here at the end.
	void update_cur_LL(double cur_LL) { _curLL = cur_LL; }

	// --- MCMC moves on gamma / epsilon ---
	// Both recompute the full LOTUS likelihood: gamma and epsilon enter every cell. Each builds a
	// candidate reporting model and installs it; on rejection the candidate is simply dropped, so
	// the only thing to restore is the cached likelihood.

	[[nodiscard]] double ll_ratio_after_parameter_move(const TStorageYMatrix &Y);
	void revert_parameter_move();

	// --- simulation ---

	/// Sizes L and the parameter storages for simulation, and sets gamma / epsilon to the values
	/// the data should be simulated under. Must run before simulate_L_from_Y.
	void prepare_for_simulation(TDataModel *box);
	/// Draws every cell of L given the simulated Y.
	void simulate_L_from_Y(const TStorageYMatrix &Y);
	/// Writes the simulated L as <prefix>_simulated_lotus.tsv: a header naming the kept trees, then
	/// one row of leaf node ids per cell whose state is 1.
	void write_simulated_L(const std::string &prefix) const;

	// --- accessors ---

	[[nodiscard]] const TStorageYMatrix &get_L() const { return _L; }
	[[nodiscard]] std::vector<std::string> kept_tree_names() const;
	[[nodiscard]] std::vector<TNtfyNotifier::ParamStats> gamma_stats() const;
	[[nodiscard]] TNtfyNotifier::ParamStats error_rate_stats() const;

	static std::string get_filename_lotus() { return ProgramOptions::LOTUS_FILENAME; }
};

#endif // USE_LOTUS
