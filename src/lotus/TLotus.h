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

#include "TCurrentState.h"
#include "cli.h"
#include "constants.h"
#include "lotus/TLotusMath.h"
#include "ntfy/TNtfyNotifier.h"
#include "stattools/ParametersObservations/TParameter.h"
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
	TFieldStorage _L;

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
	TCurrentState _tmp_state_along_last_dim;

	// private functions
	/// The number of leaves of each tree, which is the shape of L: records are indexed on every
	/// tree, so L has the same dimensions as the field.
	[[nodiscard]] std::vector<size_t> _leaf_counts() const;
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
	[[nodiscard]] double cur_LL() const { return _curLL; }

	// --- hooks used by the field update (see TMarkovField::_update_Y) ---

	void fill_tmp_state_along_last_dim(const IndexArray &start_index_clique_along_last_dim,
	                                   size_t K);
	void calculate_LL_update_Y(const IndexArray &index_in_leaves_space, size_t index_for_tmp_state,
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
	/// simulate_L_from_Y.
	void prepare_for_simulation();
	/// Draws every cell of L given the simulated Y.
	void simulate_L_from_Y(const TFieldStorage &Y);
	/// Writes the simulated L as <prefix>_simulated_lotus.tsv: a header naming every tree, then
	/// one row of leaf node ids per cell whose state is 1.
	void write_simulated_L(const std::string &prefix) const;

	// --- accessors ---

	[[nodiscard]] const TFieldStorage &get_L() const { return _L; }
	[[nodiscard]] std::vector<std::string> tree_names() const;
	[[nodiscard]] std::vector<TNtfyNotifier::ParamStats> gamma_stats() const;
	[[nodiscard]] TNtfyNotifier::ParamStats error_rate_stats() const;

	static std::string get_filename_lotus() { return ProgramOptions::LOTUS_FILENAME; }
};

#endif // USE_LOTUS
