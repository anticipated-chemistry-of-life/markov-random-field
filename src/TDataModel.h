//
// The likelihood box that anchors the model in the stattools DAG.
//
// TDataModel owns the Markov random field (the latent Y and the per-tree Z) and every source of
// information that is compiled in. Each source is an independent likelihood term over the same Y:
//
//   USE_LOTUS              TLotus            reported occurrences, research effort
//   USE_SIMPLE_ERROR_MODEL TSimpleErrorModel a flat-error-rate noisy copy of Y
//   USE_MS_DATA            TMSMSData         mass spectrometry (still dormant)
//
// At least one must be compiled in; Types.h enforces that with a static_assert. Sources are guarded
// with #ifdef rather than `if constexpr` because a discarded `if constexpr` branch must still
// name-resolve, which it cannot once the members are gone.
//
// This class is deliberately the *only* stattools box in the model. A second box would get its own
// _simulateUnderPrior call, with no ordering guarantee relative to the single Markov field
// simulation that produces Y -- so all sources must derive their simulated data from the same Y
// here.
//

#pragma once

#include "lotus/TLotus.h"
#include "TMarkovField.h"
#include "Types.h"
#include "ntfy/TNtfyNotifier.h"
#include "simple_error_model/TSimpleErrorModel.h"
#include "tree/TTree.h"
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

/// The parameters each compiled-in data source needs. Bundled so that TDataModel's constructor does
/// not have to change shape with the build configuration.
struct TDataSources {
#ifdef USE_LOTUS
	TLotus::TypeParamGamma *gamma          = nullptr;
	TLotus::TypeParamErrorRate *error_rate = nullptr;
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
	TSimpleErrorModel::TypeParamEpsilon *epsilon_simple_model = nullptr;
#endif
};

class TDataModel : public stattools::prior::TBaseLikelihoodPrior<stattools::TObservationBase,
                                                                 TypeDataObs, NumDimDataObs> {
public:
	// some type aliases, for better readability
	using BoxType = TDataModel;
	using Base    = stattools::prior::TBaseLikelihoodPrior<stattools::TObservationBase, TypeDataObs,
	                                                       NumDimDataObs>;
	using typename Base::Storage;
	using typename Base::UpdatedStorage;

private:
	// trees should be a const ref because we don't want to change the trees and don't want to copy
	// them
	const std::vector<std::unique_ptr<TTree>> &_trees;

	// the latent field: Y and, per tree, Z
	TMarkovField _markov_field;

	// the error probability standing between the two tree fields and the field; owned by TModel,
	// moved by stattools
	TMarkovField::TypeParamErrorProbability *_omega = nullptr;

	// Markov field parameter (only needed for stattools purposes to build a valid DAG)
	const MarkovFieldParams &_markov_field_stattools_param;

	// --- data sources ---
#ifdef USE_LOTUS
	TLotus _lotus;
	TLotus::TypeParamGamma *_gamma          = nullptr;
	TLotus::TypeParamErrorRate *_error_rate = nullptr;
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
	TSimpleErrorModel _simple_error_model;
	TSimpleErrorModel::TypeParamEpsilon *_epsilon_simple_model = nullptr;
#endif

	// output file
	std::string _prefix;

	// simulate or infer?
	bool _simulate = false;

	// notifications
	TNtfyNotifier _notifier;

	// counts how many burnin rounds have finished (oneBurninHasFinished no longer
	// receives the round number from stattools)
	size_t _burnin_round = 0;

	// monotonic counter of markov field updates (the update function registered via
	// addFuncToUpdate is no longer passed the iteration number by stattools)
	size_t _mrf_update_iteration = 0;

	void _simulateUnderPrior(Storage *) override;

	/// Everything the ntfy notifications report, collected once for all three hooks. Which entries
	/// exist depends on the compiled-in sources; a build without LOTUS reports no gamma at all.
	struct TNotifierStats {
		std::vector<std::string> dim_names;
		std::vector<TNtfyNotifier::ParamStats> gamma_stats;
		std::vector<TNtfyNotifier::NamedStats> scalar_stats;
	};
	[[nodiscard]] TNotifierStats _collect_notifier_stats() const;

public:
	/// `omega` is the field's own parameter, not a data source's, so it is passed apart from
	/// `sources`: every data source stands above the field, and the error probability stands
	/// inside it.
	TDataModel(std::vector<std::unique_ptr<TTree>> &trees, const TDataSources &sources,
	           TMarkovField::TypeParamErrorProbability *omega, size_t n_iterations,
	           const MarkovFieldParams &markov_field_stattools_param, std::string prefix,
	           bool simulate);
	~TDataModel() override = default;

	[[nodiscard]] std::string name() const override;
	void initialize() override;
	void guessInitialValues() override;

	void burninHasFinished() override;
	void oneBurninHasFinished() override;
	void MCMCHasFinished() override;

	/// Sum of the log-likelihoods of every compiled-in data source. They are independent terms over
	/// the same Y, so they simply add.
	[[nodiscard]] double getSumLogPriorDensity(const Storage &) const override;

	void update_markov_field();

	// --- MCMC moves, dispatched by stattools on the parameter pointer type ---
	// Each move only recomputes the likelihood of the source that parameter belongs to: the other
	// sources do not depend on it and cancel in the ratio.

#ifdef USE_LOTUS
	[[nodiscard]] double calculateLLRatio(TLotus::TypeParamGamma *, size_t /*Index*/);
	[[nodiscard]] double calculateLLRatio(TLotus::TypeParamErrorRate *, size_t /*Index*/);
	void updateTempVals(TLotus::TypeParamGamma *, size_t /*Index*/, bool Accepted);
	void updateTempVals(TLotus::TypeParamErrorRate *, size_t /*Index*/, bool Accepted);
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
	[[nodiscard]] double calculateLLRatio(TSimpleErrorModel::TypeParamEpsilon *, size_t /*Index*/);
	void updateTempVals(TSimpleErrorModel::TypeParamEpsilon *, size_t /*Index*/, bool Accepted);
#endif

	// The error probability is not a data source's parameter, so it is never behind an #ifdef: the
	// link stands between the tree fields and the field in every build.
	[[nodiscard]] double calculateLLRatio(TMarkovField::TypeParamErrorProbability *,
	                                      size_t /*Index*/);
	void updateTempVals(TMarkovField::TypeParamErrorProbability *, size_t /*Index*/, bool Accepted);

	// --- accessors ---

	[[nodiscard]] const TMarkovField &get_markov_field() const { return _markov_field; }
#ifdef USE_LOTUS
	[[nodiscard]] const TLotus &get_lotus() const { return _lotus; }
	[[nodiscard]] TLotus &get_lotus() { return _lotus; }
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
	[[nodiscard]] const TSimpleErrorModel &get_simple_error_model() const {
		return _simple_error_model;
	}
	[[nodiscard]] TSimpleErrorModel &get_simple_error_model() { return _simple_error_model; }
#endif
};
