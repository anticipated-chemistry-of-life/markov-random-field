//
// The likelihood box that anchors the model in the stattools DAG.
//
// TDataModel owns the Markov random field (the latent Y and the per-tree Z) and every data source
// that is compiled in. Each per-cell source (LOTUS, the simple error model) is an independent
// likelihood term over the same Y, and satisfies the DataSource concept in
// data_sources/data_source.h; ACOL_DATA_SOURCES(X) there lists whichever are compiled in and
// drives six of this class's per-source sites (initialize, guess_initial_values, the likelihood
// sum, notifier-stats collection, simulate + write, the accessor pair) from that one list. Mass
// spectrometry (USE_MS_DATA, TMSMSData) is not this shape -- it scores an assignment problem, not
// a per-cell state -- and stays dormant, hand-wired outside TDataModel.
//
// Three things per source stay hand-written rather than generated: the member declaration, the
// constructor's forwarding of that source's parameters, and the params-vector push in the
// constructor body below. A source's parameter count varies source to source (LOTUS has two,
// the simple error model one), so those three sites cannot be expressed as one uniform macro row
// without losing that. The parameter-move dispatch (calculateLLRatio / updateTempVals) is
// likewise hand-written per parameter type, because stattools dispatches a move by the
// parameter's own type, not by its owning source.
//
// At least one per-cell source must be compiled in; Types.h enforces that with a static_assert.
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
#include "data_sources/data_source.h"
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
	/// TNotifierStats itself lives in data_sources/data_source.h, so DataSource::contribute_stats
	/// can name it.
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

	/// `log p(L, D | Y)`: the same sum, under the name the joint density trace asks for it by. It
	/// is the data factor of the ADR-0005 factorisation, and stattools' own accessor takes a
	/// storage the field has nothing to hand it.
	///
	/// The LOTUS records and the simple error model, and nothing else. The mass spectrometry
	/// source is dormant -- nothing builds it, and it hangs off the field rather than this class --
	/// so it has no term to add. A build that wakes it adds it here, or the joint density stops
	/// being the whole of the model.
	[[nodiscard]] double data_log_likelihood() const;

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
	// Accessor pair generated for every source in ACOL_DATA_SOURCES.
#define ACOL_DATA_SOURCE_ACCESSOR(member, Type)                                                    \
	[[nodiscard]] const Type &get_##member() const { return _##member; }                             \
	[[nodiscard]] Type &get_##member() { return _##member; }
	ACOL_DATA_SOURCES(ACOL_DATA_SOURCE_ACCESSOR)
#undef ACOL_DATA_SOURCE_ACCESSOR
};
