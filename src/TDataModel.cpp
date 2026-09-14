//
// Created for the split of the data sources (LOTUS / simple error model / MS data).
//
#include "TDataModel.h"
#include "Types.h"
#include "cli.h"
#include "coretools/Main/TParameters.h"
#include <cstddef>
#include <string>
#include <vector>

TDataModel::TDataModel(std::vector<std::unique_ptr<TTree>> &trees, const TDataSources &sources,
                       TMarkovField::TypeParamErrorProbability *omega, size_t n_iterations,
                       const MarkovFieldParams &markov_field_stattools_param, std::string prefix,
                       bool simulate)
    : _trees(trees), _markov_field(n_iterations, trees, omega, prefix, simulate), _omega(omega),
      _markov_field_stattools_param(markov_field_stattools_param),
#ifdef USE_LOTUS
      _lotus(trees, sources.gamma, sources.error_rate), _gamma(sources.gamma),
      _error_rate(sources.error_rate),
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
      _simple_error_model(trees, sources.epsilon_simple_model),
      _epsilon_simple_model(sources.epsilon_simple_model),
#endif
      _prefix(std::move(prefix)), _simulate(simulate) {

	// Tell stattools which parameters hang off this box. Only the compiled-in sources contribute.
	std::vector<stattools::TNodeBase *> params;
#ifdef USE_LOTUS
	params.push_back(_gamma);
	params.push_back(_error_rate);
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
	params.push_back(_epsilon_simple_model);
#endif
	params.push_back(_omega);
	for (const auto &it : _markov_field_stattools_param) { params.push_back(it.get()); }
	this->addPriorParameter(params);
}

std::string TDataModel::name() const { return "data_model"; }

void TDataModel::initialize() {
#define ACOL_INIT_DATA_SOURCE(member, Type) _##member.initialize(this, _simulate);
	ACOL_DATA_SOURCES(ACOL_INIT_DATA_SOURCE)
#undef ACOL_INIT_DATA_SOURCE

	_omega->initStorage(this, {1});

	for (auto &it : _markov_field_stattools_param) { it->initStorage(this, {0}); }

	if (!_simulate) {
		std::vector<std::string> tree_names;
		std::vector<size_t> leaf_counts;
		tree_names.reserve(_trees.size());
		leaf_counts.reserve(_trees.size());
		for (const auto &tree : _trees) {
			tree_names.push_back(tree->get_tree_name());
			leaf_counts.push_back(tree->get_number_of_leaves());
		}
		const auto n_iter   = coretools::instances::parameters().get<size_t>("iterations", 100000);
		const auto n_burnin = coretools::instances::parameters().get<size_t>("numBurnin", 10);
		const auto n_burnin_iter = coretools::instances::parameters().get<size_t>("burnin", 1000);
		_notifier.notify_start(tree_names, leaf_counts, n_iter, n_burnin, n_burnin_iter);
	}
}

void TDataModel::guessInitialValues() {
	// Note: stattools only calls this when inferring, never when simulating.
	const auto &Y = _markov_field.get_Y_matrix();
#define ACOL_GUESS_DATA_SOURCE(member, Type) _##member.guess_initial_values(Y);
	ACOL_DATA_SOURCES(ACOL_GUESS_DATA_SOURCE)
#undef ACOL_GUESS_DATA_SOURCE

	_omega->set(TypeErrorProbability(ProgramOptions::ERROR_PROBABILITY));
	coretools::instances::logfile().list(
	    "Starting the error probability omega at ", ProgramOptions::ERROR_PROBABILITY,
	    ", under an exponential prior of rate ", ProgramOptions::ERROR_PROBABILITY_PRIOR_RATE,
	    " truncated to (0, 0.5).");
}

double TDataModel::getSumLogPriorDensity(const Storage &) const { return data_log_likelihood(); }

double TDataModel::data_log_likelihood() const {
	double sum = 0.0;
#define ACOL_SUM_DATA_SOURCE(member, Type) sum += _##member.log_likelihood();
	ACOL_DATA_SOURCES(ACOL_SUM_DATA_SOURCE)
#undef ACOL_SUM_DATA_SOURCE
	return sum;
}

void TDataModel::update_markov_field() { _markov_field.update(*this, _mrf_update_iteration++); }

#ifdef USE_LOTUS
double TDataModel::calculateLLRatio(TLotus::TypeParamGamma *, size_t /*Index*/) {
	return _lotus.ll_ratio_after_parameter_move(_markov_field.get_Y_matrix());
}

double TDataModel::calculateLLRatio(TLotus::TypeParamErrorRate *, size_t /*Index*/) {
	return _lotus.ll_ratio_after_parameter_move(_markov_field.get_Y_matrix());
}

void TDataModel::updateTempVals(TLotus::TypeParamGamma *, size_t /*Index*/, bool Accepted) {
	if (!Accepted) { _lotus.revert_parameter_move(); }
}

void TDataModel::updateTempVals(TLotus::TypeParamErrorRate *, size_t /*Index*/, bool Accepted) {
	if (!Accepted) { _lotus.revert_parameter_move(); }
}
#endif

#ifdef USE_SIMPLE_ERROR_MODEL
double TDataModel::calculateLLRatio(TSimpleErrorModel::TypeParamEpsilon *, size_t /*Index*/) {
	// O(1): the likelihood depends on the data only through the disagreement count, which does not
	// depend on epsilon. stattools has already proposed, so value() is the proposal.
	return _simple_error_model.log_likelihood_ratio((double)_epsilon_simple_model->oldValue(),
	                                                (double)_epsilon_simple_model->value());
}

void TDataModel::updateTempVals(TSimpleErrorModel::TypeParamEpsilon *, size_t /*Index*/,
                                bool /*Accepted*/) {
	// Nothing to undo: the likelihood is recomputed in O(1) from the disagreement count and the
	// current epsilon, and stattools restores the value itself when a proposal is rejected. The
	// overload still has to exist -- stattools throws at runtime if it is missing.
}
#endif

double TDataModel::calculateLLRatio(TMarkovField::TypeParamErrorProbability *, size_t /*Index*/) {
	// O(1): the link's likelihood depends on the whole field only through the six counters, and
	// they do not move with the error probability. stattools has already proposed, so the
	// parameter holds the proposal and remembers the value it replaces.
	return _markov_field.link_log_likelihood_ratio();
}

void TDataModel::updateTempVals(TMarkovField::TypeParamErrorProbability *, size_t /*Index*/,
                                bool /*Accepted*/) {
	// Nothing to undo: the link's likelihood is rebuilt in O(1) from the counters and the current
	// error probability, and stattools restores the value itself when a proposal is rejected. The
	// overload still has to exist -- stattools throws at runtime if it is missing.
}

void TDataModel::_simulateUnderPrior(Storage *) {
#ifdef USE_LOTUS
	// Sizes L and sets gamma / epsilon to the values L should be simulated under. Must happen
	// before the field is simulated, because it sizes L.
	_lotus.prepare_for_simulation();
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
	_epsilon_simple_model->set(TypeEpsilonSimpleModel(ProgramOptions::EPSILON_SIMPLE_MODEL));
#endif
	// The field is simulated at the error probability the run was given, not at a draw from the
	// prior, so a simulated data set says which value produced it.
	_omega->set(TypeErrorProbability(ProgramOptions::ERROR_PROBABILITY));

	// first simulate Markov random field
	_markov_field.simulate(*this);

	// then derive each data source from that one simulated Y, and write it out
	const auto &Y = _markov_field.get_Y_matrix();
#define ACOL_SIMULATE_DATA_SOURCE(member, Type)                                                    \
	_##member.simulate_from_Y(Y);                                                                     \
	_##member.write_simulated(_prefix);
	ACOL_DATA_SOURCES(ACOL_SIMULATE_DATA_SOURCE)
#undef ACOL_SIMULATE_DATA_SOURCE
}

TNotifierStats TDataModel::_collect_notifier_stats() const {
	TNotifierStats stats;
#define ACOL_STATS_DATA_SOURCE(member, Type) _##member.contribute_stats(stats);
	ACOL_DATA_SOURCES(ACOL_STATS_DATA_SOURCE)
#undef ACOL_STATS_DATA_SOURCE
	stats.scalar_stats.push_back({"omega", {_omega->mean(0), _omega->var(0), _omega->sd(0)}});
	return stats;
}

void TDataModel::burninHasFinished() {
	_markov_field.burninHasFinished();
	const auto stats = _collect_notifier_stats();
	_notifier.notify_burnin_finished(stats.dim_names, stats.gamma_stats, stats.scalar_stats);
}

void TDataModel::oneBurninHasFinished() {
	_markov_field.oneBurninHasFinished();
	const size_t round      = ++_burnin_round;
	const auto total_rounds = coretools::instances::parameters().get<size_t>("numBurnin", 10);
	const auto stats        = _collect_notifier_stats();
	_notifier.notify_burnin_round(round, total_rounds, stats.dim_names, stats.gamma_stats,
	                              stats.scalar_stats);
}

void TDataModel::MCMCHasFinished() {
	_markov_field.MCMCHasFinished();
	const auto stats = _collect_notifier_stats();
	_notifier.notify_mcmc_finished(stats.dim_names, stats.gamma_stats, stats.scalar_stats);
}
