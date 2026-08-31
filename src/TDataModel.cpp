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
                       size_t n_iterations, const MarkovFieldParams &markov_field_stattools_param,
                       std::string prefix, bool simulate)
    : _trees(trees), _markov_field(n_iterations, trees, prefix),
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
	for (const auto &it : _markov_field_stattools_param) { params.push_back(it.get()); }
	this->addPriorParameter(params);
}

std::string TDataModel::name() const { return "data_model"; }

void TDataModel::initialize() {
#ifdef USE_LOTUS
	_lotus.initialize(this, _simulate);
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
	// D always spans the full leaf space, so it can be sized from the trees alone -- in both
	// inference (where load_from_file then fills it) and simulation (where simulate_D_from_Y does).
	_simple_error_model.initialize_storage();
	_epsilon_simple_model->initStorage(this, {1});
	if (!_simulate) { _simple_error_model.load_from_file(ProgramOptions::SIMPLE_DATA_FILENAME); }
#endif

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
#ifdef USE_LOTUS
	_lotus.guess_initial_values(Y);
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
	_epsilon_simple_model->set(TypeEpsilonSimpleModel(ProgramOptions::EPSILON_SIMPLE_MODEL));
	_simple_error_model.guess_initial_values(Y);
	coretools::instances::logfile().list(
	    "Using simple error model data with initial epsilon_simple_model = ",
	    ProgramOptions::EPSILON_SIMPLE_MODEL, ".");
#endif
}

double TDataModel::getSumLogPriorDensity(const Storage &) const {
	double sum = 0.0;
#ifdef USE_LOTUS
	sum += _lotus.cur_LL();
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
	sum += _simple_error_model.log_likelihood();
#endif
	return sum;
}

void TDataModel::fill_tmp_state_along_last_dim(const IndexArray &start_index_clique_along_last_dim,
                                               size_t K) {
#ifdef USE_LOTUS
	_lotus.fill_tmp_state_along_last_dim(start_index_clique_along_last_dim, K);
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
	_simple_error_model.fill_tmp_state_along_last_dim(start_index_clique_along_last_dim, K);
#endif
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

void TDataModel::_simulateUnderPrior(Storage *) {
#ifdef USE_LOTUS
	// Sizes L and sets gamma / epsilon to the values L should be simulated under. Must happen
	// before the field is simulated, because it sizes L.
	_lotus.prepare_for_simulation();
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
	_epsilon_simple_model->set(TypeEpsilonSimpleModel(ProgramOptions::EPSILON_SIMPLE_MODEL));
#endif

	// first simulate Markov random field
	_markov_field.simulate(*this);

	// then derive each data source from that one simulated Y, and write it out
	const auto &Y = _markov_field.get_Y_matrix();
#ifdef USE_LOTUS
	_lotus.simulate_L_from_Y(Y);
	_lotus.write_simulated_L(_prefix);
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
	_simple_error_model.simulate_D_from_Y(Y);
	_simple_error_model.write_simulated_D(_prefix);
#endif
}

TDataModel::TNotifierStats TDataModel::_collect_notifier_stats() const {
	TNotifierStats stats;
#ifdef USE_LOTUS
	stats.dim_names   = _lotus.tree_names();
	stats.gamma_stats = _lotus.gamma_stats();
	stats.scalar_stats.push_back({"epsilon", _lotus.error_rate_stats()});
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
	stats.scalar_stats.push_back({"epsilon_simple_model",
	                              {_epsilon_simple_model->mean(0), _epsilon_simple_model->var(0),
	                               _epsilon_simple_model->sd(0)}});
#endif
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
