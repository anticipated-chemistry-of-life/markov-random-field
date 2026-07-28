/*
 * TExampleTask.cpp
 *
 *  Created on: May 15, 2023
 *      Author: mvisani
 */

#include "TCore.h"
#include "Types.h"
#include "cli.h"
#include "coretools/Main/TError.h"
#include "stattools/MCMC/TMCMC.h"
#include <memory>

using namespace coretools::instances;

//--------------------------------------
// TModel
//--------------------------------------

void TModel::_create_tree(size_t dimension, const std::string &filename,
                          const std::string &tree_name) {
	const std::string &prefix = _prefix;

	// create mean log nu
	_mean_log_nu.push_back(std::make_unique<stattools::TParameter<SpecMeanLogNu, PriorOnLogNu>>(
	    tree_name + "_mean_log_nu", &_prior_on_mean_log_nu,
	    stattools::TParameterDefinition{prefix, ProgramOptions::FIXED_PRIOR_ON_MEAN_LOG_NU}));

	// create var log nu
	_var_log_nu.push_back(std::make_unique<stattools::TParameter<SpecVarLogNu, PriorOnLogNu>>(
	    tree_name + "_var_log_nu", &_prior_on_var_log_nu,
	    stattools::TParameterDefinition{prefix, ProgramOptions::FIXED_PRIOR_ON_VAR_LOG_NU}));

	// create prior on log nu
	_prior_on_log_nu.push_back(
	    std::make_unique<PriorOnLogNu>(_mean_log_nu.back().get(), _var_log_nu.back().get()));

	// create log nu
	stattools::TParameterDefinition def_log_nu(prefix);
	def_log_nu.excludeFromDAGUpdates(true); // never update
	_log_nu.push_back(std::make_unique<stattools::TParameter<SpecLogNu, TTree>>(
	    tree_name + "_log_nu", _prior_on_log_nu.back().get(), def_log_nu));

	// create alpha
	stattools::TParameterDefinition def_alpha(prefix);
	def_alpha.excludeFromDAGUpdates(true); // never update
	_alpha.push_back(std::make_unique<stattools::TParameter<SpecAlpha, TTree>>(
	    tree_name + "_alpha", &_prior_on_alpha, def_alpha));

	// create branch lengths
	stattools::TParameterDefinition def_branch_lengths(prefix + "_" + tree_name);
	def_branch_lengths.excludeFromDAGUpdates(true); // never update
	_binned_branch_lengths.push_back(
	    std::make_unique<stattools::TParameter<SpecBinnedBranches, TTree>>(
	        tree_name + "_branch_lengths", &_prior_on_binned_branch_lengths, def_branch_lengths));

	// create tree
	_trees.emplace_back(std::make_unique<TTree>(dimension, filename, tree_name, _alpha.back().get(),
	                                            _log_nu.back().get(),
	                                            _binned_branch_lengths.back().get()));

	// create markov field (only for stattools purposes such that a valid DAG can be built)
	stattools::TParameterDefinition def_markov_field;
	def_markov_field.excludeFromDAGUpdates(true); // never update
	_markov_field_stattools_param.emplace_back(std::make_unique<ParamMarkovField>(
	    tree_name + "_MRF", _trees.back().get(), def_markov_field));
}

void TModel::_create_trees() {
	using namespace coretools::instances;
	// read filenames
	std::string filename_tree_species   = parameters().get("tree_species");
	std::string filename_tree_molecules = parameters().get("tree_molecules");
	std::vector<std::string> filenames_tree_others;
	if (parameters().exists("tree_others")) {
		parameters().fill("tree_others", filenames_tree_others);
	}

	size_t num_trees = 2 + filenames_tree_others.size();
	_trees.reserve(num_trees);

	// first tree: molecules
	_create_tree(0, filename_tree_species, "species");

	// middle trees: all others (e.g. tissues)
	for (size_t i = 1; i < num_trees - 1; ++i) {
		std::string name = coretools::str::split(filenames_tree_others[i - 1], ':');
		if (name.empty()) {
			throw coretools::TUserError("Argument 'tree_others': Please provide a name for each "
			                            "other tree, separated by a : from the "
			                            "filename (e.g. myTreeName:pathToFile)");
		}
		_create_tree(i, filenames_tree_others[i - 1], name);
	}
	// last tree: species
	_create_tree(num_trees - 1, filename_tree_molecules, "molecules");

	for (auto &tree : _trees) { tree->initialize_cliques_and_Z(_trees); }
}

TModel::TModel(size_t n_iterations, const std::string &prefix, bool simulate)
    : _prefix(prefix)
#ifdef USE_LOTUS
      ,
      _gamma("gamma", &_prior_on_gamma, {_prefix, ProgramOptions::FIXED_PRIOR_ON_GAMMA}),
      _error_rate("epsilon", &_prior_on_error_rate,
                  {_prefix, ProgramOptions::FIXED_PRIOR_ON_EPSILON})
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
      ,
      _epsilon_simple_model("epsilon_simple_model", &_prior_on_epsilon_simple_model,
                            {_prefix, ProgramOptions::FIXED_PRIOR_ON_EPSILON_SIMPLE_MODEL})
#endif
#ifdef USE_MS_DATA
      ,
      _mass_spec_filters("filter_proba", &_prior_on_mass_spec_filter, {_prefix}),
      _contamination_proba("contamination_proba", &_prior_contamination_proba,
                           {_prefix, ProgramOptions::FIXED_PRIOR_ON_MASS_SPEC_CONTAMINATION_PROBA})
#endif
{
	// create trees (including mu_0, mu_1 and binned branch lengths)
	_create_trees();

	// collect the parameters of every compiled-in data source
	TDataSources sources;
#ifdef USE_LOTUS
	sources.gamma      = &_gamma;
	sources.error_rate = &_error_rate;
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
	sources.epsilon_simple_model = &_epsilon_simple_model;
#endif

	// create the likelihood box that owns the Markov field and all data sources
	_data_model = std::make_unique<TDataModel>(_trees, sources, n_iterations,
	                                           _markov_field_stattools_param, _prefix, simulate);

	// note: fixed prior parameters are now passed through the TParameterDefinition at
	// construction time, because TNodeTyped's constructor immediately forwards them to the
	// prior via setFixedPriorParameters(). Setting them afterwards would be too late.

	// create (fake) observation for stattools
	_obs = std::make_unique<SpecDataObs>("data_obs", _data_model.get(), StorageDataObs(),
	                                     stattools::TObservationDefinition{});

#ifdef USE_MS_DATA
	_msms_data  = std::make_unique<TMSMSData>(_trees, _markov_field_stattools_param,
	                                          &_mass_spec_filters, &_contamination_proba);
	_msdata_obs = std::make_unique<SpecMSData>("msdata_obs", _msms_data.get(), StorageMSData(),
	                                           stattools::TObservationDefinition{});
#endif

	// define function that is called when updating
	_fun_update_mrf = &TDataModel::update_markov_field;
	stattools::instances::dagBuilder().addFuncToUpdate(*_data_model, _fun_update_mrf);
};

//--------------------------------------
// TCore
//--------------------------------------

TCore::TCore() { ProgramOptions::parse(); }

void TCore::infer() {
	_started = true;

	// read output prefix
	std::string prefix;
	if (parameters().exists("out")) {
		prefix = parameters().get<std::string>("out");
		logfile().list("Writing output to prefix '", prefix, "' (argument 'out').");
	} else {
		// falls back to whichever data file the build actually reads
		prefix =
		    coretools::str::readBeforeLast(ProgramOptions::default_prefix_source_filename(), ".");
		logfile().list("Writing output to default prefix '", prefix, "' (use '--out' to change).");
	}

	// create MCMC object and get number of iterations that will be run
	stattools::TMCMC mcmc;
	size_t n_iterations = ProgramOptions::NUM_ITERATIONS;

	// build model
	TModel model(n_iterations, prefix, false);

	// run MCMC
	mcmc.runMCMC(prefix);
}

void TCore::simulate() {
	_started            = true;
	std::string prefix  = parameters().get("out", "acol");
	size_t n_iterations = TMarkovField::get_num_iterations_simulation();

	// build model
	TModel model(n_iterations, prefix, true);

	// simulate
	stattools::TSimulator::simulate();
}
