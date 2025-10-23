/*
 * TExampleTask.cpp
 *
 *  Created on: May 15, 2023
 *      Author: mvisani
 */

#include "TCore.h"
#include "Types.h"
#include "coretools/Main/TError.h"
#include "stattools/MCMC/TMCMC.h"

using namespace coretools::instances;

//--------------------------------------
// TModel
//--------------------------------------

void TModel::_create_tree(size_t dimension, const std::string &filename, const std::string &tree_name,
                          const std::string &prefix) {

	// create mean log nu
	stattools::TRuntimeConfigParam config_mean_log_nu;
	config_mean_log_nu.set_name(tree_name + "_mean_log_nu");
	_mean_log_nu.push_back(std::make_unique<stattools::TParameter<SpecMeanLogNu, PriorOnLogNu>>(&_prior_on_mean_log_nu,
	                                                                                            config_mean_log_nu));

	// create var log nu
	stattools::TRuntimeConfigParam config_var_log_nu;
	config_var_log_nu.set_name(tree_name + "_var_log_nu");
	_var_log_nu.push_back(
	    std::make_unique<stattools::TParameter<SpecVarLogNu, PriorOnLogNu>>(&_prior_on_var_log_nu, config_var_log_nu));

	// create prior on log nu
	_prior_on_log_nu.push_back(std::make_unique<PriorOnLogNu>(_mean_log_nu.back().get(), _var_log_nu.back().get()));

	// create log nu
	stattools::TRuntimeConfigParam config_log_nu;
	config_log_nu.set_name(tree_name + "_log_nu");
	config_log_nu.excludeFromDAGUpdates(true); // never update
	_log_nu.push_back(
	    std::make_unique<stattools::TParameter<SpecLogNu, TTree>>(_prior_on_log_nu.back().get(), config_log_nu));

	// create alpha
	stattools::TRuntimeConfigParam config_alpha;
	config_alpha.set_name(tree_name + "_alpha");
	config_alpha.excludeFromDAGUpdates(true); // never update
	_alpha.push_back(std::make_unique<stattools::TParameter<SpecAlpha, TTree>>(&_prior_on_alpha, config_alpha));

	// create branch lengths
	stattools::TRuntimeConfigParam config_branch_lengths;
	config_branch_lengths.set_name(tree_name + "_branch_lengths");
	config_branch_lengths.excludeFromDAGUpdates(true); // never update
	config_branch_lengths.setPrefix(prefix + "_" + tree_name);
	_binned_branch_lengths.push_back(std::make_unique<stattools::TParameter<SpecBinnedBranches, TTree>>(
	    &_prior_on_binned_branch_lengths, config_branch_lengths));

	// create tree
	_trees.emplace_back(std::make_unique<TTree>(dimension, filename, tree_name, _alpha.back().get(),
	                                            _log_nu.back().get(), _binned_branch_lengths.back().get()));

	// create markov field (only for stattools purposes such that a valid DAG can be built)
	stattools::TRuntimeConfigParam config_markov_field;
	config_markov_field.set_name(tree_name + "_MRF");
	config_markov_field.excludeFromDAGUpdates(true); // never update
	_markov_field_stattools_param.emplace_back(
	    std::make_unique<stattools::TParameter<SpecMarkovField, TLotus>>(_trees.back().get(), config_markov_field));
}

void TModel::_create_trees(const std::string &prefix) {
	using namespace coretools::instances;
	// read filenames
	std::string filename_tree_species   = parameters().get("tree_species");
	std::string filename_tree_molecules = parameters().get("tree_molecules");
	std::vector<std::string> filenames_tree_others;
	if (parameters().exists("tree_others")) { parameters().fill("tree_others", filenames_tree_others); }

	size_t num_trees = 2 + filenames_tree_others.size();
	_trees.reserve(num_trees);

	// first tree: molecules
	_create_tree(0, filename_tree_molecules, "molecules", prefix);

	// middle trees: all others (e.g. tissues)
	for (size_t i = 1; i < num_trees - 1; ++i) {
		std::string name = coretools::str::split(filenames_tree_others[i - 1], ':');
		if (name.empty()) {
			UERROR("Argument 'tree_others': Please provide a name for each other tree, separated by a : from the "
			       "filename (e.g. myTreeName:pathToFile)");
		}
		_create_tree(i, filenames_tree_others[i - 1], name, prefix);
	}
	// last tree: species
	_create_tree(num_trees - 1, filename_tree_species, "species", prefix);

	for (auto &tree : _trees) { tree->initialize_cliques_and_Z(_trees); }
}

TModel::TModel(size_t n_iterations, const std::string &prefix, bool simulate) {
	// create trees (including mu_0, mu_1 and binned branch lengths)
	_create_trees(prefix);

	// create lotus
	_lotus = std::make_unique<TLotus>(_trees, &_gamma, &_error_rate, n_iterations, _markov_field_stattools_param,
	                                  prefix, simulate);
	_error_rate.getConfig().setPriorParameters("1.5");
	for (auto &it : _var_log_nu) { it->getConfig().setPriorParameters("2.0"); }

	// create (fake) observation for stattools
	_obs = std::make_unique<SpecLotus>(_lotus.get(), StorageLotus(), stattools::TRuntimeConfigObs());

	// define function that is called when updating
	_fun_update_mrf = &TLotus::update_markov_field;
	stattools::instances::dagBuilder().addFuncToUpdate(*_lotus, _fun_update_mrf);
};

//--------------------------------------
// TCore
//--------------------------------------

TCore::TCore() {
	NUMBER_OF_THREADS              = coretools::getNumThreads();
	SIMULATION_NO_Z_INITIALIZATION = coretools::instances::parameters().exists("simulation_no_Z_initilisation");
	SIMULATION_NO_Y_INITIALIZATION = coretools::instances::parameters().exists("simulation_no_Y_initilisation");
	WRITE_Y                        = coretools::instances::parameters().exists("write_Y");
	WRITE_Y_TRACE                  = coretools::instances::parameters().exists("write_Y_trace");
	WRITE_Z                        = coretools::instances::parameters().exists("write_Z");
	WRITE_Z_TRACE                  = coretools::instances::parameters().exists("write_Z_trace");
	WRITE_JOINT_LOG_PROB_DENSITY   = coretools::instances::parameters().exists("write_joint_log_prob_density");
	WRITE_BRANCH_LENGTHS           = coretools::instances::parameters().exists("write_branch_lengths");
}

void TCore::infer() {
	// read filename of lotus
	std::string filename_lotus = TLotus::get_filename_lotus();

	// read output prefix
	std::string prefix;
	if (parameters().exists("out")) {
		prefix = parameters().get<std::string>("out");
		logfile().list("Writing output to prefix '", prefix, "' (argument 'out').");
	} else {
		prefix = coretools::str::readBeforeLast(filename_lotus, ".");
		logfile().list("Writing output to default prefix '", prefix, "' (use 'out' to change).");
	}

	// create MCMC object and get number of iterations that will be run
	stattools::TMCMC mcmc;
	size_t n_iterations = mcmc.get_num_iterations();

	// build model
	TModel model(n_iterations, prefix, false);

	// run MCMC
	mcmc.runMCMC(prefix);
}

void TCore::simulate() {
	std::string prefix  = parameters().get("out", "acol");
	size_t n_iterations = TMarkovField::get_num_iterations_simulation();

	// build model
	TModel model(n_iterations, prefix, true);

	// simulate
	stattools::TSimulator::simulate(prefix);
}
