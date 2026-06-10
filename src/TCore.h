#ifndef TEXAMPLETASK_H_
#define TEXAMPLETASK_H_

#include "TLotus.h"
#include "TMarkovField.h"
#include "Types.h"
#include "cli.h"
#include "coretools/Main/TTask.h"
#include "stattools/ParametersObservations/TParameter.h"
#include "tree/TTree.h"
#include <memory>

//--------------------------------------
// TModel
//--------------------------------------

class TModel {
private:
	// mean log nu
	PriorOnMeanLogNu _prior_on_mean_log_nu{};
	std::vector<std::unique_ptr<stattools::TParameter<SpecMeanLogNu, PriorOnLogNu>>> _mean_log_nu;

	// var log nu
	PriorOnVarLogNu _prior_on_var_log_nu{};
	std::vector<std::unique_ptr<stattools::TParameter<SpecVarLogNu, PriorOnLogNu>>> _var_log_nu;

	// Now we can make the log nu
	std::vector<std::unique_ptr<PriorOnLogNu>> _prior_on_log_nu{};
	std::vector<std::unique_ptr<stattools::TParameter<SpecLogNu, TTree>>> _log_nu;

	// alpha
	PriorOnAlpha _prior_on_alpha{};
	std::vector<std::unique_ptr<stattools::TParameter<SpecAlpha, TTree>>> _alpha{};

	// binned branch lengths
	PriorOnBinnedBranches _prior_on_binned_branch_lengths{};
	std::vector<std::unique_ptr<stattools::TParameter<SpecBinnedBranches, TTree>>>
	    _binned_branch_lengths;

	// trees
	std::vector<std::unique_ptr<TTree>> _trees;

	// Markov field parameters (only needed for stattools)
	std::vector<std::unique_ptr<stattools::TParameter<SpecMarkovField, TLotus>>>
	    _markov_field_stattools_param;

	// gamma
	PriorOnGamma _prior_on_gamma{};
	stattools::TParameter<SpecGamma, TLotus> _gamma{&_prior_on_gamma};

	// error rate
	PriorOnErrorRate _prior_on_error_rate{};
	stattools::TParameter<SpecErrorRate, TLotus> _error_rate{&_prior_on_error_rate};

	// observation
	std::unique_ptr<TLotus> _lotus;
	std::unique_ptr<SpecLotus> _obs; // "fake" observation, only needed for stattools

#ifdef USE_MS_DATA
	// Mass spec stuff
	// The probability to pass a filter
	PriorOnMassSpecFilter _prior_on_mass_spec_filter{};
	stattools::TParameter<SpecMassSpecFilter, TMSMSData> _mass_spec_filters{
	    &_prior_on_mass_spec_filter};

	// Mass spec contamination parameters
	PriorOnContaminationProba _prior_contamination_proba{};
	stattools::TParameter<SpecContaminationProba, TMSMSData> _contamination_proba{
	    &_prior_contamination_proba};

	std::unique_ptr<TMSMSData> _msms_data;
	std::unique_ptr<SpecMSData> _msdata_obs; // "fake" observation, only needed for stattools
#endif

	// functions that are called when updating
	void (TLotus::*_fun_update_mrf)(size_t);

	void _create_tree(size_t dimension, const std::string &filename, const std::string &tree_name,
	                  const std::string &prefix);
	void _create_trees(const std::string &prefix);

public:
	TModel(size_t n_iterations, const std::string &prefix, bool simulate);
};

//--------------------------------------
// TCore
//--------------------------------------

class TCore {
private:
public:
	TCore();
	~TCore() = default;
	void infer();
	void simulate();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_infer : public coretools::TTask {
public:
	TTask_infer()
	    : coretools::TTask(
	          "Inferring metabolite presence in all the species there are on this planet!!") {};

	void run() override {
		TCore core;
		if (coretools::instances::parameters().exists("help")) {
			ProgramOptions::printHelp();
			return;
		}
		core.infer();
	};
};

class TTask_simulate : public coretools::TTask {
public:
	TTask_simulate()
	    : coretools::TTask(
	          "Simulating metabolite presence in all the species there are on this planet!!") {};

	void run() override {
		TCore core;
		if (coretools::instances::parameters().exists("help")) {
			ProgramOptions::printHelp();
			return;
		}
		core.simulate();
	};
};

#endif /* TEXAMPLETASK_H_ */
