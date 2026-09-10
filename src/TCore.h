#ifndef TEXAMPLETASK_H_
#define TEXAMPLETASK_H_

#include "TDataModel.h"
#include "TMarkovField.h"
#include "Types.h"
#include "cli.h"
#include "coretools/Main/TTask.h"
#include "stattools/ParametersObservations/TParameter.h"
#include "tree/TTree.h"
#include <memory>
#include <string>

//--------------------------------------
// TModel
//--------------------------------------

class TModel {
private:
	// Declared first so it can be used to initialize the parameters below: their construction order
	// follows declaration order, and the #ifdef'd parameter blocks must not have to carry it.
	const std::string _prefix;

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
	MarkovFieldParams _markov_field_stattools_param;

	// The error probability standing between the two tree fields and the field. Not behind an
	// #ifdef: the link stands there in every build, whichever data sources were compiled in.
	PriorOnErrorProbability _prior_on_error_probability{};
	TMarkovField::TypeParamErrorProbability _error_probability;

#ifdef USE_LOTUS
	// gamma
	PriorOnGamma _prior_on_gamma{};
	// note: initialized in the TModel constructor so the definition can be given the global
	// output prefix (needed to enable mean/var storage and trace files)
	TLotus::TypeParamGamma _gamma;

	// error rate
	PriorOnErrorRate _prior_on_error_rate{};
	TLotus::TypeParamErrorRate _error_rate;
#endif

#ifdef USE_SIMPLE_ERROR_MODEL
	// error rate of the simple error model
	PriorOnEpsilonSimpleModel _prior_on_epsilon_simple_model{};
	TSimpleErrorModel::TypeParamEpsilon _epsilon_simple_model;
#endif

	// observation
	std::unique_ptr<TDataModel> _data_model;
	std::unique_ptr<SpecDataObs> _obs; // "fake" observation, only needed for stattools

#ifdef USE_MS_DATA
	// Mass spec stuff
	// The probability to pass a filter
	PriorOnMassSpecFilter _prior_on_mass_spec_filter{};
	stattools::TParameter<SpecMassSpecFilter, TMSMSData> _mass_spec_filters;

	// Mass spec contamination parameters
	PriorOnContaminationProba _prior_contamination_proba{};
	stattools::TParameter<SpecContaminationProba, TMSMSData> _contamination_proba;

	std::unique_ptr<TMSMSData> _msms_data;
	std::unique_ptr<SpecMSData> _msdata_obs; // "fake" observation, only needed for stattools
#endif

	// functions that are called when updating
	void (TDataModel::*_fun_update_mrf)();

	void _create_tree(size_t dimension, const std::string &filename, const std::string &tree_name);
	void _create_trees();

public:
	TModel(size_t n_iterations, const std::string &prefix, bool simulate);
};

//--------------------------------------
// TCore
//--------------------------------------

class TCore {
private:
	bool _started = false;

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
