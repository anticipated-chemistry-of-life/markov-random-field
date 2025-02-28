/*
 * TExampleTask.h
 *
 *  Created on: Nov 30, 2020
 *      Author: phaentu
 */

#ifndef TEXAMPLETASK_H_
#define TEXAMPLETASK_H_

#include "TMarkovField.h"
#include "TStorageYVector.h"
#include "TTree.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TTask.h"
#include "coretools/Types/TStringHash.h"
#include "stattools/DAG/TDAGBuilder.h"
#include "stattools/ParametersObservations/TNodeBase.h"
#include "stattools/ParametersObservations/TObservation.h"
#include "stattools/ParametersObservations/TParameter.h"
#include "stattools/ParametersObservations/spec.h"
#include "stattools/Priors/TPriorBernoulli.h"
#include "stattools/Priors/TPriorExponential.h"
#include "stattools/Priors/TPriorNormal.h"
#include "stattools/Priors/TPriorPoisson.h"
#include <memory>

//--------------------------------------
// TModel
//--------------------------------------

template<bool SimpleErrorModel> class TModel {
private:
	// mu_0 and mu_1
	PriorOnMu _prior_on_mu{};
	std::vector<std::unique_ptr<stattools::TParameter<SpecMu_0, TTree>>> _mu_0;
	std::vector<std::unique_ptr<stattools::TParameter<SpecMu_1, TTree>>> _mu_1;

	// binned branch lengths
	PriorOnBinnedBranches _prior_on_binned_branch_lengths{};
	std::vector<std::unique_ptr<stattools::TParameter<SpecBinnedBranches, TTree>>> _binned_branch_lengths;

	// trees
	std::vector<std::unique_ptr<TTree>> _trees;

	// Markov field
	std::unique_ptr<TMarkovField> _markov_field;

	// gamma
	PriorOnGamma _prior_on_gamma{};
	stattools::TParameter<SpecGamma, TLotus<SimpleErrorModel>> _gamma{&_prior_on_gamma};

	// observation
	std::unique_ptr<TLotus<SimpleErrorModel>> _lotus;
	using SpecLotus =
	    stattools::TObservation<TypeLotus, stattools::name("lotus"), NumDimLotus, TLotus<SimpleErrorModel>>;
	std::unique_ptr<SpecLotus> _obs;

	// functions that are called when updating
	void (TMarkovField::*_fun_update_mrf)();

	void _create_tree(size_t dimension, std::string_view filename, std::string_view name) {
		// TODO: add possibility to provide custom names such that mu's of different trees have a different name!
		// create mu_0, mu_1 and binned branch lengths
		_mu_0.push_back(std::make_unique<stattools::TParameter<SpecMu_0, TTree>>(&_prior_on_mu));
		_mu_1.push_back(std::make_unique<stattools::TParameter<SpecMu_1, TTree>>(&_prior_on_mu));
		_binned_branch_lengths.push_back(
		    std::make_unique<stattools::TParameter<SpecBinnedBranches, TTree>>(&_prior_on_binned_branch_lengths));

		// create tree
		_trees.emplace_back(std::make_unique<TTree>(dimension, filename, name, _mu_0.back().get(), _mu_1.back().get(),
		                                            _binned_branch_lengths.back().get()));
	}

	void _create_trees() {
		using namespace coretools::instances;
		// read filenames
		std::string filename_tree_species   = parameters().get("tree_species");
		std::string filename_tree_molecules = parameters().get("tree_molecules");
		std::vector<std::string> filenames_tree_others;
		if (parameters().exists("tree_others")) { parameters().fill("tree_others", filenames_tree_others); }

		size_t num_trees = 2 + filenames_tree_others.size();
		_trees.reserve(num_trees);

		// first tree: molecules
		_create_tree(0, filename_tree_molecules, "molecules");

		// middle trees: all others (e.g. tissues)
		for (size_t i = 1; i < num_trees - 1; ++i) {
			std::string name = coretools::str::split(filenames_tree_others[i - 1], ':');
			if (name.empty()) {
				UERROR("Argument 'tree_others': Please provide a name for each other tree, separated by a : from the "
				       "filename (e.g. myTreeName:pathToFile)");
			}
			_create_tree(i, filenames_tree_others[i - 1], name);
		}
		// last tree: species
		_create_tree(num_trees - 1, filename_tree_species, "species");
	}

public:
	TModel(size_t n_iterations, const std::string& prefix) {
		// create trees (including mu_0, mu_1 and binned branch lengths)
		_create_trees();

		// create markov field
		_markov_field = std::make_unique<TMarkovField>(n_iterations, _trees);

		// create lotus
		_lotus = std::make_unique<TLotus<SimpleErrorModel>>(_trees, &_gamma, _markov_field->get_Y(), prefix);

		// create (fake) observation for stattools
		_obs = std::make_unique<SpecLotus>(_lotus.get(), StorageLotus(), {});

		// define function that is called when updating
		_fun_update_mrf = &TMarkovField::update_markov_field;
		stattools::instances::dagBuilder().addFuncToUpdate(*_lotus, _fun_update_mrf);
	};
};

//--------------------------------------
// TCore
//--------------------------------------

class TCore {
private:
public:
	void infer();
	void simulate();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_infer : public coretools::TTask {
public:
	// constructor must fill explanation shown to users
	TTask_infer() : coretools::TTask("Inferring metabolite presence in all the species there are on this planet!!") {};

	// a task must overload the run function that takes two arguments:
	// coretools::TParameters & Parameters, coretools::TLog* Logfile Usually, a
	// task creates an object and calls its function
	void run() override {
		TCore core;
		core.infer();
	};
};

class TTask_simulate : public coretools::TTask {
public:
	// constructor must fill explanation shown to users
	TTask_simulate()
	    : coretools::TTask("Simulating metabolite presence in all the species there are on this planet!!") {};

	// a task must overload the run function that takes two arguments:
	// coretools::TParameters & Parameters, coretools::TLog* Logfile Usually, a
	// task creates an object and calls its function
	void run() override {
		TCore core;
		core.simulate();
	};
};

#endif /* TEXAMPLETASK_H_ */
