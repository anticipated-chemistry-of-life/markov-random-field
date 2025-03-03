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
#include "stattools/DAG/TDAGBuilder.h"
#include "stattools/MCMC/TMCMC.h"
#include <memory>

//--------------------------------------
// TModel
//--------------------------------------

class TModel {
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
	stattools::TParameter<SpecGamma, TLotus> _gamma{&_prior_on_gamma};

	// observation
	std::unique_ptr<TLotus> _lotus;
	std::unique_ptr<SpecLotus> _obs;

	// functions that are called when updating
	void (TMarkovField::*_fun_update_mrf)(size_t);

	void _create_tree(size_t dimension, const std::string &filename, const std::string &tree_name);
	void _create_trees();

public:
	TModel(size_t n_iterations, const std::string &prefix);
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
	TTask_infer() : coretools::TTask("Inferring metabolite presence in all the species there are on this planet!!") {};

	void run() override {
		TCore core;
		core.infer();
	};
};

class TTask_simulate : public coretools::TTask {
public:
	TTask_simulate()
	    : coretools::TTask("Simulating metabolite presence in all the species there are on this planet!!") {};

	void run() override {
		TCore core;
		core.simulate();
	};
};

#endif /* TEXAMPLETASK_H_ */
