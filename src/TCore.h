/*
 * TExampleTask.h
 *
 *  Created on: Nov 30, 2020
 *      Author: phaentu
 */

#ifndef TEXAMPLETASK_H_
#define TEXAMPLETASK_H_

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
#include "TMarkovField.h"

class TCore {
private:
	void (TMarkovField:: *_fun_update_mrf)();

public:
	void infer(){
		// define function that is called when updating
		//_fun_update_mrf = &TMarkovField::update_markov_field;
		//stattools::instances::dagBuilder().addFuncToUpdate(*_lotus, _funUpdateLotus);
	};
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_infer : public coretools::TTask {
public:
	// constructor must fill explanation shown to users
	TTask_infer() : coretools::TTask("Inferring metabolite presence in all the species there are on this planet!!"){};

	// a task must overload the run function that takes two arguments:
	// coretools::TParameters & Parameters, coretools::TLog* Logfile Usually, a
	// task creates an object and calls its function
	void run() override {
		TCore core;
		core.infer();
	};
};

#endif /* TEXAMPLETASK_H_ */
