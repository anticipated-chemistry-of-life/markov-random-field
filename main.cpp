/*
 * main.cpp
 */

#include "coretools/Main/TMain.h"

//---------------------------------------------------------------------------
// Includes for tasks
//---------------------------------------------------------------------------

#include "TCore.h"
#include "TTree.h"

//---------------------------------------------------------------------------
// Existing Tasks
//---------------------------------------------------------------------------

void addTask(coretools::TMain &main) {
	// Tasks consist of a name and a pointer to a TTask object.
	// Use main.addRegularTask() to add a regular task (shown in
	// list of available tasks) Use main.addDebugTask() to add a
	// debug task (not shown in list of available tasks)

	// TODO: add all your tasks here (infer and simulate)
	main.addRegularTask("infer", new TTask_infer());
};

//---------------------------------------------------------------------------
// Existing Integration tests
//---------------------------------------------------------------------------

void addTests(coretools::TMain &){
    // Use main.addTest to add integration tests

    // Use main.addTestSuite to add test suites
};

//---------------------------------------------------------------------------
// Main function
//---------------------------------------------------------------------------

int main(int argc, char *argv[]) {
	// Create main by providing a program name, a version, an
	// affiliation, link to repo and contact email
	coretools::TMain main("acol", "0.1", "University of Fribourg",
	                      "https://bitbucket.org/wegmannlab/metabolite_inference", "marco.visani@unifr.ch");

	// add existing tasks and tests
	addTask(main);
	addTests(main);

	// now run program
	return main.run(argc, argv);
};
