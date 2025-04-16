/*
 * main.cpp
 */

#include "coretools/Main/TMain.h"

//---------------------------------------------------------------------------
// Includes for tasks
//---------------------------------------------------------------------------

#include "TCore.h"

//---------------------------------------------------------------------------
// Existing Tasks
//---------------------------------------------------------------------------

void addTask(coretools::TMain &main) {
	// Tasks consist of a name and a pointer to a TTask object.
	// Use main.addRegularTask() to add a regular task (shown in
	// list of available tasks) Use main.addDebugTask() to add a
	// debug task (not shown in list of available tasks)

	main.addRegularTask("infer", new TTask_infer());
	main.addRegularTask("simulate", new TTask_simulate());
};

//---------------------------------------------------------------------------
// Existing Integration tests
//---------------------------------------------------------------------------

void addTests(coretools::TMain &) {
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
	                      "https://github.com/anticipated-chemistry-of-life/markov-random-field",
	                      "marco.visani@unifr.ch");

	// add existing tasks and tests
	addTask(main);
	addTests(main);

	// now run program
	return main.run(argc, argv);
};
