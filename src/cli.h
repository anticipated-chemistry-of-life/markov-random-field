#pragma once

#include "coretools/Main/TParameters.h"
#include <cstddef>
#include <iostream>

class ProgramOptions {
public:
	static inline size_t NUMBER_OF_THREADS = 1;

	static inline bool SIMULATION_NO_Z_INITIALIZATION = false;
	static inline bool SIMULATION_NO_Y_INITIALIZATION = false;

	static inline bool WRITE_Y       = false;
	static inline bool WRITE_Y_TRACE = false;

	static inline bool WRITE_Z       = false;
	static inline bool WRITE_Z_TRACE = false;

	static inline bool WRITE_JOINT_LOG_PROB_DENSITY = false;

	static inline bool WRITE_BRANCH_LENGTHS = false;

	static inline double EPSILON = 0.0001;

	static inline size_t BRANCH_LENGTHS_BINS = 100;

	static inline size_t SHEET_SIZE_K = (size_t)(1e7);

	static inline bool FIX_Y = false;

	static inline bool FIX_Z = false;

	static inline std::string LOTUS_FILENAME = "lotus.tsv";

	static inline double GAMMA = 1.5;

	static void parse() {

		auto &params = coretools::instances::parameters();

		NUMBER_OF_THREADS = coretools::getNumThreads();

		SIMULATION_NO_Z_INITIALIZATION = params.exists("simulation_no_Z_initilisation");

		SIMULATION_NO_Y_INITIALIZATION = params.exists("simulation_no_Y_initilisation");

		WRITE_Y       = params.exists("write_Y");
		WRITE_Y_TRACE = params.exists("write_Y_trace");

		WRITE_Z       = params.exists("write_Z");
		WRITE_Z_TRACE = params.exists("write_Z_trace");

		WRITE_JOINT_LOG_PROB_DENSITY = params.exists("write_joint_log_prob_density");

		WRITE_BRANCH_LENGTHS = params.exists("write_branch_lengths");

		EPSILON = params.get<double>("epsilon", 0.0001);

		BRANCH_LENGTHS_BINS = params.get("n_bins", 100);

		SHEET_SIZE_K = params.get("K", size_t(1e7));

		FIX_Y = !params.get("Y.update", true);
		FIX_Z = !params.get("Z.update", true);

		LOTUS_FILENAME = params.get("lotus");

		GAMMA = params.get<double>("gamma", 1.5);
	}

	static void printHelp() {
		std::cout << "--numThreads                   Set the number of threads you want to use\n";
		std::cout << "--write_Y                      Write Y output\n";
		std::cout << "--write_Y_trace                Write Y trace\n";
		std::cout << "--write_Z                      Write Z output\n";
		std::cout << "--write_Z_trace                Write Z trace\n";
		std::cout << "--write_branch_lengths         Output branch lengths\n";
	}
};
