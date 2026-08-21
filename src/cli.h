#pragma once

#include "coretools/Main/TError.h"
#include "coretools/Main/TParameters.h"
#include <cstddef>
#include <iostream>
#include <string_view>

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

	static inline double EPSILON = 0.001;

	static inline size_t BRANCH_LENGTHS_BINS = 100;

	static inline size_t SHEET_SIZE_K = static_cast<size_t>(1e7);

	static inline bool FIX_Y = false;

	static inline bool FIX_Z = false;

	static inline std::string LOTUS_FILENAME = "lotus.tsv";

	static inline std::string SIMPLE_DATA_FILENAME = "simple_data.tsv";

	/// Error rate of the simple error model: the probability that a cell of the observed matrix D
	/// reports the opposite of the latent state Y. Used as the simulated truth during simulation
	/// and as the starting value of the inferred epsilon_simple_model parameter during inference.
	static inline double EPSILON_SIMPLE_MODEL = 0.1;

	static inline double GAMMA = 1.1;

	static inline double ALPHA = 0.5;

	static inline double LOG_NU_C = 0.0;

	/// Standard deviation of the Normal that jitters the per-clique log_nu initial values around
	/// LOG_NU_C. Must be > 0: if every log_nu starts identical, the MLE that seeds var_log_nu is 0,
	/// producing a degenerate (zero-variance) prior that freezes log_nu, mean_log_nu and var_log_nu
	/// (proposals keep being made but are never accepted).
	static inline double LOG_NU_C_INIT_SD = 0.1;

	/// the index of the probability to pass the MS filter. Actual probabilities are found in
	/// constants.h
	static inline size_t INDEX_PROBA_TO_PASS_MS_FILTER = 128;

	static inline double PROBA_OF_MS_CONTAMINATION = 0.001;

	/// Beta: the probability with which the assignment update proposes to move an already assigned
	/// feature back to the unknown molecule. With probability 1 - beta the update instead picks a
	/// molecule among the candidates of that feature. Must be strictly inside (0, 1): at 0 no
	/// feature ever returns to the unknown molecule (and the reverse move has an infinite Hastings
	/// ratio), at 1 no feature ever moves between two real molecules.
	static inline double MS_PROBA_MOVE_TO_UNKNOWN = 0.1;

	static inline std::string_view FIXED_PRIOR_ON_MASS_SPEC_CONTAMINATION_PROBA = "0.3,5.0";

	static inline std::string_view FIXED_PRIOR_ON_VAR_LOG_NU = "1.0,2.0";

	static inline size_t NUM_ITERATIONS = 20000;

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

		EPSILON = params.get<double>("epsilon", EPSILON);

		BRANCH_LENGTHS_BINS = params.get<size_t>("n_bins", BRANCH_LENGTHS_BINS);

		SHEET_SIZE_K = params.get("K", SHEET_SIZE_K);

		FIX_Y = !params.get("Y.update", true);
		FIX_Z = !params.get("Z.update", true);

		LOTUS_FILENAME = params.get("lotus", LOTUS_FILENAME);

		SIMPLE_DATA_FILENAME = params.get("simple_data", SIMPLE_DATA_FILENAME);

		EPSILON_SIMPLE_MODEL = params.get<double>("epsilon_simple_model", EPSILON_SIMPLE_MODEL);
		// The parameter is ZeroOneOpen; letting an out-of-range value through here would surface as
		// an opaque throw from deep inside stattools instead of a message the user can act on.
		if (EPSILON_SIMPLE_MODEL <= 0.0 || EPSILON_SIMPLE_MODEL >= 1.0) {
			throw coretools::TUserError("--epsilon_simple_model must be strictly between 0 and 1, "
			                            "but got ",
			                            EPSILON_SIMPLE_MODEL, ".");
		}

		GAMMA = params.get<double>("gamma", GAMMA);

		ALPHA = params.get<double>("alpha", ALPHA);

		LOG_NU_C = params.get<double>("log_nu_c", LOG_NU_C);

		LOG_NU_C_INIT_SD = params.get<double>("log_nu_c_init_sd", LOG_NU_C_INIT_SD);

		NUM_ITERATIONS = params.get<size_t>("iterations", NUM_ITERATIONS);

		MS_PROBA_MOVE_TO_UNKNOWN =
		    params.get<double>("ms_proba_move_to_unknown", MS_PROBA_MOVE_TO_UNKNOWN);
		if (MS_PROBA_MOVE_TO_UNKNOWN <= 0.0 || MS_PROBA_MOVE_TO_UNKNOWN >= 1.0) {
			throw coretools::TUserError("--ms_proba_move_to_unknown must be strictly between 0 and "
			                            "1, but got ",
			                            MS_PROBA_MOVE_TO_UNKNOWN, ".");
		}
	}

	/// The file the default output prefix is derived from when '--out' is not given. Which data
	/// source that is depends on what was compiled in.
	[[nodiscard]] static std::string default_prefix_source_filename() {
#ifdef USE_LOTUS
		return LOTUS_FILENAME;
#else
		return SIMPLE_DATA_FILENAME;
#endif
	}

	static void printHelp() {
		std::cout << "--numThreads                   Set the number of threads you want to use\n";
		std::cout << "--write_Y                      Write Y output\n";
		std::cout << "--write_Y_trace                Write Y trace\n";
		std::cout << "--write_Z                      Write Z output\n";
		std::cout << "--write_Z_trace                Write Z trace\n";
		std::cout << "--write_branch_lengths         Output branch lengths\n";
		std::cout << "--lotus                        LOTUS data file (only with the 'lotus' build "
		             "option)\n";
		std::cout << "--simple_data                  Simple error model data file (only with the "
		             "'simple_data' build option)\n";
		std::cout << "--epsilon_simple_model         Error rate of the simple error model, in "
		             "(0,1). Simulated truth when simulating, starting value when inferring\n";
		std::cout << "--ms_proba_move_to_unknown     Probability of proposing to move an assigned "
		             "MS feature back to the unknown molecule (in (0,1))\n";
	}
};
