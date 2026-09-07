#pragma once

#include "tree/TTree.h"

class TDataModel; // forward declaration

/// Per-cell outcome of one Y update. Each data source keeps its own likelihood bookkeeping (they
/// are independent terms), so the results are handed back separately instead of merged.
struct TYUpdateResult {
	int diff_counter_1_in_last_dim = 0;
#ifdef USE_LOTUS
	/// P(L_cell | x = new_state). Neutral value 1.0 (log 0), used when collapsing makes the LOTUS
	/// term identical for both Y states and calculate_LL_update_Y leaves it untouched.
	double prob_lotus_new_state = 1.0;
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
	/// Whether the observed D cell contradicts the state Y was just set to.
	bool simple_model_disagrees = false;
#endif
};

class TMarkovField {
private:
	// trees and Y
	std::vector<std::unique_ptr<TTree>> &_trees;
	TFieldStorage _Y;
	std::string _prefix;

	// fix values?
	bool _fix_Y = false;
	bool _fix_Z = false;
	// stuff for updating Y
	IndexArray _num_leaves_per_dim_except_last{};

private:
	template<bool IsSimulation, bool initYFromData>
	TYUpdateResult _update_Y(const IndexArray &index_in_leaves_space,
	                         const TDataModel &data_model) {
		// prepare log probabilities for the two possible states
		std::array<coretools::TSumLogProbability, 2> sum_log;

		// calculate probabilities in Markov random field
		if constexpr (!initYFromData) { _calculate_log_prob_field(index_in_leaves_space, sum_log); }
		std::array<coretools::TSumLogProbability, 2> sum_log_field = sum_log;

		// Declared outside the IsSimulation branch so the simulation instantiation does not warn
		// about an unused variable. 1.0 is the neutral value: calculate_LL_update_Y leaves prob
		// untouched when collapsing makes the LOTUS term identical for both states, and adding
		// log(1) = 0 is exactly the no-op that case needs.
#ifdef USE_LOTUS
		std::array<double, 2> prob_lotus{1.0, 1.0};
#endif
		if constexpr (!IsSimulation) {
			// calculate log likelihood (lotus)
#ifdef USE_LOTUS
			_calc_lotus_LL(index_in_leaves_space, prob_lotus, data_model);
			for (size_t i = 0; i < 2; ++i) { sum_log[i].add(prob_lotus[i]); }
#endif
			// calculate log likelihood (simple error model)
#ifdef USE_SIMPLE_ERROR_MODEL
			std::array<double, 2> prob_simple{};
			_calc_simple_error_model_LL(index_for_tmp_state, prob_simple, data_model);
			for (size_t i = 0; i < 2; ++i) { sum_log[i].add(prob_simple[i]); }
#endif
			// calculate log likelihood mass spec data
			if (_ms_data.has_value()) { _ms_data->add_log_likelihood(index_copy, sum_log); }
		}

		// sample state
		const bool new_state = sample(sum_log);
	}

public:
	TMarkovField(size_t n_iterations, std::vector<std::unique_ptr<TTree>> &Trees,
	             std::string _prefix);
	~TMarkovField() = default;

	// updates
	void update(TDataModel &data_model, size_t iteration);

	// simulation
	void simulate(TDataModel &data_model);

	// get Y
	[[nodiscard]] const TFieldStorage &get_Y_matrix() const;

	// functions to perform stuff on Y after burnin / MCMC finished
	void burninHasFinished();
	void MCMCHasFinished();
	void oneBurninHasFinished();

	static size_t get_num_iterations_simulation() { return ProgramOptions::NUM_ITERATIONS; }
};
