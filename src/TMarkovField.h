#pragma once

#include "storages/storage_concepts.h"
#include "tree/TTree.h"

class TDataModel; // forward declaration

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
