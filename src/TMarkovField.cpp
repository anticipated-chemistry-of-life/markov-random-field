#include "TMarkovField.h"
#include "TDataModel.h"

TMarkovField::TMarkovField(size_t n_iterations, std::vector<std::unique_ptr<TTree>> &Trees,
                           std::string _prefix)
    : _trees(Trees), _prefix(std::move(_prefix)) {
	using namespace coretools::instances;

	// find molecule and species dimensions; construct mass spec data if both trees are present
	// _ms_data.emplace(_trees); // TODO: once we have data, we can remove this

	// read: fix Y or Z?
	_fix_Y = ProgramOptions::FIX_Y;
	if (_fix_Y) { logfile().list("Will fix Y during the MCMC."); }
	_fix_Z = ProgramOptions::FIX_Z;
	if (_fix_Z) { logfile().list("Will fix Z during the MCMC."); }

	// set number of leaves per dimension (set the last dimension to one)
	_num_leaves_per_dim_except_last.fill(1);
	for (size_t i = 0; i < _trees.size() - 1; ++i) {
		_num_leaves_per_dim_except_last[i] = _trees[i]->get_number_of_leaves();
	}

	// initialize Y
	auto num_leaves_per_dim   = _num_leaves_per_dim_except_last;
	num_leaves_per_dim.back() = _trees.back()->get_number_of_leaves();
	_Y.initialize(n_iterations, num_leaves_per_dim);
	if (parameters().exists("set_Y")) {
		std::string filename = parameters().get("set_Y", "acol_simulated_Y.txt");
		// _read_Y_from_file(filename);
	}
}

void TMarkovField::oneBurninHasFinished() { _Y.remove_zeros(); }
void TMarkovField::burninHasFinished() {
	_Y.reset_counts();
	_Y.remove_zeros();
}

void TMarkovField::MCMCHasFinished() {
	// write function to write the posterior state of Y to file
	// _write_Y_to_file<false>(_prefix + "_Y_posterior.txt");
}
