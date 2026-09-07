#include "TMarkovField.h"
#include "TDataModel.h"
#include <stdexcept>

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
	_num_leaves_per_dim.fill(1);
	for (size_t i = 0; i < _trees.size(); ++i) {
		_num_leaves_per_dim[i] = _trees[i]->get_number_of_leaves();
	}

	// initialize Y
	_Y.initialize(n_iterations, _num_leaves_per_dim);
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
	_write_Y_to_file<false>(_prefix + "_Y_posterior.txt");
}

void TMarkovField::update(TDataModel &data_model, size_t iteration) {
	if (ProgramOptions::WRITE_JOINT_LOG_PROB_DENSITY && iteration == 0 &&
	    !_joint_density_file.isOpen()) {
		_joint_density_file.open(_prefix + "_simulated_joint_density.txt",
		                         {
		                             "joint_density",
		                         },
		                         "\t");
	}

	if (iteration == 0 && !_z_initialized_from_children) {
		_update_all_Y<false, true>(data_model, iteration);
		for (auto &tree : _trees) { tree->initialize_Z_from_children(_Y); }
		_z_initialized_from_children = true;
	} else {
		_update_all_Y<false, false>(data_model, iteration);
	}
	if (_fix_Z) {
		_update_all_Z<false, true>(iteration);
	} else {
		_update_all_Z<false, false>(iteration);
	}
	_Y.add_to_counter(iteration);
	// calculate joint density
	if (ProgramOptions::WRITE_JOINT_LOG_PROB_DENSITY && iteration % _Y.get_thinning_factor() == 0) {
		auto sum_log_field = _calculate_complete_joint_density();
		_joint_density_file.writeln(sum_log_field);
	}
}

void TMarkovField::simulate([[maybe_unused]] TDataModel &data_model) {
	throw std::runtime_error("Not implemented yet");
}

const TFieldStorage &TMarkovField::get_Y_matrix() const { return _Y; }

void TMarkovField::_calculate_log_prob_field(
    const IndexArray &index_in_leaves_space,
    std::array<coretools::TSumLogProbability, 2> &sum_log) const {
	for (size_t dim = 0; dim < _trees.size(); ++dim) {
		// get relevant clique
		const auto &tree = _trees[dim].get();

		// translate index in leaves to the index in tree
		const size_t index_in_tree =
		    tree->get_node_index_from_leaf_index(index_in_leaves_space[dim]);

		const size_t clique_index = tree->clique_index(index_in_leaves_space);

		tree->calculate_log_prob_parent_to_node(
		    index_in_tree, clique_index, tree->get_binned_branch_length(index_in_tree), sum_log);
	}
}

#ifdef USE_LOTUS
void TMarkovField::_calc_lotus_LL(const IndexArray &index_in_leaves_space,
                                  std::array<double, 2> &prob, const TDataModel &data_model) {
	data_model.get_lotus().calculate_LL_update_Y(index_in_leaves_space, prob);
}
#endif

#ifdef USE_SIMPLE_ERROR_MODEL
void TMarkovField::_calc_simple_error_model_LL(const IndexArray &multidim_index,
                                               std::array<double, 2> &prob,
                                               const TDataModel &data_model) {
	data_model.get_simple_error_model().probabilities_for_Y_update(multidim_index, prob);
}

bool TMarkovField::_simple_error_model_disagrees(const IndexArray &multidim_index, bool new_state,
                                                 const TDataModel &data_model) {
	return data_model.get_simple_error_model().disagrees_with(multidim_index, new_state);
}
#endif

double TMarkovField::_calculate_complete_joint_density() {

	// we can initialize the sum_log_field for the joint probability of the Markov random field

	// Easy case: Y
	double sum_log_field = coretools::containerSum(_complete_log_density);

	// now we loop over all Z to get the joint probability
	for (const auto &tree : _trees) { sum_log_field += tree->get_complete_joint_density(); }

	return sum_log_field;
};

void TDataSweepAccumulator::commit(TDataModel &data_model) {
#ifdef USE_LOTUS
	double sum_new_LL = 0.0;
	for (auto &i : _lotus_LL) {
		// loop over all LL (stored per thread) and sum
		sum_new_LL += i.getSum();
	}
	data_model.get_lotus().update_cur_LL(sum_new_LL);
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
	size_t total_disagree = 0;
	for (const auto &i : _n_disagree) { total_disagree += i; }
	// The sweep visits every cell of Y exactly once, so this is the complete disagreement count.
	data_model.get_simple_error_model().set_n_disagree(total_disagree);
#endif
}
