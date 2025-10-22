//
// Created by madleina on 10.01.25.
//

#include "TMarkovField.h"
#include "TClique.h"
#include "TCurrentState.h"
#include "TLotus.h"
#include "TStorageYVector.h"
#include "TTree.h"
#include "Types.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/progressTools.h"
#include "coretools/algorithms.h"
#include "coretools/devtools.h"
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

TMarkovField::TMarkovField(size_t n_iterations, std::vector<std::unique_ptr<TTree>> &Trees, std::string _prefix)
    : _trees(Trees), _prefix(std::move(_prefix)), _clique_last_dim(*_trees.back().get(), 1) {
	using namespace coretools::instances;

	// read K (sheet size for updating Y)
	_K = parameters().get("K", 100000);

	// read: fix Y or Z?
	_fix_Y = !parameters().get("Y.update", true);
	if (_fix_Y) { logfile().list("Will fix Y during the MCMC."); }
	_fix_Z = !parameters().get("Z.update", true);
	if (_fix_Z) { logfile().list("Will fix Z during the MCMC."); }

	// number of outer loops = the number of times to repeat K such that all leaves of the last dimension are parsed
	_num_outer_loops = std::ceil((double)_trees.back()->get_number_of_leaves() / (double)_K);

	// set number of leaves per dimension (set the last dimension to one)
	_num_leaves_per_dim_except_last.resize(_trees.size(), 1);
	for (size_t i = 0; i < _trees.size() - 1; ++i) {
		_num_leaves_per_dim_except_last[i] = _trees[i]->get_number_of_leaves();
	}

	// create sheets: one per dimension except the last dimension
	_sheets.reserve(_trees.size() - 1);
	for (size_t i = 0; i < _trees.size() - 1; ++i) { _sheets.emplace_back(i, *_trees[i].get(), *_trees.back().get()); }

	// initialize Y
	auto num_leaves_per_dim   = _num_leaves_per_dim_except_last;
	num_leaves_per_dim.back() = _trees.back()->get_number_of_leaves();
	_Y.initialize(n_iterations, num_leaves_per_dim);
	if (parameters().exists("set_Y")) {
		std::string filename = parameters().get("set_Y", "acol_simulated_Y.txt");
		_read_Y_from_file(filename);
	}
}

bool TMarkovField::_need_to_update_sheet(size_t sheet_ix, const std::vector<size_t> &start_index_in_leaves_space,
                                         const std::vector<size_t> &previous_ix) const {
	for (size_t j = 0; j < _sheets.size(); ++j) { // loop over all sheets
		if (j == sheet_ix) { continue; }          // ignore current sheet
		if (start_index_in_leaves_space[j] != previous_ix[j]) { return true; }
	}
	return false;
}

void TMarkovField::_update_sheets(bool first, const std::vector<size_t> &start_index_in_leaves_space,
                                  const std::vector<size_t> &previous_ix, size_t K_cur_sheet) {
	for (size_t j = 0; j < _sheets.size(); ++j) {
		if (first || _need_to_update_sheet(j, start_index_in_leaves_space, previous_ix)) {
			// first iteration or different index than before -> re-compute sheet
			_sheets[j].fill(start_index_in_leaves_space, K_cur_sheet, _Y);
		}
	}
}

void TMarkovField::_fill_clique_along_last_dim(std::vector<size_t> start_index_in_leaves_space) {
	start_index_in_leaves_space.back() = 0; // set last dimension to zero
	// fill all Y along last dimension
	_clique_last_dim.fill_Y_along_last_dim(start_index_in_leaves_space, _trees.back()->get_number_of_leaves(), _Y);
	// fill all Z along last dimension
	_clique_last_dim.fill_Z_along_last_dim(start_index_in_leaves_space, _trees.back()->get_number_of_internal_nodes(),
	                                       _trees.back()->get_Z());
}

void TMarkovField::_calculate_log_prob_field(const std::vector<size_t> &index_in_leaves_space,
                                             std::array<coretools::TSumLogProbability, 2> &sum_log) const {
	for (size_t dim = 0; dim < _trees.size(); ++dim) {
		// get relevant clique
		const auto &clique = _trees[dim]->get_clique(index_in_leaves_space);

		// translate index in leaves to the index in tree
		const size_t index_in_tree = _trees[dim]->get_node_index_from_leaf_index(index_in_leaves_space[dim]);

		// get leaf index in tree of last dimension
		const size_t leaf_index_in_tree_of_last_dim = index_in_leaves_space.back();

		// calculate P(parent | node = 0) and P(parent | node = 1)
		// Note: leaves can never be roots -> they always have a parent (no need to bother with stationary)
		if (dim == _trees.size() - 1) { // last dim -> use _clique_last_dim
			clique.calculate_log_prob_parent_to_node(
			    index_in_tree, _trees[dim]->get_binned_branch_length(index_in_tree), _trees[dim].get(),
			    leaf_index_in_tree_of_last_dim, _clique_last_dim, sum_log);
		} else { // use sheet
			clique.calculate_log_prob_parent_to_node(
			    index_in_tree, _trees[dim]->get_binned_branch_length(index_in_tree), _trees[dim].get(),
			    leaf_index_in_tree_of_last_dim, _sheets[dim], sum_log);
		}
	}
}

void TMarkovField::_update_counter_1_cliques(bool new_state, bool old_state,
                                             const std::vector<size_t> &index_in_leaves_space) {
	// update counter of leaves with value 1 for all dimensions except the last one
	// reason: we parallelize over the last dimension -> can not update the counter there, as this would result in race
	// condition
	for (size_t dim = 0; dim < _trees.size() - 1; ++dim) {
		_trees[dim]->get_clique(index_in_leaves_space).update_counter_leaves_state_1(new_state, old_state);
	}
}

void TMarkovField::_update_cur_LL_lotus(TLotus &lotus, std::vector<coretools::TSumLogProbability> &new_LL) {
	double sum_new_LL = 0.0;
	for (auto &i : new_LL) {
		// loop over all LL (stored per thread) and sum
		sum_new_LL += i.getSum();
	}
	lotus.update_cur_LL(sum_new_LL);
}

int TMarkovField::_set_new_Y(bool new_state, const std::vector<size_t> &index_in_leaves_space,
                             std::vector<TStorageY> &linear_indices_in_Y_space_to_insert) {
	const size_t leaf_index_in_tree_of_last_dim = index_in_leaves_space.back();

	// get current state, exists and index in TStorageYVector
	// get it from _clique_last_dim, as this is re-computed for each update and is thus reliable with respect to indices
	// of where TStorageYVector is a one
	const auto [cur_state, exists_in_TStorageYVector, index_in_TStorageYVector] =
	    _clique_last_dim.get_state_exist_ix_TStorageYVector(leaf_index_in_tree_of_last_dim);

	if (cur_state && !new_state) { // 1 -> 0: set Y to zero
		_Y.set_to_zero(index_in_TStorageYVector);
	} else if (!cur_state && new_state) { // 0 -> 1
		if (exists_in_TStorageYVector) {
			_Y.set_to_one(index_in_TStorageYVector);
		} else { // does not yet exist -> remember linear index in Y to insert later
			const size_t linear_index_in_Y_space = _Y.get_linear_index_in_Y_space(index_in_leaves_space);
			linear_indices_in_Y_space_to_insert.emplace_back(linear_index_in_Y_space);
		}
	}

	// update values of 1-valued leaves in clique
	_update_counter_1_cliques(new_state, cur_state, index_in_leaves_space);
	int diff_counter_1_in_last_dim = (int)new_state - (int)cur_state;

	// update value in _sheets
	for (size_t dim = 0; dim < _trees.size() - 1; ++dim) {
		// translate leaf index to node index
		const size_t node_index_in_tree_of_dim =
		    _trees[dim]->get_node_index_from_leaf_index(index_in_leaves_space[dim]);
		_sheets[dim].set(node_index_in_tree_of_dim, leaf_index_in_tree_of_last_dim, new_state);
	}
	// Note: no need to update in _clique_last_dim, as this will anyway be overwritten for next Y

	return diff_counter_1_in_last_dim;
}

void TMarkovField::_read_Y_from_file(const std::string &filename) {
	coretools::TInputFile file(filename, coretools::FileType::Header);
	if (file.numCols() != 5) { UERROR("Simulated Y is expected to have 5 columns, but has ", file.numCols(), " !"); }

	// read each line of the file
	for (; !file.empty(); file.popFront()) {
		auto linear_index_in_Y_space = file.get<uint64_t>(0);
		bool state                   = file.get<bool>(1);
		if (state) { _Y.insert_one(linear_index_in_Y_space); }
	}
}

void TMarkovField::update(TLotus &lotus, size_t iteration) {
	if (WRITE_JOINT_LOG_PROB_DENSITY && iteration == 0 && !_joint_density_file.isOpen()) {
		_joint_density_file.open(_prefix + "_simulated_joint_density.txt",
		                         {
		                             "joint_density",
		                         },
		                         "\t");
	}
	_update_all_Y<false>(lotus, iteration);
	if (iteration == 0 && !_z_initialized_from_children) {
		for (auto &tree : _trees) { tree->initialize_Z_from_children(_Y); }
		_z_initialized_from_children = true;
	}
	if (_fix_Z) {
		_update_all_Z<false, true>(iteration);
	} else {
		_update_all_Z<false, false>(iteration);
	}
	_Y.add_to_counter(iteration);
	// calculate joint density
	if (WRITE_JOINT_LOG_PROB_DENSITY && iteration % _Y.get_thinning_factor() == 0) {
		auto sum_log_field = _calculate_complete_joint_density();
		_joint_density_file.writeln(sum_log_field);
	}
}

void TMarkovField::_calc_lotus_LL(const std::vector<size_t> &index_in_leaves_space, size_t index_for_tmp_state,
                                  size_t leaf_index_last_dim, std::array<double, 2> &prob, TLotus &lotus) {
	const bool cur_state = _clique_last_dim.get_Y(leaf_index_last_dim);
	lotus.calculate_LL_update_Y(index_in_leaves_space, index_for_tmp_state, cur_state, prob);
}

void TMarkovField::_prepare_lotus_LL(const std::vector<size_t> &start_index_in_leaves_space, size_t K_cur_sheet,
                                     TLotus &lotus) {
	lotus.fill_tmp_state_along_last_dim(start_index_in_leaves_space, K_cur_sheet);
}

void TMarkovField::simulate(TLotus &lotus) {
	// For simulation we always draw from the prior. (top-down)
	// 1. Draw branch len -> draw mus
	// 2. For every tree draw the root from those mus and the we BFS sample all the internal nodes based on the state of
	// the parent. At the end also draw the state of the leaves. For each tree we know which leaves and which internal
	// node are 0 or 1.
	// 2. draw Z since its unique per tree
	// 3. Draw Y (this function)
	// 4. Draw Lotus.
	// 5. From that simulated Lotus can we go back and infer the prior params
	//
	for (size_t tree_index = 0; tree_index < _trees.size(); ++tree_index) {
		auto &tree = _trees[tree_index];
		tree->simulate_Z(tree_index);
	}

	_simulate_Y();

	// for iteration in 1->max_iteration, (max_iteration should be passed from CLI)
	// we use tree.update_Z(); and then
	// update Y where likelihood of data is always one so it doesn't matter.
	size_t max_iteration = get_num_iterations_simulation();

	// create the Markov field density file
	if (WRITE_JOINT_LOG_PROB_DENSITY) {
		_joint_density_file.open(_prefix + "_simulated_joint_density.txt",
		                         {
		                             "joint_density",
		                         },
		                         "\t");
	}

	std::string report = "Running an MCMC chain of " + coretools::str::toString(max_iteration) + " iterations";
	coretools::TProgressReporter prog(max_iteration, report);
	for (size_t iteration = 0; iteration < max_iteration; ++iteration) {

		_update_all_Y<true>(lotus, iteration);

		if (_fix_Z) {
			_update_all_Z<true, true>(iteration);
		} else {
			_update_all_Z<true, false>(iteration);
		}
		_Y.add_to_counter(iteration);

		// calculate joint density
		if (iteration % _Y.get_thinning_factor() == 0 && WRITE_JOINT_LOG_PROB_DENSITY) {
			auto sum_log_field = _calculate_complete_joint_density();
			_joint_density_file.writeln(sum_log_field);
		}
		prog.next();
	}
	prog.done();
	if (WRITE_Y) { _write_Y_to_file<true>(_prefix + "_simulated_Y.txt"); }
	if (WRITE_Z) {
		for (size_t tree_idx = 0; tree_idx < _trees.size(); ++tree_idx) {
			const auto &tree = _trees[tree_idx];
			tree->write_Z_to_file<true>(_prefix + "_simulated_Z_" + tree->get_tree_name() + ".txt", _trees, tree_idx);
		}
	}
}

void TMarkovField::_simulate_Y() {
	// to sample Y we need to know the state of the parent for each leaf that is represented in a Y entry.
	// We are going to iterate over all possible Y and sample givent the product of probabilities of the child given the
	// parent.
	// set number of leaves per dimension (set the last dimension to one)
	if (SIMULATION_NO_Y_INITIALIZATION) { return; }
	for (size_t linear_index_in_leaves_space = 0; linear_index_in_leaves_space < _Y.total_size_of_container_space();
	     ++linear_index_in_leaves_space) {
		std::vector<size_t> multidim_index_in_Y = _Y.get_multi_dimensional_index(linear_index_in_leaves_space);
		std::array<coretools::TSumLogProbability, 2> sum_log;
		for (size_t dim = 0; dim < _trees.size(); ++dim) {
			// get relevant clique
			const auto &clique = _trees[dim]->get_clique(multidim_index_in_Y);
			TCurrentState current_state(*_trees[dim].get(), clique.get_increment(), _trees[dim]->get_number_of_leaves(),
			                            _trees[dim]->get_number_of_internal_nodes());
			// translate index in leaves to the index in tree
			const size_t index_in_tree = _trees[dim]->get_node_index_from_leaf_index(multidim_index_in_Y[dim]);
			// calculate P(parent | node = 0) and P(parent | node = 1)
			// Note: leaves can never be roots -> they always have a parent (no need to bother with stationary)
			clique.calculate_log_prob_parent_to_node(index_in_tree,
			                                         _trees[dim]->get_binned_branch_length(index_in_tree),
			                                         _trees[dim].get(), 0, current_state, sum_log);
		}
		bool y_state = sample(sum_log);
		if (y_state) {
			_Y.insert_one(linear_index_in_leaves_space);
			for (auto &_tree : _trees) {
				// get relevant clique
				auto &clique = _tree->get_clique(multidim_index_in_Y);
				clique.update_counter_leaves_state_1(true, false);
			}
		}
	}
}

void TMarkovField::burninHasFinished() {
	_Y.reset_counts();
	_Y.remove_zeros();
}

void TMarkovField::MCMCHasFinished() {
	// write function to write the posterior state of Y to file
	_write_Y_to_file<false>(_prefix + "_Y_posterior.txt");
}

const TStorageYVector &TMarkovField::get_Y_vector() const { return _Y; }
const TStorageY &TMarkovField::get_Y(size_t index_in_TStorageYVector) const { return _Y[index_in_TStorageYVector]; }
size_t TMarkovField::size_Y() const { return _Y.size(); }

double TMarkovField::_calculate_complete_joint_density() {

	// we can initialize the sum_log_field for the joint probability of the Markov random field

	// Easy case: Y
	double sum_log_field = coretools::containerSum(_complete_log_density);

	// now we loop over all Z to get the joint probability
	for (const auto &tree : _trees) { sum_log_field += tree->get_complete_joint_density(); }

	return sum_log_field;
};
