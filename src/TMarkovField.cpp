//
// Created by madleina on 10.01.25.
//

#include "TMarkovField.h"
#include "TClique.h"
#include "smart_binary_search.h"

TMarkovField::TMarkovField(size_t n_iterations) : _trees(_make_trees()), _clique_last_dim(_trees.back(), 1) {

	// read K (sheet size for updating Y)
	_K = coretools::instances::parameters().get("K", 100);

	// number of outer loops = the number of times to repeat K such that all leaves of the last dimension are parsed
	_num_outer_loops = std::ceil((double)_trees.back().get_number_of_leaves() / (double)_K);

	// set number of leaves per dimension (set the last dimension to one)
	_num_leaves_per_dim_except_last.resize(_trees.size(), 1);
	for (size_t i = 0; i < _trees.size() - 1; ++i) {
		_num_leaves_per_dim_except_last[i] = _trees[i].get_number_of_leaves();
	}

	// create sheets: one per dimension except the last dimension
	_sheets.reserve(_trees.size() - 1);
	for (size_t i = 0; i < _trees.size() - 1; ++i) { _sheets.emplace_back(i, _trees[i], _trees.back()); }

	// initialize Y
	auto num_leaves_per_dim   = _num_leaves_per_dim_except_last;
	num_leaves_per_dim.back() = _trees.back().get_number_of_leaves();
	_Y.initialize(n_iterations, num_leaves_per_dim);
}

std::vector<TTree> TMarkovField::_make_trees() {
	using namespace coretools::instances;
	// read filenames
	std::string filename_tree_species   = parameters().get("tree_species");
	std::string filename_tree_molecules = parameters().get("tree_molecules");
	std::vector<std::string> filenames_tree_others;
	if (parameters().exists("tree_others")) { parameters().fill("tree_others", filenames_tree_others); }

	size_t num_trees = 2 + filenames_tree_others.size();
	std::vector<TTree> trees;
	trees.reserve(num_trees);

	// first tree: molecules
	trees.emplace_back(0, filename_tree_molecules, "molecules");
	// middle trees: all others (e.g. tissues)
	for (size_t i = 1; i < num_trees - 1; ++i) {
		std::string name = coretools::str::split(filenames_tree_others[i - 1], ':');
		if (name.empty()) {
			UERROR("Argument 'tree_others': Please provide a name for each other tree, separated by a : from the "
			       "filename (e.g. myTreeName:pathToFile)");
		}
		trees.emplace_back(i, filenames_tree_others[i - 1], name);
	}
	// last tree: species
	trees.emplace_back(num_trees - 1, filename_tree_species, "species");

	return trees;
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

void TMarkovField::_fill_clique_along_last_dim(const std::vector<size_t> &start_index_in_leaves_space) {
	// fill all Y along last dimension
	_clique_last_dim.fill_Y_along_last_dim(start_index_in_leaves_space, _trees.back().get_number_of_leaves(), _Y);
	// fill all Z along last dimension
	_clique_last_dim.fill_Z_along_last_dim(start_index_in_leaves_space, _trees.back().get_number_of_internal_nodes(),
	                                       _trees.back().get_Z());
}

void TMarkovField::_calculate_log_prob_field(const std::vector<size_t> &index_in_leaves_space,
                                             std::array<coretools::TSumLogProbability, 2> &sum_log) const {
	for (size_t dim = 0; dim < _trees.size(); ++dim) {
		// get relevant clique
		const auto &clique         = _trees[dim].get_clique(index_in_leaves_space);
		// translate index in leaves to the index in tree
		const size_t index_in_tree = _trees[dim].get_node_index_from_leaf_index(index_in_leaves_space[dim]);
		// get leaf index in tree of last dimension
		const size_t leaf_index_in_tree_of_last_dim = index_in_leaves_space.back();
		// calculate P(parent | node = 0) and P(parent | node = 1)
		// Note: leaves can never be roots -> they always have a parent (no need to bother with stationary)
		if (dim == _trees.size() - 1) { // last dim -> use _clique_last_dim
			clique.calculate_log_prob_parent_to_node(index_in_tree, _trees[dim], leaf_index_in_tree_of_last_dim,
			                                         _clique_last_dim, sum_log);
		} else { // use sheet
			clique.calculate_log_prob_parent_to_node(index_in_tree, _trees[dim], leaf_index_in_tree_of_last_dim,
			                                         _sheets[dim], sum_log);
		}
	}
}

int TMarkovField::_update_Y(std::vector<size_t> index_in_leaves_space, size_t leaf_index_last_dim,
                             std::vector<TStorageY> &linear_indices_in_Y_space_to_insert) {
	index_in_leaves_space.back() = leaf_index_last_dim;

	// prepare log probabilities for the two possible states
	std::array<coretools::TSumLogProbability, 2> sum_log;

	// calculate probabilities in Markov random field
	_calculate_log_prob_field(index_in_leaves_space, sum_log);

	// calculate log likelihood (lotus)
	// ...

	// calculate log likelihood (virtual mass spec)...

	// sample state
	bool new_state = sample(sum_log);

	// update Y accordingly
	return _set_new_Y(new_state, index_in_leaves_space, linear_indices_in_Y_space_to_insert);
}

void TMarkovField::_update_counter_1_cliques(bool new_state, bool old_state,
                                             const std::vector<size_t> &index_in_leaves_space) {
	// update counter of leaves with value 1 for all dimensions except the last one
	// reason: we parallelize over the last dimension -> can not update the counter there, as this would result in race
	// condition
	for (size_t dim = 0; dim < _trees.size() - 1; ++dim) {
		_trees[dim].get_clique(index_in_leaves_space).update_counter_leaves_state_1(new_state, old_state);
	}
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
		const size_t node_index_in_tree_of_dim = _trees[dim].get_node_index_from_leaf_index(index_in_leaves_space[dim]);
		_sheets[dim].set(node_index_in_tree_of_dim, leaf_index_in_tree_of_last_dim, new_state);
	}
	// Note: no need to update in _clique_last_dim, as this will anyway be overwritten for next Y

	return diff_counter_1_in_last_dim;
}

void TMarkovField::update_Y() {
	// loop over sheets in last dimension
	for (size_t k = 0; k < _num_outer_loops; ++k) {
		const size_t start_ix_in_leaves_last_dim = k * _K; // 0, _K, 2*_K, ...

		// loop over all dimensions except last (linearized)
		size_t num_inner_loops = coretools::containerProduct(_num_leaves_per_dim_except_last);
		std::vector<size_t> previous_ix;
		for (size_t i = 0; i < num_inner_loops; ++i) {
			// get multi-dimensional index from linear coordinate and set the start of the last dimension
			auto start_index_in_leaves_space   = coretools::getSubscripts(i, _num_leaves_per_dim_except_last);
			start_index_in_leaves_space.back() = start_ix_in_leaves_last_dim;
			// calculate size of current sheet (make sure not to overshoot)
			const size_t K_cur_sheet = std::min(_K, _trees.back().get_number_of_leaves() - start_ix_in_leaves_last_dim);
			// update sheet(s), if necessary
			_update_sheets(i == 0, start_index_in_leaves_space, previous_ix, K_cur_sheet);

			// fill clique along last dimension
			start_index_in_leaves_space.back() = 0; // start from the beginning
			_fill_clique_along_last_dim(start_index_in_leaves_space);

			// now loop along all leaves of the last dimension for updating (only K leaves for which we have everything)
			const size_t end_ix_in_leaves_last_dim = start_ix_in_leaves_last_dim + K_cur_sheet;
			std::vector<std::vector<TStorageY>> linear_indices_in_Y_space_to_insert(NUMBER_OF_THREADS);
			int diff_counter_1_in_last_dim = 0;
#pragma omp parallel for num_threads(NUMBER_OF_THREADS) reduction(+:diff_counter_1_in_last_dim)
			for (size_t j = start_ix_in_leaves_last_dim; j < end_ix_in_leaves_last_dim; ++j) {
				diff_counter_1_in_last_dim += _update_Y(start_index_in_leaves_space, j, linear_indices_in_Y_space_to_insert[omp_get_thread_num()]);
			}

			// insert new 1-valued indices into Y
			// Note: indices of where Y is one in _sheets is not accurate anymore, but we don't use them, so it's ok
			_Y.insert_in_Y(linear_indices_in_Y_space_to_insert);
			_trees.back()
			    .get_clique(start_index_in_leaves_space)
			    .update_counter_leaves_state_1(diff_counter_1_in_last_dim);

			previous_ix = start_index_in_leaves_space;
		}
	}
}

void TMarkovField::update_Z() {
	for (auto &_tree : _trees) { _tree.update_Z(_Y); }
}
