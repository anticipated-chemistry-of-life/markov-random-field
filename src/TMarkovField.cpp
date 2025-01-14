//
// Created by madleina on 10.01.25.
//

#include "TMarkovField.h"
#include "TClique.h"
#include "smart_binary_search.h"

TMarkovField::TMarkovField()
    : _number_of_threads(coretools::getNumThreads()), _trees(_make_trees()), _clique_last_dim(_trees.back(), 1) {

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
	trees.emplace_back(0, _number_of_threads, filename_tree_molecules, "molecules");
	// middle trees: all others (e.g. tissues)
	for (size_t i = 1; i < num_trees - 1; ++i) {
		trees.emplace_back(i, _number_of_threads, filenames_tree_others[i - 1], filenames_tree_others[i - 1]);
	}
	// last tree: species
	trees.emplace_back(num_trees - 1, _number_of_threads, filename_tree_species, "species");

	return trees;
}

void TMarkovField::_update_sheets(bool first, const std::vector<size_t> &start_index_in_leaves_space,
                                  const std::vector<size_t> &previous_ix, size_t K_cur_sheet) {
	for (size_t j = 0; j < _sheets.size(); ++j) {
		if (first || previous_ix[j] != start_index_in_leaves_space[j]) {
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
		if (dim == _trees.size() - 1) { // last dim -> use _clique_last_dim
			clique.calculate_log_prob_parent_to_node(index_in_tree, _trees[dim], leaf_index_in_tree_of_last_dim,
			                                         _clique_last_dim, sum_log);
		} else {
			clique.calculate_log_prob_parent_to_node(index_in_tree, _trees[dim], leaf_index_in_tree_of_last_dim,
			                                         _sheets[dim], sum_log);
		}
	}
}

void TMarkovField::_update_Y(const std::vector<size_t> &index_in_leaves_space) {
	// prepare log probabilities for the two possible states
	std::array<coretools::TSumLogProbability, 2> sum_log;

	// calculate probabilities in random Markov field
	_calculate_log_prob_field(index_in_leaves_space, sum_log);

	// calculate likelihood (lotus)
	// ...

	// sample state
	bool new_state = sample(sum_log);
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
#pragma omp parallel for num_threads(this->_number_of_threads) schedule(static)
			for (size_t j = start_ix_in_leaves_last_dim; j < end_ix_in_leaves_last_dim; ++j) {
				start_index_in_leaves_space.back() = j;
				_update_Y(start_index_in_leaves_space);
			}

			previous_ix = start_index_in_leaves_space;
		}
	}
}