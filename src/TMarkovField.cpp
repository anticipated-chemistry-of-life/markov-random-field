//
// Created by madleina on 10.01.25.
//

#include "TMarkovField.h"
#include "smart_binary_search.h"




void TMarkovField::initialize() {
	// TODO: initialize trees

	// read K (sheet size for updating Y)
	_K = coretools::instances::parameters().get("K", 100);

	// number of outer loops = the number of times we need to repeat K in order to parse the first dimension completely
	_num_outer_loops = std::ceil((double)_trees[0].size() / (double)_K);

	// set number of leaves per dimension (omit the last dimension)
	_num_leaves_per_dim_except_last.resize(_trees.size() - 1);
	for (size_t i = 0; i < _trees.size() - 1; ++i) {
		_num_leaves_per_dim_except_last[i] = _trees[i].get_number_of_leaves();
	}

	// resize sheets: one per dimension except the last dimension
	_sheets.resize(_trees.size() - 1);
}

void TMarkovField::_update_sheets(bool first, const std::vector<size_t> &ix, const std::vector<size_t> &previous_ix) {
	for (size_t j = 0; j < ix.size(); ++j) {
		if (first || previous_ix[j] != ix[j]) {
			// first iteration or different dimension than before -> re-compute sheet
			_sheets[j].fill();
		}
	}
}

void TMarkovField::_fill_clique_along_last_dim(){
	_clique_last_dim.fill();
}

void TMarkovField::_update_Y(){
	for (size_t i = 0; i < _K; ++i){

	}
}

void TMarkovField::update_Y() {
	// loop over sheets in last dimension
	for (size_t k = 0; k < _num_outer_loops; ++k) {
		const size_t start_last_dim = k * _K; // 0, _K, 2*_K, ...

		// loop over all dimensions except last (linearized)
		size_t num_inner_loops = coretools::containerProduct(_num_leaves_per_dim_except_last);
		std::vector<size_t> previous_ix;
		for (size_t i = 0; i < num_inner_loops; ++i) {
			// get multi-dimensional index from linear coordinate (index in all dimensions except the last dimension)
			const std::vector<size_t> ix = coretools::getSubscripts(i, _num_leaves_per_dim_except_last);
			// update sheet(s), if necessary
			const size_t end_last_dim = std::min(start_last_dim + _K, _trees.back().get_number_of_leaves());
			_update_sheets(i == 0, ix, previous_ix);

			// fill clique along last dimension
			_fill_clique_last_dim();

			// loop along last dimension (start until start + K - 1) for updating
			for (size_t j = start_last_dim; j < end_last_dim; ++j) {
				// update Y for which we have all data
				_update_Y();
			}
		}
	}
}