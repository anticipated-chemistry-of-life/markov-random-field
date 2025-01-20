//
// Created by madleina on 10.01.25.
//

#ifndef ACOL_TMARKOVFIELD_H
#define ACOL_TMARKOVFIELD_H

#include "TCurrentState.h"
#include "TTree.h"

//-----------------------------------
// TMarkovField
//-----------------------------------

class TMarkovField {
private:
	// trees and Y
	std::vector<TTree> _trees;
	TStorageYVector _Y;

	// stuff for updating Y
	size_t _K;
	size_t _num_outer_loops;
	std::vector<size_t> _num_leaves_per_dim_except_last;
	std::vector<TSheet> _sheets;
	TCurrentState _clique_last_dim;

	// functions for initializing
	std::vector<TTree> _make_trees();

	// functions for updating Y
	void _update_sheets(bool first, const std::vector<size_t> &start_index_in_leaves_space,
	                    const std::vector<size_t> &previous_ix, size_t K_cur_sheet);
	void _fill_clique_along_last_dim(const std::vector<size_t> &start_index_in_leaves_space);
	void _calculate_log_prob_field(const std::vector<size_t> &index_in_leaves_space,
	                               std::array<coretools::TSumLogProbability, 2> &sum_log) const;
	void _update_Y(std::vector<size_t> index_in_leaves_space, size_t leaf_index,
	               std::vector<TStorageY> &linear_indices_in_Y_space_to_insert);
	bool _need_to_update_sheet(size_t sheet_ix, const std::vector<size_t> &start_index_in_leaves_space,
	                           const std::vector<size_t> &previous_ix);
	void _set_new_Y(bool new_state, const std::vector<size_t> &index_in_leaves_space,
	                std::vector<TStorageY> &linear_indices_in_Y_space_to_insert);

public:
	TMarkovField();
	~TMarkovField() = default;

	// updates
	void update_Y();
	void update_Z();
};

#endif // ACOL_TMARKOVFIELD_H
