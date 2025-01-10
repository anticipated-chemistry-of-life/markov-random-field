//
// Created by madleina on 10.01.25.
//

#ifndef ACOL_TMARKOVFIELD_H
#define ACOL_TMARKOVFIELD_H

#include "TTree.h"

//-----------------------------------
// TMarkovField
//-----------------------------------

class TMarkovField {
private:
	std::vector<TTree> _trees;
	TStorageYVector _Y;

	// for updating Y
	size_t _K;
	size_t _num_outer_loops;
	std::vector<size_t> _num_leaves_per_dim_except_last;
	std::vector<TSheet> _sheets;
	TCurrentState _clique_last_dim;

public:
	TMarkovField()  = default;
	~TMarkovField() = default;

	// updates
	void update_Y();
	void update_Z();
};

#endif // ACOL_TMARKOVFIELD_H
