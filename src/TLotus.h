#pragma once

#include "TStorageYVector.h"
#include "TTree.h"
#include <cstddef>
#include <string>
#include <vector>

// For this class, we need to enforce that the first tree will always be the one of the species and
// the second tree will always be the one of the molecules. The rest of the trees, we won't care.
class TLotus {
private:
	// trees should be a const ref because we don't want to change the trees and don't want to copy them
	const std::vector<TTree> &_trees;
	TStorageYVector _L_sm;
	TStorageYVector _x_sm;
	std::vector<size_t> _species_counter;
	std::vector<size_t> _molecules_counter;

	// private functions
	static std::vector<size_t> _get_dimensions_Lotus_space(const std::vector<TTree> &trees) {
		return {trees[0].get_number_of_leaves(), trees[1].get_number_of_leaves()};
	}

	/// To construct the _data_X_of_Lotus we will need to have as input a vector of TStorageYVector but this time with
	/// as many dimensions as there are trees.
	void _initialize_x_sm(const TStorageYVector &Y);

	double _calculate_probability_of_L_sm(bool x_sm, bool L_sm, size_t linear_index_in_L_space) const;

public:
	// we add the trees when we construct the object.
	explicit TLotus(const std::vector<TTree> &trees) : _trees(trees) {
		if (trees[0].get_tree_name() != "species") {
			UERROR("The species tree was not found in the provided trees.");
		} else if (trees[1].get_tree_name() != "molecules") {
			UERROR("The molecules tree was not found in the provided trees.");
		}

		_L_sm.initialize(0, _get_dimensions_Lotus_space(trees));
		_x_sm.initialize(0, _get_dimensions_Lotus_space(trees));
	}

	// default destructor
	~TLotus() = default;

	void load_from_file(const std::string &filename);

	[[nodiscard]] const TStorageYVector &get_TStorageYVector() const { return _L_sm; }

	double calculate_research_effort(size_t linear_index_in_L_space) const;

	double calculate_log_likelihood_of_L() const;
};
