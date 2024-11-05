#pragma once

#include "TStorageY.h"
#include "TStorageYVector.h"
#include "TTree.h"
#include "coretools/Main/TError.h"
#include "coretools/algorithms.h"
#include <cmath>
#include <variant>
#include <vector>

template<typename T, typename Underlying = typename T::value_type>
std::tuple<bool, size_t, Underlying> binary_search(const T &vec, const std::vector<size_t> &multi_dim_index,
                                                   std::variant<size_t, typename T::const_iterator> lower_interval,
                                                   std::variant<size_t, typename T::const_iterator> upper_interval) {
	auto coordinate = vec.get_linear_coordinate(multi_dim_index);

	// Determine the upper bound based on the type of upper_interval
	auto begin_it = std::holds_alternative<size_t>(lower_interval)
	                    ? vec.begin() + std::get<size_t>(lower_interval)
	                    : std::get<typename T::const_iterator>(lower_interval);
	auto end_it   = std::holds_alternative<size_t>(upper_interval) ? vec.begin() + std::get<size_t>(upper_interval)
	                                                               : std::get<typename T::const_iterator>(upper_interval);

	// lower_bound return the first element that is not less than the value
	auto it = std::lower_bound(begin_it, end_it, coordinate);

	// if our coordinate is bigger than the biggest element in the vector
	// we say that we haven't found our element and that if we want to
	// insert it, we should insert it at the end of the vector
	if (it == vec.end()) { return {false, vec.size(), vec.size()}; }

	// else our coordinate is in the range of the coordinates in the vector
	// meaning that if we haven't found it, we will insert it at that position
	// to keep the vector sorted
	if (it->get_coordinate() != coordinate) { return {false, vec.size(), it->get_coordinate()}; }

	// if we found the coordinate we return the index and true
	return {true, std::distance(vec.begin(), it), it->get_coordinate()};
};

std::vector<bool> fill_current_state_easy(const TStorageYVector &Y, const std::vector<size_t> &start_multidim_index,
                                          size_t n_nodes_in_Y_clique) {

	// NOTE : This is valid only when the dimension we are in is the last dimension. This allows us to increment the
	// index in the Y vector by 1 and get the next element in the Y vector.
	std::vector<bool> current_state(n_nodes_in_Y_clique, false);

	auto [found, index_in_Y, index_in_full_Y] =
	    binary_search(Y, start_multidim_index, Y.begin(), Y.get_linear_coordinate(start_multidim_index) + 1);
	if (found) {
		current_state[0] = true;
		index_in_Y += 1;
		index_in_full_Y = Y[index_in_Y].get_coordinate();
	}
	auto start_index_linear = Y.get_linear_coordinate(start_multidim_index);
	for (size_t i = 1; i < n_nodes_in_Y_clique; ++i) {
		auto linear_index = start_index_linear + i;
		if (linear_index < index_in_full_Y) {
			continue;
		} else if (linear_index == index_in_full_Y) {
			current_state[i] = true;
			index_in_Y += 1;
		} else {
			index_in_full_Y = Y[index_in_Y].get_coordinate();
			UERROR("The linear index can't be bigger than the upper bound ! That means that there are more elements in "
			       "Y than in the total possible combinations of Y !");
		}
	}
	return current_state;
}

std::vector<bool> fill_current_state_hard(size_t n_nodes_in_Y_clique, const TStorageYVector &Y,
                                          const std::vector<size_t> &number_of_leaves_in_each_dimension,
                                          const std::vector<size_t> &linear_indices_of_clique) {
	if (n_nodes_in_Y_clique != linear_indices_of_clique.size()) {
		UERROR("The number of nodes in the clique is different from the number of linear indices of the clique !");
	}
	std::vector<bool> current_state(n_nodes_in_Y_clique, false);
	size_t lower_bound = 0;
	size_t upper_bound = linear_indices_of_clique[0];
	auto [found, index_in_Y, index_in_full_Y] =
	    binary_search(Y, linear_indices_of_clique, Y.begin(), linear_indices_of_clique[0] + 1);
	if (found) {
		current_state[0] = true;
		index_in_Y += 1;
		index_in_full_Y = Y[index_in_Y].get_coordinate();
		lower_bound     = upper_bound;
	}
	for (size_t i = 1; i < n_nodes_in_Y_clique; ++i) {
		if (linear_indices_of_clique[i] < index_in_full_Y) {
			continue;
		} else if (linear_indices_of_clique[i] == index_in_full_Y) {
			current_state[i] = true;
			index_in_Y += 1;
			index_in_full_Y = Y[index_in_Y].get_coordinate();
			lower_bound     = linear_indices_of_clique[i];
		}

		return current_state;
	}
