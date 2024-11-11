
#ifndef TUPDATE_CURRENT_STATE_H
#define TUPDATE_CURRENT_STATE_H
#include "TStorageY.h"
#include "TStorageYVector.h"
#include "coretools/Main/TError.h"
#include <cmath>
#include <cstddef>
#include <tuple>
#include <variant>
#include <vector>

template<typename T, typename Underlying = typename T::value_type>
inline std::tuple<bool, size_t, Underlying, bool>
binary_search(const T &vec, size_t linear_index, std::variant<size_t, typename T::const_iterator> lower_interval,
              std::variant<size_t, typename T::const_iterator> upper_interval) {

	// Determine the upper bound based on the type of upper_interval
	auto begin_it = std::holds_alternative<size_t>(lower_interval)
	                    ? vec.begin() + std::get<size_t>(lower_interval)
	                    : std::get<typename T::const_iterator>(lower_interval);
	auto end_it   = std::holds_alternative<size_t>(upper_interval) ? vec.begin() + std::get<size_t>(upper_interval)
	                                                               : std::get<typename T::const_iterator>(upper_interval);

	if (end_it > vec.end()) { end_it = vec.end(); }

	// lower_bound return the first element that is not less than the value
	auto it = std::lower_bound(begin_it, end_it, linear_index);

	// if our coordinate is bigger than the biggest element in the vector
	// we say that we haven't found our element and that if we want to
	// insert it, we should insert it at the end of the vector
	if (it == vec.end()) { return {false, vec.size(), vec.size(), true}; }

	// else our coordinate is in the range of the coordinates in the vector
	// meaning that if we haven't found it, we will insert it at that position
	// to keep the vector sorted
	auto distance = std::distance(vec.begin(), it);
	if (it->get_coordinate() != linear_index) { return {false, distance, it->get_coordinate(), false}; }

	// if we found the coordinate we return the index and true
	return {true, distance, it->get_coordinate(), distance == vec.size() - 1};
};

inline std::vector<bool> fill_current_state_easy(const TStorageYVector &Y,
                                                 const std::vector<size_t> &start_multidim_index,
                                                 size_t n_nodes_in_Y_clique) {

	// NOTE : This is valid only when the dimension we are in is the last dimension. This allows us to increment the
	// index in the Y vector by 1 and get the next element in the Y vector.
	std::vector<bool> current_state(n_nodes_in_Y_clique, false);
	auto start_linear_index = Y.get_linear_coordinate(start_multidim_index);

	auto [found, index_in_Y, index_in_full_Y, is_last_element] =
	    binary_search(Y, start_linear_index, Y.begin(), start_linear_index + 1);
	if (found) {
		current_state[0] = true;
		if (is_last_element) { return current_state; }
		index_in_Y += 1;
		index_in_full_Y = Y[index_in_Y].get_coordinate();
	}
	if (is_last_element) { return current_state; }

	for (size_t i = 1; i < n_nodes_in_Y_clique; ++i) {
		auto linear_index = start_linear_index + i;
		if (linear_index < index_in_full_Y) {
			continue;
		} else if (linear_index == index_in_full_Y) {
			current_state[i] = true;
			index_in_Y += 1;
			if (index_in_Y == Y.size()) { return current_state; }
			index_in_full_Y = Y[index_in_Y].get_coordinate();
		} else {
			UERROR("The linear index can't be bigger than the upper bound ! That means that there are more elements in "
			       "Y than in the total possible combinations of Y !");
		}
	}
	return current_state;
}

inline std::tuple<std::vector<bool>, int, int, int>
fill_current_state_hard(size_t n_nodes_in_Y_clique, const TStorageYVector &Y,
                        const std::vector<size_t> &multi_dim_start_index, size_t increment, size_t total_size_Y) {
	std::vector<bool> current_state(n_nodes_in_Y_clique, false);
	auto linear_start_index = Y.get_linear_coordinate(multi_dim_start_index);

	auto [found, index_in_Y, index_in_full_Y, is_last_element] =
	    binary_search(Y, linear_start_index, Y.begin(), linear_start_index + 1);
	if (found) { current_state[0] = true; }

	const double p                      = (double)Y.size() / (double)total_size_Y;
	const double increment_p            = increment * p;
	const double two_standard_deviation = 2 * std::sqrt(increment_p * (1 - p));
	const size_t jump_right             = static_cast<size_t>(std::ceil(increment_p + two_standard_deviation));
	const size_t jump_left = static_cast<size_t>(std::max(0.0, std::floor(increment_p - two_standard_deviation)));

	int lower_chunk  = 0;
	int upper_chunk  = 0;
	int middle_chunk = 0;
	for (size_t i = 1; i < n_nodes_in_Y_clique; ++i) {
		auto curr_index_in_full_Y = linear_start_index + i * increment;

		if (curr_index_in_full_Y < index_in_full_Y) { continue; }
		if (is_last_element) { return {current_state, lower_chunk, middle_chunk, upper_chunk}; }

		auto upper_bound = index_in_Y + jump_right;
		if (upper_bound >= Y.size()) { upper_bound = Y.size() - 1; }
		auto upper_index_in_full_Y = Y[upper_bound].get_coordinate();

		if (curr_index_in_full_Y == upper_index_in_full_Y) {
			found            = true;
			index_in_full_Y  = upper_index_in_full_Y;
			current_state[i] = true;
			continue;
		}
		if (curr_index_in_full_Y > upper_index_in_full_Y) {
			upper_chunk++;
			auto [found, index_in_Y, index_in_full_Y, is_last_element] =
			    binary_search(Y, curr_index_in_full_Y, upper_bound, Y.end());
			if (found) {
				current_state[i] = true;
				continue;
			}
		}

		auto lower_bound = index_in_Y + jump_left;

		auto lower_index_in_full_Y = Y[lower_bound].get_coordinate();
		if (curr_index_in_full_Y == lower_index_in_full_Y) {
			found            = true;
			index_in_full_Y  = lower_index_in_full_Y;
			current_state[i] = true;
			continue;
		}

		if (curr_index_in_full_Y > lower_index_in_full_Y) {
			middle_chunk++;
			auto [found, index_in_Y, index_in_full_Y, is_last_element] =
			    binary_search(Y, curr_index_in_full_Y, lower_bound, upper_bound);
			if (found) {
				current_state[i] = true;
				continue;
			}
		} else {
			lower_chunk++;
			auto [found, index_in_Y, index_in_full_Y, is_last_element] =
			    binary_search(Y, curr_index_in_full_Y, Y.begin(), lower_bound);
			if (found) {
				current_state[i] = true;
				continue;
			}
		}
	}
	return {current_state, lower_chunk, middle_chunk, upper_chunk};
}

#endif // TUPDATE_CURRENT_STATE_H
