
#ifndef TUPDATE_CURRENT_STATE_H
#define TUPDATE_CURRENT_STATE_H

#include "TStorageYVector.h"
#include "TStorageZVector.h"
#include "TTree.h"
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
	return {true, distance, it->get_coordinate(), (unsigned long)distance == vec.size() - 1};
};

template<typename Container>
std::vector<bool> fill_current_state_easy(const Container &container, const std::vector<size_t> &start_multidim_index,
                                          size_t n_nodes_in_container_clique) {

	// NOTE : This is valid only when the dimension we are in is the last dimension. This allows us to increment the
	// index in the Y vector by 1 and get the next element in the Y vector.
	std::vector<bool> current_state(n_nodes_in_container_clique, false);
	auto start_linear_index = container.get_linear_coordinate(start_multidim_index);

	auto [found, index_in_Y, index_in_full_Y, is_last_element] =
	    binary_search(container, start_linear_index, container.begin(), start_linear_index + 1);
	if (found) {
		current_state[0] = container[index_in_Y].is_one();
		if (is_last_element) { return current_state; }
		index_in_Y += 1;
		index_in_full_Y = container[index_in_Y].get_coordinate();
	}
	if (is_last_element) { return current_state; }

	for (size_t i = 1; i < n_nodes_in_container_clique; ++i) {
		auto linear_index = start_linear_index + i;
		if (linear_index < index_in_full_Y) {
			continue;
		} else if (linear_index == index_in_full_Y) {
			current_state[i] = true;
			index_in_Y += 1;
			if (index_in_Y == container.size()) { return current_state; }
			index_in_full_Y = container[index_in_Y].get_coordinate();
		} else {
			UERROR("The linear index can't be bigger than the upper bound ! That means that there are more elements in "
			       "the container than in the total possible combinations of container !");
		}
	}
	return current_state;
}

template<typename Container>
std::vector<bool> fill_current_state_hard(size_t n_nodes_in_container_clique, const Container &container,
                                          const std::vector<size_t> &multi_dim_start_index, size_t increment,
                                          size_t total_size_of_container) {
	std::vector<bool> current_state(n_nodes_in_container_clique, false);
	auto linear_start_index = container.get_linear_coordinate(multi_dim_start_index);

	auto [found, index_in_Y, index_in_full_Y, is_last_element] =
	    binary_search(container, linear_start_index, container.begin(), linear_start_index + 1);
	if (found) { current_state[0] = container[index_in_Y].is_one(); }

	const double p                      = (double)container.size() / (double)total_size_of_container;
	const double increment_p            = increment * p;
	const double two_standard_deviation = 2 * std::sqrt(increment_p * (1 - p));
	const size_t jump_right             = static_cast<size_t>(std::ceil(increment_p + two_standard_deviation));
	const size_t jump_left = static_cast<size_t>(std::max(0.0, std::floor(increment_p - two_standard_deviation)));

	for (size_t i = 1; i < n_nodes_in_container_clique; ++i) {
		auto curr_index_in_full_Y = linear_start_index + i * increment;

		if (curr_index_in_full_Y < index_in_full_Y) { continue; }
		if (is_last_element) { return current_state; }

		// calculate the upper bound
		auto upper_bound = index_in_Y + jump_right;
		if (upper_bound >= container.size()) { upper_bound = container.size() - 1; }
		auto upper_index_in_full_Y = container[upper_bound].get_coordinate();

		// calculate the lower bound
		auto lower_bound           = index_in_Y + jump_left;
		auto lower_index_in_full_Y = container[lower_bound].get_coordinate();

		if (curr_index_in_full_Y == upper_index_in_full_Y) {
			index_in_Y       = upper_bound;
			index_in_full_Y  = upper_index_in_full_Y;
			current_state[i] = true;
			continue;
		}
		if (curr_index_in_full_Y == lower_index_in_full_Y) {
			index_in_Y       = lower_bound;
			index_in_full_Y  = lower_index_in_full_Y;
			current_state[i] = true;
			continue;
		}

		if (curr_index_in_full_Y > upper_index_in_full_Y) {

			std::tie(found, index_in_Y, index_in_full_Y, is_last_element) =
			    binary_search(container, curr_index_in_full_Y, upper_bound, container.end());
			if (found) {
				current_state[i] = container[index_in_Y].is_one();
				continue;
			}
		} else if (curr_index_in_full_Y > lower_index_in_full_Y && curr_index_in_full_Y < upper_index_in_full_Y) {
			std::tie(found, index_in_Y, index_in_full_Y, is_last_element) =
			    binary_search(container, curr_index_in_full_Y, lower_bound, upper_bound);
			if (found) {
				current_state[i] = container[index_in_Y].is_one();
				continue;
			}
		} else {
			std::tie(found, index_in_Y, index_in_full_Y, is_last_element) =
			    binary_search(container, curr_index_in_full_Y, index_in_Y, lower_bound);
			if (found) {
				current_state[i] = container[index_in_Y].is_one();
				continue;
			}
		}
	}
	return current_state;
}

class TCurrentState {
private:
	std::vector<bool> _current_state_Y;
	std::vector<bool> _current_state_Z;
	const std::vector<size_t> _start_index;
	size_t _increment;
	const TTree &_tree;

public:
	TCurrentState(const std::vector<size_t> &start_index, size_t increment, const TStorageYVector &Y,
	              const TStorageZVector &Z, const TTree &tree)
	    : _start_index(start_index), _increment(increment), _tree(tree) {
		if (increment == 1) {
			_current_state_Y = fill_current_state_easy(Y, _start_index, tree.get_number_of_leaves());
			_current_state_Z = fill_current_state_easy(Z, _start_index, tree.get_number_of_internal_nodes());
		} else {
			_current_state_Y =
			    fill_current_state_hard(tree.get_number_of_leaves(), Y, _start_index, _increment, Y.size());
			_current_state_Z =
			    fill_current_state_hard(tree.get_number_of_internal_nodes(), Z, _start_index, _increment, Z.size());
		}
	}

	bool get(size_t index_in_clique) const {
		if (_tree.get_node(index_in_clique).isLeaf()) {
			return _current_state_Y[_tree.get_index_within_leaves(index_in_clique)];
		}
		return _current_state_Z[_tree.get_index_within_internal_nodes(index_in_clique)];
	}
	void set(size_t index_in_clique, bool value) {
		if (_tree.get_node(index_in_clique).isLeaf()) {
			_current_state_Y[_tree.get_index_within_leaves(index_in_clique)] = value;
		} else {
			_current_state_Z[_tree.get_index_within_internal_nodes(index_in_clique)] = value;
		}
	}
};

#endif // TUPDATE_CURRENT_STATE_H
