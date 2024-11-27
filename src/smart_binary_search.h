#pragma once

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
binary_search(const T &vec, size_t linear_index_in_container_space,
              std::variant<size_t, typename T::const_iterator> lower_interval,
              std::variant<size_t, typename T::const_iterator> upper_interval) {

	// Determine the upper bound based on the type of upper_interval
	auto begin_it = std::holds_alternative<size_t>(lower_interval)
	                    ? vec.begin() + std::get<size_t>(lower_interval)
	                    : std::get<typename T::const_iterator>(lower_interval);
	auto end_it   = std::holds_alternative<size_t>(upper_interval) ? vec.begin() + std::get<size_t>(upper_interval)
	                                                               : std::get<typename T::const_iterator>(upper_interval);

	if (end_it > vec.end()) { end_it = vec.end(); }

	// lower_bound return the first element that is not less than the value
	auto it = std::lower_bound(begin_it, end_it, linear_index_in_container_space);

	// if our coordinate is bigger than the biggest element in the vector
	// we say that we haven't found our element and that if we want to
	// insert it, we should insert it at the end of the vector
	if (it == vec.end()) { return {false, vec.size(), vec.size(), true}; }

	// else our coordinate is in the range of the coordinates in the vector
	// meaning that if we haven't found it, we will insert it at that position
	// to keep the vector sorted
	auto distance = std::distance(vec.begin(), it);
	if (it->get_linear_index_in_container_space() != linear_index_in_container_space) {
		return {false, distance, it->get_linear_index_in_container_space(), false};
	}

	// if we found the coordinate we return the index and true
	return {true, distance, it->get_linear_index_in_container_space(), (unsigned long)distance == vec.size() - 1};
};

template<typename Container>
std::tuple<std::vector<bool>, std::vector<bool>, std::vector<size_t>>
fill_current_state_easy(const Container &container,
                        const std::vector<size_t> &multidim_index_of_first_in_container_space,
                        size_t n_nodes_in_clique_of_container) {

	// NOTE : This is valid only when the dimension we are in is the last dimension. This allows us to increment the
	// index in the Y vector by 1 and get the next element in the Y vector.
	std::vector<bool> current_state(n_nodes_in_clique_of_container, false);
	std::vector<bool> exists_in_container(n_nodes_in_clique_of_container, false);
	std::vector<size_t> index_in_TStorageVector(n_nodes_in_clique_of_container, container.size());
	auto start_linear_index = container.get_linear_index_in_container_space(multidim_index_of_first_in_container_space);

	auto [found, index_in_TStorage, linear_index_in_container_space, is_last_element] =
	    binary_search(container, start_linear_index, container.begin(), start_linear_index + 1);
	index_in_TStorageVector[0] = index_in_TStorage;
	if (found) {
		current_state[0]       = container.is_one(index_in_TStorage);
		exists_in_container[0] = true;
		if (is_last_element) { return std::make_tuple(current_state, exists_in_container, index_in_TStorageVector); }
		index_in_TStorage += 1;
		linear_index_in_container_space = container[index_in_TStorage].get_linear_index_in_container_space();
	}
	if (is_last_element) { return std::make_tuple(current_state, exists_in_container, index_in_TStorageVector); }

	for (size_t i = 1; i < n_nodes_in_clique_of_container; ++i) {
		auto linear_index_in_container_space_of_i = start_linear_index + i;
		if (linear_index_in_container_space_of_i < linear_index_in_container_space) {
			index_in_TStorageVector[i] = index_in_TStorage;
			continue;
		} else if (linear_index_in_container_space_of_i == linear_index_in_container_space) {
			current_state[i]           = container.is_one(index_in_TStorage);
			exists_in_container[i]     = true;
			index_in_TStorageVector[i] = index_in_TStorage;
			index_in_TStorage += 1;
			if (index_in_TStorage == container.size()) {
				return std::make_tuple(current_state, exists_in_container, index_in_TStorageVector);
			}
			linear_index_in_container_space = container[index_in_TStorage].get_linear_index_in_container_space();
		} else {
			UERROR("The linear index can't be bigger than the upper bound ! That means that there are more elements in "
			       "the container than in the total possible combinations of container !");
		}
	}
	return std::make_tuple(current_state, exists_in_container, index_in_TStorageVector);
}

template<typename Container>
std::tuple<std::vector<bool>, std::vector<bool>, std::vector<size_t>>
fill_current_state_hard(const Container &container, size_t n_nodes_in_clique_of_container,
                        const std::vector<size_t> &multidim_index_of_first_in_container_space, size_t increment,
                        size_t total_size_of_container) {
	std::vector<bool> current_state(n_nodes_in_clique_of_container, false);
	std::vector<bool> exists_in_container(n_nodes_in_clique_of_container, false);
	std::vector<size_t> index_in_TStorageVector(n_nodes_in_clique_of_container, container.size());
	auto linear_start_index = container.get_linear_index_in_container_space(multidim_index_of_first_in_container_space);

	auto [found, index_in_TStorage, linear_index_in_container_space, is_last_element] =
	    binary_search(container, linear_start_index, container.begin(), linear_start_index + 1);
	index_in_TStorageVector[0] = index_in_TStorage;
	if (found) {
		current_state[0]       = container[index_in_TStorage].is_one();
		exists_in_container[0] = true;
	}

	const double p                      = (double)container.size() / (double)total_size_of_container;
	const double increment_p            = (double)increment * p;
	const double two_standard_deviation = 2 * std::sqrt(increment_p * (1 - p));
	const size_t jump_right             = static_cast<size_t>(std::ceil(increment_p + two_standard_deviation));
	const size_t jump_left = static_cast<size_t>(std::max(0.0, std::floor(increment_p - two_standard_deviation)));

	for (size_t i = 1; i < n_nodes_in_clique_of_container; ++i) {
		auto linear_index_in_container_space_of_i = linear_start_index + i * increment;

		if (linear_index_in_container_space_of_i < linear_index_in_container_space) {
			index_in_TStorageVector[i] = index_in_TStorage;
			continue;
		}
		if (is_last_element) { return std::make_tuple(current_state, exists_in_container, index_in_TStorageVector); }

		// calculate the upper bound
		auto upper_bound = index_in_TStorage + jump_right;
		if (upper_bound >= container.size()) { upper_bound = container.size() - 1; }
		auto upper_linear_index_in_container_space = container[upper_bound].get_linear_index_in_container_space();

		// calculate the lower bound
		auto lower_bound                           = index_in_TStorage + jump_left;
		auto lower_linear_index_in_container_space = container[lower_bound].get_linear_index_in_container_space();

		if (linear_index_in_container_space_of_i == upper_linear_index_in_container_space) {
			index_in_TStorage               = upper_bound;
			linear_index_in_container_space = upper_linear_index_in_container_space;
			current_state[i]                = container[index_in_TStorage].is_one();
			exists_in_container[i]          = true;
			index_in_TStorageVector[i]      = index_in_TStorage;
			continue;
		}
		if (linear_index_in_container_space_of_i == lower_linear_index_in_container_space) {
			index_in_TStorage               = lower_bound;
			linear_index_in_container_space = lower_linear_index_in_container_space;
			current_state[i]                = container[index_in_TStorage].is_one();
			exists_in_container[i]          = true;
			index_in_TStorageVector[i]      = index_in_TStorage;
			continue;
		}

		if (linear_index_in_container_space_of_i > upper_linear_index_in_container_space) {

			std::tie(found, index_in_TStorage, linear_index_in_container_space, is_last_element) =
			    binary_search(container, linear_index_in_container_space_of_i, upper_bound, container.end());
			index_in_TStorageVector[i] = index_in_TStorage;
			if (found) {
				current_state[i]       = container[index_in_TStorage].is_one();
				exists_in_container[i] = true;
				continue;
			}
		} else if (linear_index_in_container_space_of_i > lower_linear_index_in_container_space &&
		           linear_index_in_container_space_of_i < upper_linear_index_in_container_space) {
			std::tie(found, index_in_TStorage, linear_index_in_container_space, is_last_element) =
			    binary_search(container, linear_index_in_container_space_of_i, lower_bound, upper_bound);
			index_in_TStorageVector[i] = index_in_TStorage;
			if (found) {
				current_state[i]       = container[index_in_TStorage].is_one();
				exists_in_container[i] = true;
				continue;
			}
		} else {
			std::tie(found, index_in_TStorage, linear_index_in_container_space, is_last_element) =
			    binary_search(container, linear_index_in_container_space_of_i, index_in_TStorage, lower_bound);
			index_in_TStorageVector[i] = index_in_TStorage;
			if (found) {
				current_state[i]       = container[index_in_TStorage].is_one();
				exists_in_container[i] = true;
				continue;
			}
		}
	}
	return std::make_tuple(current_state, exists_in_container, index_in_TStorageVector);
}

class TCurrentState {
private:
	std::vector<bool> _current_state_Y;
	std::vector<bool> _exists_in_Y;
	std::vector<size_t> _index_in_TStorageYVector;
	std::vector<bool> _current_state_Z;
	std::vector<bool> _exists_in_Z;
	std::vector<size_t> _index_in_TStorageZVector;
	const std::vector<size_t> &_multidim_index_of_first_in_container_space;
	size_t _increment;
	const TTree &_tree;

public:
	TCurrentState(const std::vector<size_t> &multidim_index_of_first_in_container_space, size_t increment,
	              const TStorageYVector &Y, const TStorageZVector &Z, const TTree &tree)
	    : _multidim_index_of_first_in_container_space(multidim_index_of_first_in_container_space),
	      _increment(increment), _tree(tree) {
		if (increment == 1) {
			std::tie(_current_state_Y, _exists_in_Y, _index_in_TStorageYVector) =
			    fill_current_state_easy(Y, _multidim_index_of_first_in_container_space, tree.get_number_of_leaves());
			std::tie(_current_state_Z, _exists_in_Z, _index_in_TStorageZVector) = fill_current_state_easy(
			    Z, _multidim_index_of_first_in_container_space, tree.get_number_of_internal_nodes());
		} else {
			std::tie(_current_state_Y, _exists_in_Y, _index_in_TStorageYVector) = fill_current_state_hard(
			    Y, tree.get_number_of_leaves(), _multidim_index_of_first_in_container_space, _increment, Y.size());
			std::tie(_current_state_Z, _exists_in_Z, _index_in_TStorageZVector) =
			    fill_current_state_hard(Z, tree.get_number_of_internal_nodes(),
			                            _multidim_index_of_first_in_container_space, _increment, Z.size());
		}
	}

	bool get(size_t index_in_tree) const {
		if (_tree.get_node(index_in_tree).isLeaf()) {
			return _current_state_Y[_tree.get_index_within_leaves(index_in_tree)];
		}
		return _current_state_Z[_tree.get_index_within_internal_nodes(index_in_tree)];
	}
	void set(size_t index_in_tree, bool value) {
		if (_tree.get_node(index_in_tree).isLeaf()) {
			_current_state_Y[_tree.get_index_within_leaves(index_in_tree)] = value;
		} else {
			_current_state_Z[_tree.get_index_within_internal_nodes(index_in_tree)] = value;
		}
	}

	size_t get_index_in_TStorageVector(size_t index_in_tree) const {
		if (_tree.get_node(index_in_tree).isLeaf()) {
			return _index_in_TStorageYVector[_tree.get_index_within_leaves(index_in_tree)];
		}
		return _index_in_TStorageZVector[_tree.get_index_within_internal_nodes(index_in_tree)];
	}

	bool exists_in_TStorageVector(size_t index_in_tree) const {
		if (_tree.get_node(index_in_tree).isLeaf()) {
			return _exists_in_Y[_tree.get_index_within_leaves(index_in_tree)];
		}
		return _exists_in_Z[_tree.get_index_within_internal_nodes(index_in_tree)];
	}
};
