#pragma once

#include "coretools/Main/TError.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <tuple>
#include <vector>
class TTree;

// lower_idx and upper_idx are indices into vec (upper_idx is exclusive, like end()).
// Passing upper_idx > vec.size() is safe — it is clamped to vec.size().
template<typename Container>
inline std::tuple<bool, size_t, size_t, bool> binary_search(const Container &vec, size_t target,
                                                            size_t lower_idx, size_t upper_idx) {
	const auto begin_it = vec.begin() + lower_idx;
	const auto end_it   = vec.begin() + std::min(upper_idx, vec.size());

	auto it = std::lower_bound(begin_it, end_it, target);

	if (it == vec.end()) { return {false, vec.size(), vec.size(), true}; }

	const size_t pos = static_cast<size_t>(it - vec.begin());
	if (it->get_linear_index_in_container_space() != target) {
		return {false, pos, static_cast<size_t>(it->get_linear_index_in_container_space()), false};
	}
	return {true, pos, static_cast<size_t>(it->get_linear_index_in_container_space()),
	        pos == vec.size() - 1};
}

struct CurrentStateResult {
	std::vector<uint8_t> current_state;
	std::vector<uint8_t> exists_in_container;
	std::vector<size_t> index_in_TStorageVector;
};

template<typename Container>
CurrentStateResult fill_current_state_easy(const Container &container,
                                           const std::vector<size_t> &start_index_in_leaves_space,
                                           size_t n_nodes_in_clique_of_container) {

	// NOTE : This is valid only when the dimension we are in is the last dimension. This allows us
	// to increment the index in the Y vector by 1 and get the next element in the Y vector.
	std::vector<uint8_t> current_state(n_nodes_in_clique_of_container, false);
	std::vector<uint8_t> exists_in_container(n_nodes_in_clique_of_container, false);
	std::vector<size_t> index_in_TStorageVector(n_nodes_in_clique_of_container, container.size());
	size_t start_linear_index =
	    container.get_linear_index_in_container_space(start_index_in_leaves_space);

	auto [found, index_in_TStorage, linear_index_in_container_space, is_last_element] =
	    binary_search(container, start_linear_index, size_t{0}, start_linear_index + 1);
	index_in_TStorageVector[0] = index_in_TStorage;
	if (found) {
		current_state[0]       = container.is_one(index_in_TStorage);
		exists_in_container[0] = true;
		if (is_last_element) {
			return {current_state, exists_in_container, index_in_TStorageVector};
		}
		index_in_TStorage += 1;
		linear_index_in_container_space =
		    container[index_in_TStorage].get_linear_index_in_container_space();
	}
	if (is_last_element) { return {current_state, exists_in_container, index_in_TStorageVector}; }

	for (size_t i = 1; i < n_nodes_in_clique_of_container; ++i) {
		const auto linear_index_in_container_space_of_i = start_linear_index + i;
		if (linear_index_in_container_space_of_i < linear_index_in_container_space) {
			index_in_TStorageVector[i] = index_in_TStorage;
			continue;
		} else if (linear_index_in_container_space_of_i == linear_index_in_container_space) {
			current_state[i]           = container.is_one(index_in_TStorage);
			exists_in_container[i]     = true;
			index_in_TStorageVector[i] = index_in_TStorage;
			index_in_TStorage += 1;
			if (index_in_TStorage == container.size()) {
				return {current_state, exists_in_container, index_in_TStorageVector};
			}
			linear_index_in_container_space =
			    container[index_in_TStorage].get_linear_index_in_container_space();
		} else {
			throw coretools::TUserError(
			    "The linear index can't be bigger than the upper bound ! That means that there are "
			    "more elements in "
			    "the container than in the total possible combinations of container !");
		}
	}
	return {current_state, exists_in_container, index_in_TStorageVector};
}

template<typename Container>
CurrentStateResult
fill_current_state_hard_with_linear_window(const Container &container,
                                           size_t n_nodes_in_clique_of_container,
                                           const std::vector<size_t> &start_index_in_leaves_space,
                                           size_t increment, size_t total_size_of_container) {
	std::vector<uint8_t> current_state(n_nodes_in_clique_of_container, false);
	std::vector<uint8_t> exists_in_container(n_nodes_in_clique_of_container, false);
	std::vector<size_t> index_in_TStorageVector(n_nodes_in_clique_of_container, container.size());
	auto linear_start_index =
	    container.get_linear_index_in_container_space(start_index_in_leaves_space);

	auto [found, index_in_TStorage, linear_index_in_container_space, is_last_element] =
	    binary_search(container, linear_start_index, size_t{0}, linear_start_index + 1);
	index_in_TStorageVector[0] = index_in_TStorage;
	if (found) {
		current_state[0]       = container[index_in_TStorage].is_one();
		exists_in_container[0] = true;
	}

	const double p =
	    static_cast<double>(container.size()) / static_cast<double>(total_size_of_container);
	const double increment_p            = static_cast<double>(increment) * p;
	const double two_standard_deviation = 2 * std::sqrt(increment_p * (1 - p));
	const auto jump_right = static_cast<size_t>(std::ceil(increment_p + two_standard_deviation));
	const size_t jump_left =
	    static_cast<size_t>(std::max(0.0, std::floor(increment_p - two_standard_deviation)));

	for (size_t i = 1; i < n_nodes_in_clique_of_container; ++i) {
		auto linear_index_in_container_space_of_i = linear_start_index + i * increment;

		if (linear_index_in_container_space_of_i < linear_index_in_container_space) {
			index_in_TStorageVector[i] = index_in_TStorage;
			continue;
		}
		if (is_last_element) {
			return {current_state, exists_in_container, index_in_TStorageVector};
		}

		// calculate the upper bound
		auto upper_bound = index_in_TStorage + jump_right;
		if (upper_bound >= container.size()) { upper_bound = container.size() - 1; }
		auto upper_linear_index_in_container_space =
		    container[upper_bound].get_linear_index_in_container_space();

		// calculate the lower bound
		auto lower_bound = index_in_TStorage + jump_left;
		if (lower_bound >= container.size()) { lower_bound = container.size() - 1; }
		auto lower_linear_index_in_container_space =
		    container[lower_bound].get_linear_index_in_container_space();

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
			    binary_search(container, linear_index_in_container_space_of_i, upper_bound,
			                  container.size());
		} else if (linear_index_in_container_space_of_i > lower_linear_index_in_container_space &&
		           linear_index_in_container_space_of_i < upper_linear_index_in_container_space) {
			size_t pos = lower_bound;
			while (pos < upper_bound && container[pos].get_linear_index_in_container_space() <
			                                linear_index_in_container_space_of_i) {
				++pos;
			}
			const auto pos_lin              = container[pos].get_linear_index_in_container_space();
			found                           = (pos_lin == linear_index_in_container_space_of_i);
			index_in_TStorage               = pos;
			linear_index_in_container_space = pos_lin;
			is_last_element                 = (pos == container.size() - 1);
		} else {
			std::tie(found, index_in_TStorage, linear_index_in_container_space, is_last_element) =
			    binary_search(container, linear_index_in_container_space_of_i, index_in_TStorage,
			                  lower_bound);
		}
		index_in_TStorageVector[i] = index_in_TStorage;
		if (found) {
			current_state[i]       = container[index_in_TStorage].is_one();
			exists_in_container[i] = true;
			continue;
		}
	}
	return {current_state, exists_in_container, index_in_TStorageVector};
}

// Linear scan variant: after the initial position, scan forward element-by-element.
// Costs O(increment * density) per query on average — wins for sparse data.
template<typename Container>
CurrentStateResult hard_linear_scan(const Container &container, size_t n_nodes,
                                    const std::vector<size_t> &start_index, size_t increment,
                                    size_t /*total_size*/) {
	std::vector<uint8_t> state(n_nodes, false);
	std::vector<uint8_t> exists(n_nodes, false);
	std::vector<size_t> idx_vec(n_nodes, container.size());
	const size_t N = container.size();

	const auto lin_start         = container.get_linear_index_in_container_space(start_index);
	auto [found, idx, lin, last] = binary_search(container, lin_start, size_t{0}, lin_start + 1);
	idx_vec[0]                   = idx;
	if (found) {
		state[0]  = container[idx].is_one();
		exists[0] = true;
	}

	for (size_t i = 1; i < n_nodes; ++i) {
		if (idx >= N) return {state, exists, idx_vec};
		const auto target = lin_start + i * increment;
		if (target < lin) {
			idx_vec[i] = idx;
			continue;
		}
		while (idx < N && container[idx].get_linear_index_in_container_space() < target) { ++idx; }
		if (idx >= N) return {state, exists, idx_vec};
		lin        = container[idx].get_linear_index_in_container_space();
		last       = (idx == N - 1);
		idx_vec[i] = idx;
		if (lin == target) {
			state[i]  = container[idx].is_one();
			exists[i] = true;
		}
	}
	return {state, exists, idx_vec};
}

template<typename Container>
CurrentStateResult fill_current_state_hard(const Container &container,
                                           size_t n_nodes_in_clique_of_container,
                                           const std::vector<size_t> &start_index_in_leaves_space,
                                           size_t increment, size_t total_size_of_container) {
	const auto density =
	    static_cast<double>(container.size()) / static_cast<double>(total_size_of_container);

	if (density <= 0.001) {
		return hard_linear_scan(container, n_nodes_in_clique_of_container,
		                        start_index_in_leaves_space, increment, total_size_of_container);
	} else {
		return fill_current_state_hard_with_linear_window(container, n_nodes_in_clique_of_container,
		                                                  start_index_in_leaves_space, increment,
		                                                  total_size_of_container);
	}
}

template<bool AlongLastDim, typename Container>
CurrentStateResult fill_current_state(const Container &container,
                                      size_t n_nodes_in_clique_of_container,
                                      const std::vector<size_t> &start_index_in_leaves_space,
                                      size_t increment, size_t total_size_of_container) {
	if constexpr (AlongLastDim) {
		return fill_current_state_easy(container, start_index_in_leaves_space,
		                               n_nodes_in_clique_of_container);
	}
	return fill_current_state_hard(container, n_nodes_in_clique_of_container,
	                               start_index_in_leaves_space, increment, total_size_of_container);
}

template<typename Container>
CurrentStateResult fill_current_state(const Container &container,
                                      size_t n_nodes_in_clique_of_container,
                                      const std::vector<size_t> &start_index_in_leaves_space,
                                      size_t increment, size_t total_size_of_container) {
	if (increment == 1) {
		return fill_current_state<true>(container, n_nodes_in_clique_of_container,
		                                start_index_in_leaves_space, increment,
		                                total_size_of_container);
	}
	return fill_current_state<false>(container, n_nodes_in_clique_of_container,
	                                 start_index_in_leaves_space, increment,
	                                 total_size_of_container);
}
