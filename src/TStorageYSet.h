#pragma once

#include "TStorageY.h" // Your existing class
#include "coretools/algorithms.h"
#include <cstddef>
#include <cstdint>
#include <set>
#include <vector>

class TStorageYHolder {
public:
	mutable TStorageY _item;

	explicit TStorageYHolder(uint64_t linear_index) : _item(linear_index) {}

	// Needed for emplace and find
	bool operator<(const TStorageYHolder &other) const {
		return _item.get_linear_index_in_Y_space() < other._item.get_linear_index_in_Y_space();
	}
	bool operator<(uint64_t index) const { return _item.get_linear_index_in_Y_space() < index; }

	friend bool operator<(uint64_t index, const TStorageYHolder &holder) {
		return index < holder._item.get_linear_index_in_Y_space();
	}
};

struct Comparator {
	using is_transparent = void; // Enables heterogeneous lookup

	bool operator()(const TStorageYHolder &lhs, const TStorageYHolder &rhs) const { return lhs < rhs; }
	bool operator()(const TStorageYHolder &lhs, uint64_t rhs) const { return lhs < rhs; }
	bool operator()(uint64_t lhs, const TStorageYHolder &rhs) const { return lhs < rhs; }
};

class TStorageYSet {
private:
	std::set<TStorageYHolder, Comparator> _set;
	size_t _thinning_factor;
	size_t _total_counts;

	/// _dimensions_Y_space is the number of leaf nodes in each dimension
	std::vector<size_t> _dimensions_Y_space;

public:
	using value_type     = uint64_t;
	using const_iterator = typename std::set<TStorageY, Comparator>::const_iterator;
	TStorageYSet()       = default;
	~TStorageYSet()      = default;
	TStorageYSet(const size_t n_iterations, const std::vector<size_t> &dimensions_Y_space) {
		initialize(n_iterations, dimensions_Y_space);
	};

	void initialize(const size_t n_iterations, const std::vector<size_t> &dimensions_Y_space) {
		constexpr uint16_t max_value = std::numeric_limits<uint16_t>::max();
		_thinning_factor             = std::ceil(static_cast<double>(n_iterations) / static_cast<double>(max_value));
		_total_counts                = n_iterations / _thinning_factor;
		_dimensions_Y_space          = dimensions_Y_space;
	}
	// Insert or find
	TStorageY &operator[](uint64_t linear_index) {
		auto [it, inserted] = _set.emplace(linear_index);
		return it->_item;
	}

	bool contains(uint64_t linear_index) const { return _set.find(linear_index) != _set.end(); }

	std::set<TStorageYHolder, Comparator>::const_iterator find(uint64_t linear_index) const {
		return _set.find(linear_index);
	}

	std::set<TStorageYHolder, Comparator>::iterator find(uint64_t linear_index) { return _set.find(linear_index); }

	void erase(uint64_t linear_index) {
		auto it = _set.find(linear_index);
		if (it != _set.end()) { _set.erase(it); }
	}

	auto begin() { return _set.begin(); }
	auto end() { return _set.end(); }
	auto begin() const { return _set.begin(); }
	auto end() const { return _set.end(); }

	size_t size() const { return _set.size(); }

	void clear() { _set.clear(); }

	void insert_one(uint64_t linear_index_in_Y_space) {
		auto [it, inserted] = _set.emplace(linear_index_in_Y_space);
		if (inserted) {
			it->_item.set_state(true);
			it->_item.set_linear_index_in_Y_space(linear_index_in_Y_space);
		} else {
			it->_item.set_state(true); // Ensure the state is set to true
		}
	}

	void insert_zero(uint64_t linear_index_in_Y_space) {
		auto [it, inserted] = _set.emplace(linear_index_in_Y_space);
		if (inserted) {
			it->_item.set_state(false);
			it->_item.set_linear_index_in_Y_space(linear_index_in_Y_space);
		} else {
			it->_item.set_state(false); // Ensure the state is set to false
		}
	}

	bool is_one(uint64_t linear_index_in_Y_space) const {
		auto it = _set.find(linear_index_in_Y_space);
		if (it != _set.end()) { return it->_item.is_one(); }
		return false; // Not found, hence not one
	}

	void add_to_counter(size_t iteration) {
		if (iteration % _thinning_factor == 0) {
			for (auto &holder : _set) { holder._item.update_counter(); }
		}
	}

	double get_fraction_of_ones(uint64_t linear_index_in_Y_space) const {
		auto it = _set.find(linear_index_in_Y_space);
		if (it != _set.end()) {
			return static_cast<double>(it->_item.get_counter()) / static_cast<double>(_total_counts);
		}
		return 0.0; // Not found, hence fraction is 0
	}

	void reset_counts() {
		for (auto &holder : _set) { holder._item.reset_counter(); }
	}

	void remove_zeros() {
		for (auto it = _set.begin(); it != _set.end();) {
			if (!it->_item.is_one()) {
				it = _set.erase(it); // Erase returns the next iterator
			} else {
				++it;
			}
		}
	}

	/// Given a multi-dimensional index, we want to get its linear index.
	/// @param multi_dim_index the multi-dimensional index
	/// @return the linear index
	uint64_t get_linear_index_in_Y_space(const std::vector<size_t> &multidim_index_in_Y_space) const {
		return coretools::getLinearIndex(multidim_index_in_Y_space, _dimensions_Y_space);
	}
	uint64_t get_linear_index_in_container_space(const std::vector<size_t> &multidim_index_in_Y_space) const {
		return get_linear_index_in_Y_space(multidim_index_in_Y_space);
	}

	/// Given a linear index, we want to get the multi-dimensional index.
	[[nodiscard]] std::vector<size_t> get_multi_dimensional_index(uint64_t linear_index_in_Y_space) const {
		const auto tmp = static_cast<size_t>(linear_index_in_Y_space);
		return coretools::getSubscripts(tmp, _dimensions_Y_space);
	}

	/// Returns the product of the dimensions in the container. This is the
	/// maximum number of ones that can be stored in the vector given
	/// the dimensions of the container. For example, if the container
	/// has dimension sizes [2, 3, 4], the maximum number of ones that can be
	/// stored in the vector is 24.
	[[nodiscard]]
	size_t total_size_of_container_space() const {
		return coretools::containerProduct(_dimensions_Y_space);
	}
};
