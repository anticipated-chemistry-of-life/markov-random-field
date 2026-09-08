#pragma once

#include "TStorageY.h"
#include "storages/storage_concepts.h"
#include <coretools/Main/TError.h>
#include <coretools/algorithms.h>
#include <unordered_map>

class TStorageYSparse {
private:
	size_t _thinning_factor = 1;
	size_t _total_counts    = 0;
	IndexArray _dimensions{};
	std::unordered_map<size_t, TStorageY> _states{};

public:
	using TCell      = TStorageY;
	using TCellConst = const TCell;

	void initialize(size_t n_iterations, const IndexArray &dimensions) {
		_thinning_factor = std::max<size_t>(
		    1, static_cast<size_t>(std::ceil(static_cast<double>(n_iterations) /
			                                 static_cast<double>(TStorageY::MAX_COUNTER))));

		_total_counts = 0;
		_dimensions   = dimensions;
	}

	[[nodiscard]]
	size_t total_size_of_container_space() const {
		return coretools::containerProduct(_dimensions);
	}

	[[nodiscard]] bool empty() const { return _states.empty(); }
	[[nodiscard]] auto find(size_t linear_index) const { return _states.find(linear_index); }
	[[nodiscard]] auto find(const IndexArray &multidim_index) const {
		return find(get_linear_index_in_container_space(multidim_index));
	}

	[[nodiscard]] IsOneResult<TCell> is_one(size_t linear_index) {
		auto stored = _states.find(linear_index);
		if (stored == _states.end()) { return {false, false, linear_index, nullptr}; }
		return {stored->second.is_one(), true, linear_index, &stored->second};
	}

	[[nodiscard]] IsOneResult<TCell> is_one(const IndexArray &multidim_index) {
		return is_one(get_linear_index_in_container_space(multidim_index));
	}

	[[nodiscard]] IsOneResult<TCellConst> is_one(size_t linear_index) const {
		auto stored = _states.find(linear_index);
		if (stored == _states.end()) { return {false, false, linear_index, nullptr}; }
		return {stored->second.is_one(), true, linear_index, &stored->second};
	}

	[[nodiscard]] IsOneResult<TCellConst> is_one(const IndexArray &multidim_index) const {
		return is_one(get_linear_index_in_container_space(multidim_index));
	}

	[[nodiscard]] size_t
	get_linear_index_in_container_space(const IndexArray &multidim_index) const {
		return coretools::getLinearIndex(multidim_index, _dimensions);
	}

	[[nodiscard]] IndexArray get_multi_dimensional_index(size_t linear_index) const {
		return coretools::getSubscriptsAsArray(linear_index, _dimensions);
	}

	void set_state(size_t linear_index, bool state) {
		DEBUG_ASSERT(_states.find(linear_index) != _states.end());
		_states[linear_index].set_state(state);
	}
	void insert_one(size_t linear_index) {
		DEBUG_ASSERT(_states.find(linear_index) == _states.end());
		_states[linear_index] = TStorageY(true);
	}
	void insert_one(const IndexArray &multidim_index) {
		insert_one(get_linear_index_in_container_space(multidim_index));
	}
	void insert_zero(size_t linear_index) {
		DEBUG_ASSERT(_states.find(linear_index) == _states.end());
		_states[linear_index] = TStorageY(false);
	}
	void remove_zeros() {
		for (auto it = _states.begin(); it != _states.end();) {
			if (!it->second.is_one()) {
				it = _states.erase(it);
			} else {
				++it;
			}
		}
	}

	void reset_counts() {
		for (auto &[linear_index, state] : _states) { state.set_counter(0); }
		_total_counts = 0;
	}

	[[nodiscard]] double get_fraction_of_ones(size_t linear_index) const {
		if (_total_counts == 0) { return 0.0; }
		const auto it = _states.find(linear_index);
		if (it == _states.end()) { return 0.0; }
		return static_cast<double>(it->second.get_counter()) / static_cast<double>(_total_counts);
	}

	[[nodiscard]] size_t get_total_counts() const { return _total_counts; }
	[[nodiscard]] size_t get_thinning_factor() const { return _thinning_factor; }

	void add_to_counter(size_t iteration) {
		if (iteration % _thinning_factor != 0) { return; }
		for (auto &[linear_index, state] : _states) {
			if (!state.is_one()) { continue; }
			state.update_counter();
		}
		++_total_counts;
	}

	void
	insert_ones_in_container(const std::vector<std::vector<size_t>> &linear_indices_to_insert) {
		for (const auto &indices : linear_indices_to_insert) {
			for (size_t linear_index : indices) { insert_one(linear_index); }
		}
	}
};
static_assert(FieldStorage<TStorageYSparse>);
