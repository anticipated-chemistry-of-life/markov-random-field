//
// How often each cell of one tree field was a one.
//
// The field carries its own counter, inside its storage. A tree field does not: it is the leaf
// block of a node state, and a node state has no counter at all (storage_concepts.h). This holds
// one beside it, over the leaf-pair space the tree field shares with the field.
//
// It is a counter of its own rather than a counter added to the node-state storages because only
// the leaf block is wanted. A counter inside the storage would carry every internal node too --
// the species node state alone is `n_nodes(species) x n_leaves(molecules)` -- to report a block
// the size of the field.
//
// Why it is wanted: `omega` and the two alphas trade against each other along a ridge at constant
// `a~_s * a~_m` (ADR-0005, derivation 3). The two tree fields moving in opposite directions is
// what shows up first when a chain rides that ridge, and it cannot be seen from the field, whose
// density is the product and stays put.
//
// The thinning factor is the field's own, so both files count the same iterations and their
// denominators agree.
//

#pragma once

#include "constants.h"
#include "coretools/Main/TError.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

/// The posterior of one tree field: how often each leaf pair held a one, and over how many counted
/// iterations.
class TTreeFieldPosterior {
private:
	std::vector<uint16_t> _counts;
	/// The number of iterations actually counted, and so the largest a cell's counter can be. It
	/// is counted rather than derived from the chain length, for the reasons `FieldStorage` gives
	/// in storages/storage_concepts.h.
	size_t _total_counts    = 0;
	size_t _thinning_factor = 1;

public:
	/// The largest value a counter can hold. A chain thinned by the field's thinning factor cannot
	/// reach past it: the field sizes that factor from its own 15-bit counter, so a counter of 16
	/// bits has room to spare here.
	static constexpr uint16_t MAX_COUNTER = std::numeric_limits<uint16_t>::max();

	/// @param n_cells         the leaf-pair space, which is the field's container space
	/// @param thinning_factor the field's own, so that both posteriors count the same iterations
	void initialize(size_t n_cells, size_t thinning_factor) {
		if (thinning_factor == 0) {
			throw coretools::TDevError("A tree field posterior needs a thinning factor of at least "
			                           "one.");
		}
		_counts.assign(n_cells, 0);
		_total_counts    = 0;
		_thinning_factor = thinning_factor;
	}

	/// Counts every leaf pair the tree field holds a one at, on one iteration in the thinning
	/// factor.
	///
	/// A leaf pair sits at the same `(row, column)` in the field and in either node state
	/// (ADR-0005), so the field's own index conversion names the cell in both. Reading the shape
	/// off the field rather than restating it is what keeps this right whichever storage backs
	/// either of them.
	template<typename Field, typename NodeState>
	void add_to_counter(size_t iteration, const Field &Y, const NodeState &Z) {
		if (iteration % _thinning_factor != 0) { return; }
		for (size_t cell = 0; cell < _counts.size(); ++cell) {
			if (!Z.is_one(
			        Z.get_linear_index_in_container_space(Y.get_multi_dimensional_index(cell)))) {
				continue;
			}
			if (_counts[cell] == MAX_COUNTER) {
				throw coretools::TDevError("A tree field counter exceeds the 16-bit maximum (",
				                           MAX_COUNTER, ")");
			}
			++_counts[cell];
		}
		// After the loop, so a counter that overflows does not leave the denominator counting an
		// iteration the cells never got.
		++_total_counts;
	}

	void reset_counts() {
		std::fill(_counts.begin(), _counts.end(), 0);
		_total_counts = 0;
	}

	[[nodiscard]] size_t size() const { return _counts.size(); }

	/// How often the cell was counted a one -- the numerator of the fraction below.
	[[nodiscard]] uint16_t get_counter(size_t cell) const {
		DEBUG_ASSERT(cell < _counts.size());
		return _counts[cell];
	}

	[[nodiscard]] double get_fraction_of_ones(size_t cell) const {
		DEBUG_ASSERT(cell < _counts.size());
		// Nothing counted yet: no posterior to report, and nothing to divide by.
		if (_total_counts == 0) { return 0.0; }
		return static_cast<double>(_counts[cell]) / static_cast<double>(_total_counts);
	}

	[[nodiscard]] size_t get_total_counts() const { return _total_counts; }
	[[nodiscard]] size_t get_thinning_factor() const { return _thinning_factor; }
};
