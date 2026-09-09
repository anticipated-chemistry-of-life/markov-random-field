//
// The clique's cells, and its branches' bins, as a test writes them.
//
// The density, the draw and the walk are three headers written against one clique column. Each of
// them was given a fixture of its own, and the three drifted apart while saying the same thing.
// They share this one instead, the way the property suites already share
// tests/phylogeny_generators.h and tests/written_uniforms.h.
//

#pragma once

#include "tree/TPhylogeny.h"
#include "written_uniforms.h"

#include <array>
#include <cstddef>
#include <vector>

namespace clique {

/// One clique's column, as a test writes it. It satisfies the column concept the draw and the walk
/// write through, and the state concept the density reads.
///
/// The linear index of a node is deliberately not the node index: a clique's cells are a strided
/// run of the node state, and the uniform a node draws is named by that stride and not by the
/// node.
class TColumn {
private:
	std::vector<bool> _states;
	size_t _offset;
	size_t _stride;

public:
	TColumn(size_t n_nodes, size_t offset, size_t stride)
	    : _states(n_nodes, false), _offset(offset), _stride(stride) {}
	/// A column of the node state's first cells, for a suite that is not about the stride.
	explicit TColumn(size_t n_nodes) : TColumn(n_nodes, 0, 1) {}

	[[nodiscard]] bool is_one(size_t node) const { return static_cast<bool>(_states[node]); }
	[[nodiscard]] size_t linear_index(size_t node) const { return _offset + node * _stride; }
	void set_state(size_t node, bool state) { _states[node] = state; }

	/// The configuration as one integer, bit `node` per node. What an enumeration counts by.
	[[nodiscard]] size_t mask() const {
		size_t mask = 0;
		for (size_t node = 0; node < _states.size(); ++node) {
			if (_states[node]) { mask |= size_t{1} << node; }
		}
		return mask;
	}

	/// The configuration `mask` names: bit `node` of it is that node's state.
	void set_from_mask(size_t mask) {
		for (size_t node = 0; node < _states.size(); ++node) {
			_states[node] = ((mask >> node) & 1U) != 0U;
		}
	}
};

/// A link that says nothing: both states of every leaf are equally supported, so a leaf is drawn
/// from its parent alone. It is what a clique with nothing under its leaves looks like, and it is
/// what lets a suite compare the walk against the node-state density, which has no link term.
struct TNoLink {
	[[nodiscard]] static std::array<double, 2> prob_of_leaf_states(size_t /*leaf*/) {
		return {1.0, 1.0};
	}
};

/// A link written out one leaf at a time, so that a test names what stands below each of them.
class TWrittenLink {
private:
	std::vector<std::array<double, 2>> _probability;

public:
	explicit TWrittenLink(size_t n_leaves) : _probability(n_leaves, {1.0, 1.0}) {}

	void set(size_t leaf, std::array<double, 2> probability) { _probability[leaf] = probability; }

	[[nodiscard]] std::array<double, 2> prob_of_leaf_states(size_t leaf) const {
		return _probability[leaf];
	}
};

/// Every branch in the same bin, which is what a test that is not about branch lengths wants.
struct TOneBin {
	size_t bin = 0;
	size_t operator()(size_t /*node*/) const { return bin; }
};

/// A different bin per branch, so a term taken from the wrong node's branch does not cancel.
struct TBinPerBranch {
	size_t n_bins = 1;
	size_t operator()(size_t node) const { return node % n_bins; }
};

/// A stream wide enough for the uniforms a column of this stride asks for.
inline uniforms::TWrittenUniforms uniforms_for(const TPhylogeny &topology, size_t offset,
                                               size_t stride, double value) {
	return uniforms::TWrittenUniforms(offset + topology.n_nodes() * stride, value);
}

} // namespace clique
