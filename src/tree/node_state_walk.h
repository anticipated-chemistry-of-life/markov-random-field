//
// The Gibbs walk over one clique's node state, under that tree's own process.
//
// A chain redraws every node of the clique, leaves included. A root is drawn from the stationary
// distribution and its children. Every other internal node is drawn from its parent and its
// children. A leaf is drawn from its parent and from the link the caller bound to this clique.
//
// It is the conditional form of tree/node_state_density.h. That file scores a whole configuration,
// one term per node. This one redraws one node from the terms that name it. So a test can compare
// the walk with the density rather than with a second copy of the arithmetic. The draw beside it
// is already checked that way.
//
// It reads a topology, a process, a run of bins, a stream of uniforms and a link, and it writes
// through a column. A tree stands up three live stattools parameters and a command line; a
// phylogeny and a transition grid are values, so this runs in a test.
//

#pragma once

#include "coretools/Math/TSumLog.h"
#include "random/TCellUniforms.h"
#include "random/two_state_draw.h"
#include "tree/TPhylogeny.h"
#include "tree/branch/TTransitionGrid.h"
#include "tree/node_state_density.h"
#include "tree/node_state_draw.h"

#include <array>
#include <concepts>
#include <cstddef>

namespace node_state_walk {

using node_state_density::BranchBins;
using node_state_density::CliqueStates;
// The column the walk writes through is the one the draw writes through, so the concept is stated
// once and beside the draw.
using node_state_draw::CliqueColumn;

/// What stands below one leaf of this clique, under each of the leaf's two states.
///
/// A leaf is the one node of a clique whose Markov blanket reaches outside its own tree. The walk
/// therefore takes one extra term there, and takes it from a seam its caller binds per clique.
/// Neither the walk nor the tree that calls it knows what is behind the seam, which is what keeps
/// a tree's update free of the field, of the other tree and of the error probability (ADR-0005).
///
/// `prob_of_leaf_states(leaf)` answers `{ p(what is below | leaf = 0), p(what is below | leaf = 1)
/// }`. It is asked once per leaf, at that leaf's node index -- which is also its index in leaf
/// space (ADR-0004), so a seam that addresses leaves needs no conversion.
///
/// Both values must be strictly positive. The walk sums logs, so a zero is a state ruled out
/// rather than a state made unlikely, and it has no place in a conditional the walk draws from.
/// The link satisfies this by construction: the error probability is strictly inside (0, 0.5), so
/// the link table reaches neither 0 nor 1 (field/TFieldMath.h). A seam that has to rule a state
/// out says so with a value near zero, which the draw treats the same way to every bit that
/// matters.
template<typename T>
concept LeafLink = requires(const T &link, size_t leaf) {
	{ link.prob_of_leaf_states(leaf) } -> std::same_as<std::array<double, 2>>;
};

/// Adds `p(child | node)` to both sums, for every child of `node` and under both of its states.
///
/// A child's own branch carries the term, so `bin_of` is asked for the child and never for the
/// node. The child's state is the one the column holds now, which on a bottom-up walk is the state
/// that walk has just given it.
template<CliqueStates States, BranchBins Bins>
void add_log_prob_of_children(const TPhylogeny &topology, const TTransitionGrid &process,
                              const Bins &bin_of, const States &states, size_t node,
                              std::array<coretools::TSumLogProbability, 2> &sum_log) {
	for (const size_t child : topology.children_of(node)) {
		const size_t bin       = bin_of(child);
		const bool child_state = states.is_one(child);
		for (size_t state = 0; state < 2; ++state) {
			sum_log[state].add(process.probability(bin, state != 0, child_state));
		}
	}
}

/// Redraws every node of one clique, bottom-up.
///
/// Canonical order puts every leaf first, every child below its parent and every root last
/// (ADR-0004), so one forward pass over the node indices reaches every child before its parent.
/// There is no queue and no traversal to set up, and a tree with several roots needs no special
/// case.
///
/// The walk keeps no running density. Scoring a node against its parent *and* against each of its
/// children counts every internal branch twice, which is a conditional and not a density.
/// tree/node_state_density.h answers that question over the configuration the walk leaves behind.
template<CliqueColumn Column, BranchBins Bins, CellUniforms Uniforms, LeafLink Link>
void update_clique(const TPhylogeny &topology, const TTransitionGrid &process, const Bins &bin_of,
                   const Uniforms &uniforms, const Link &link, Column &column) {
	for (size_t node = 0; node < topology.n_nodes(); ++node) {
		std::array<coretools::TSumLogProbability, 2> sum_log;
		if (topology.is_root(node)) {
			for (size_t state = 0; state < 2; ++state) {
				sum_log[state].add(process.stationary(state != 0));
			}
		} else {
			// The parent comes after this node in canonical order, so the walk has not touched it.
			// Its state is still the one this update started from.
			const size_t bin           = bin_of(node);
			const bool state_of_parent = column.is_one(topology.parent_of(node));
			for (size_t state = 0; state < 2; ++state) {
				sum_log[state].add(process.probability(bin, state_of_parent, state != 0));
			}
		}

		// A leaf has no children, and no other node reads the link. The two terms are exclusive,
		// so the walk says which one a node takes rather than adding an empty sum.
		if (topology.is_leaf(node)) {
			const std::array<double, 2> below = link.prob_of_leaf_states(node);
			for (size_t state = 0; state < 2; ++state) { sum_log[state].add(below[state]); }
		} else {
			add_log_prob_of_children(topology, process, bin_of, column, node, sum_log);
		}

		// The cell's own uniform, and not a draw from a running generator: two containers, two
		// iterations and two thread counts then give one answer (ADR-0007).
		const bool state = two_state_draw::sample(sum_log, uniforms.at(column.linear_index(node)));
		// The column decides what a write costs. It drops a write of the state the cell already
		// carries, so reading the cell again here would only pay for that read twice.
		column.set_state(node, state);
	}
}

} // namespace node_state_walk
