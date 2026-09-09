//
// The top-down draw of one clique's node state, under that tree's own process.
//
// A simulated chain draws each tree's whole node state first: every root from the stationary
// distribution, every other node from its parent, leaves included. The leaf block of what this
// leaves behind is that tree's tree field (ADR-0005).
//
// It is the exact inverse of tree/node_state_density.h: that file scores a configuration under one
// term per node, and this one draws a configuration from the same terms. A test can therefore
// check the draw against the density rather than against a second copy of the arithmetic.
//
// It reads a topology, a process, a run of bins and a stream of uniforms, and it writes through a
// column. A tree stands up three live stattools parameters and a command line; a phylogeny and a
// transition grid are values, so this runs in a test.
//

#pragma once

#include "random/TCellUniforms.h"
#include "tree/TPhylogeny.h"
#include "tree/branch/TTransitionGrid.h"
#include "tree/node_state_density.h"

#include <concepts>
#include <cstddef>

namespace node_state_draw {

using node_state_density::BranchBins;

/// The cells of one clique's column, as the draw reads and writes them. A view over the column
/// satisfies this, and so does a vector a test writes by hand.
///
/// `linear_index` is what names a cell to the stream of uniforms, so the state a node is given
/// does not move when a traversal reaches it in another order (ADR-0007).
template<typename T>
concept CliqueColumn = requires(T &column, const T &const_column, size_t node, bool state) {
	{ const_column.is_one(node) } -> std::same_as<bool>;
	{ const_column.linear_index(node) } -> std::same_as<size_t>;
	{ column.set_state(node, state) } -> std::same_as<void>;
};

/// Draws every node of one clique, top-down: a root from the stationary distribution, every other
/// node from the state its parent was just given.
///
/// Canonical order puts every child below its parent (ADR-0004), so one backward pass over the
/// node indices reaches every parent before its children. There is no queue and no traversal to
/// set up, and a tree with several roots needs no special case.
template<CliqueColumn Column, BranchBins Bins, CellUniforms Uniforms>
void draw_clique(const TPhylogeny &topology, const TTransitionGrid &process, const Bins &bin_of,
                 const Uniforms &uniforms, Column &column) {
	for (size_t node = topology.n_nodes(); node-- > 0;) {
		const double probability_of_one =
		    topology.is_root(node)
		        ? process.stationary(true)
		        : process.probability(bin_of(node), column.is_one(topology.parent_of(node)),
		                              /*to =*/true);
		// The cell's own uniform, and not a draw from a running generator: two containers, two
		// iterations and two thread counts then give one answer (ADR-0007).
		const bool state = uniforms.at(column.linear_index(node)) < probability_of_one;
		// A write of the state the cell already carries is dropped. A sparse column would
		// otherwise buffer an insert for a cell it does not hold, and then throw that insert away.
		if (column.is_one(node) != state) { column.set_state(node, state); }
	}
}

} // namespace node_state_draw
