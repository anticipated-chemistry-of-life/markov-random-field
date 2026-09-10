//
// The log density of one tree's node state, under that tree's own process.
//
// The joint the chain targets factors as `p(Z_s | theta_s) * p(Z_m | theta_m) * p(Y | Z_s, Z_m,
// omega) * p(L, D | Y)` (ADR-0005). This file is the first two factors, one clique at a time.
//
// One term per node and no more: a root is scored against the stationary distribution, and every
// other node against its parent. So each branch is counted exactly once. The trace this feeds used
// to be built as a by-product of the node-state walk, which scored a node against its parent *and*
// against each of its children, and so counted every internal edge twice.
//
// It reads a topology, a process and a run of states, and nothing else. A tree stands up three live
// stattools parameters and a command line; a phylogeny and a transition grid are values, so this
// runs in a test.
//

#pragma once

#include "tree/TPhylogeny.h"
#include "tree/branch/TTransitionGrid.h"

#include <cmath>
#include <concepts>
#include <cstddef>

namespace node_state_density {

/// The states of one clique's nodes, addressed by node index. A view over the clique's column
/// satisfies this, and so does a vector a test writes by hand.
template<typename T>
concept CliqueStates = requires(const T &states, size_t node) {
	{ states.is_one(node) } -> std::same_as<bool>;
};

/// The bin a node's branch sits in, addressed by node index. Never asked of a root, which has no
/// branch.
template<typename T>
concept BranchBins = requires(const T &bins, size_t node) {
	{ bins(node) } -> std::convertible_to<size_t>;
};

/// `log p(states | process)` for one clique: the stationary term of every root, plus the
/// parent-to-node term of every other node.
///
/// The states are the ones the clique holds now, and the bins are the ones its branches sit in
/// now, so the answer is the density of the configuration as it stands rather than of the states
/// as they were drawn.
template<CliqueStates States, BranchBins Bins>
[[nodiscard]] double log_density_of_clique(const TPhylogeny &topology,
                                           const TTransitionGrid &process, const States &states,
                                           const Bins &bin_of) {
	double sum = 0.0;
	for (size_t node = 0; node < topology.n_nodes(); ++node) {
		const bool state = states.is_one(node);
		if (topology.is_root(node)) {
			sum += std::log(process.stationary(state));
		} else {
			const bool parent_state = states.is_one(topology.parent_of(node));
			sum += std::log(process.probability(bin_of(node), parent_state, state));
		}
	}
	return sum;
}

} // namespace node_state_density
