//
// One clique's node state, reduced to what its log-likelihood under a transition grid reads.
//
// The log-likelihood has one term per node (tree/node_state_density.h): a root's stationary
// probability, and every other node's probability of its state given its parent's, across the bin
// its branch sits in. A term names a bin, a parent state and a state -- or, for a root, a state --
// and nothing else. So the log-likelihood under any grid is a sum over those kinds of term: how
// many nodes fall in each, times the log of its probability. That is `4 * n_bins + 2` integers per
// clique, whatever the size of the tree.
//
// The alpha and nu moves score every clique under three grids -- the installed one, alpha's
// candidate and nu's candidate -- over one configuration: the walk runs before alpha and not again,
// and the bins move only after nu. So the alpha pass counts each clique once, scores both of its
// grids off the counts, and the nu move scores its grid off the same counts without reading a node
// state at all (ADR-0016).
//
// The counts are integers, so they come out the same whichever order the nodes are counted in, and
// a score adds the kinds in one fixed order. A score is therefore a function of the counts and the
// grid alone: the same under either backend and at any thread count.
//

#pragma once

#include "coretools/Main/TError.h"
#include "tree/TPhylogeny.h"
#include "tree/branch/TTransitionGridTable.h"
#include "tree/node_state_density.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace transition_counts {

using node_state_density::BranchBins;
using node_state_density::CliqueStates;

/// How many kinds of term a clique's log-likelihood has under a grid of `n_bins` bins: four
/// transitions per bin, and a root in either state.
[[nodiscard]] constexpr size_t n_counts(size_t n_bins) { return 4 * n_bins + 2; }

/// Where the count of nodes in `bin` whose parent is `from` and who are `to` sits.
[[nodiscard]] constexpr size_t transition_slot(size_t bin, bool from, bool to) {
	return 4 * bin + 2 * static_cast<size_t>(from) + static_cast<size_t>(to);
}

/// Where the count of roots in `state` sits, after every bin's transitions.
[[nodiscard]] constexpr size_t root_slot(size_t n_bins, bool state) {
	return 4 * n_bins + static_cast<size_t>(state);
}

/// The bin count a block of counts is laid out for.
[[nodiscard]] inline size_t n_bins_of(std::span<const uint32_t> counts) {
	DEBUG_ASSERT(counts.size() >= 2 && (counts.size() - 2) % 4 == 0);
	return (counts.size() - 2) / 4;
}

/// Counts one clique's nodes by the kind of term each adds to its log-likelihood, into `counts`,
/// which it overwrites. `counts` is laid out for the grid the bins index: `n_counts(n_bins)` long.
///
/// A branch is named by the node below it, so every non-root node is counted once by its own
/// branch, and every root once by its state. A count fits 32 bits: it is at most the node count.
template<CliqueStates States, BranchBins Bins>
void count_transitions(const TPhylogeny &topology, const Bins &bin_of, const States &states,
                       std::span<uint32_t> counts) {
	const size_t n_bins = n_bins_of(counts);
	std::fill(counts.begin(), counts.end(), uint32_t{0});
	for (const size_t node : topology.branches()) {
		const size_t bin = bin_of(node);
		DEBUG_ASSERT(bin < n_bins);
		++counts[transition_slot(bin, states.is_one(topology.parent_of(node)), states.is_one(node))];
	}
	for (const size_t root : topology.roots()) { ++counts[root_slot(n_bins, states.is_one(root))]; }
}

/// `log p(states | process)` of the clique the counts were taken from: each kind's count times the
/// log of its probability, bin by bin and then the roots.
///
/// A kind no node falls in adds nothing, and its probability is not asked for: zero times the log
/// of a zero probability is NaN, and the clique carries no such term.
template<TransitionGridLike Process>
[[nodiscard]] double log_likelihood(std::span<const uint32_t> counts, const Process &process) {
	const size_t n_bins = n_bins_of(counts);
	double sum          = 0.0;
	for (size_t bin = 0; bin < n_bins; ++bin) {
		for (const bool from : {false, true}) {
			for (const bool to : {false, true}) {
				const uint32_t n = counts[transition_slot(bin, from, to)];
				if (n != 0) {
					sum += static_cast<double>(n) * std::log(process.probability(bin, from, to));
				}
			}
		}
	}
	for (const bool state : {false, true}) {
		const uint32_t n = counts[root_slot(n_bins, state)];
		if (n != 0) { sum += static_cast<double>(n) * std::log(process.stationary(state)); }
	}
	return sum;
}

/// One block of counts per clique, in one allocation, addressed the way `TTransitionGridTable`
/// addresses its grids. A clique's block is written by the one thread that holds its view, so the
/// blocks need no lock.
class TTransitionCountTable {
private:
	size_t _per_clique = 0;
	std::vector<uint32_t> _counts;

public:
	TTransitionCountTable(size_t n_cliques, size_t n_bins)
	    : _per_clique(n_counts(n_bins)), _counts(n_cliques * _per_clique) {}

	[[nodiscard]] std::span<uint32_t> of(size_t clique) {
		DEBUG_ASSERT((clique + 1) * _per_clique <= _counts.size());
		return {_counts.data() + clique * _per_clique, _per_clique};
	}
	[[nodiscard]] std::span<const uint32_t> of(size_t clique) const {
		DEBUG_ASSERT((clique + 1) * _per_clique <= _counts.size());
		return {_counts.data() + clique * _per_clique, _per_clique};
	}
};

} // namespace transition_counts
