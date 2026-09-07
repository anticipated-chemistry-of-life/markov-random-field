//
// What one tree's node state is worth, under that tree's own process.
//
// The property that matters is that it is a *density*: summed over every configuration of the
// clique it comes to exactly 1. That is what fails the moment a branch is scored twice, and it is
// why the joint density trace is trustworthy under ADR-0005 where it was not before. The suite
// enumerates small trees, so the sum is computed rather than estimated.
//
// A phylogeny and a transition grid are values, so nothing here builds a tree or a chain.
//

#include "phylogeny_generators.h"
#include "tree/TPhylogeny.h"
#include "tree/branch/TBinGrid.h"
#include "tree/branch/TTransitionGrid.h"
#include "tree/node_state_density.h"
#include "gtest/gtest.h"

#include <cmath>
#include <cstddef>
#include <vector>

namespace {

using node_state_density::log_density_of_clique;
using phylo::edge;

constexpr size_t N_BINS = 5;

/// The states of a clique, as a test writes them: one bit per node, indexed by node index.
class TStates {
private:
	std::vector<bool> _states;

public:
	explicit TStates(size_t n_nodes) : _states(n_nodes, false) {}

	/// The configuration `mask` names: bit `node` of it is that node's state.
	void set_from_mask(size_t mask) {
		for (size_t node = 0; node < _states.size(); ++node) {
			_states[node] = ((mask >> node) & 1U) != 0U;
		}
	}
	void set(size_t node, bool state) { _states[node] = state; }

	[[nodiscard]] bool is_one(size_t node) const { return _states[node]; }
};

/// Every branch in the same bin, which is what a test that is not about branch lengths wants.
struct TOneBin {
	size_t bin = 0;
	size_t operator()(size_t /*node*/) const { return bin; }
};

/// The total probability the clique's density assigns to every configuration of `topology`.
double total_probability(const TPhylogeny &topology, const TTransitionGrid &process,
                         const TOneBin &bins) {
	TStates states(topology.n_nodes());
	double total = 0.0;
	for (size_t mask = 0; mask < (size_t{1} << topology.n_nodes()); ++mask) {
		states.set_from_mask(mask);
		total += std::exp(log_density_of_clique(topology, process, states, bins));
	}
	return total;
}

// -------------------------------------------------------------------------
// It is a density
// -------------------------------------------------------------------------

/// The one property the old trace did not have. A term per branch and a term per root sums to 1
/// over the configurations; scoring a node against its children as well does not.
TEST(NodeStateDensity, sums_to_one_over_every_configuration_of_a_tree) {
	const TBinGrid grid(N_BINS);
	const TPhylogeny topology = build_phylogeny(
	    {edge("a", "inner"), edge("b", "inner"), edge("inner", "root"), edge("c", "root")});

	for (const double alpha : {0.1, 0.5, 0.8}) {
		for (const double nu : {0.2, 1.0, 4.0}) {
			for (const size_t bin : {size_t{0}, N_BINS - 1}) {
				SCOPED_TRACE("alpha=" + std::to_string(alpha) + " nu=" + std::to_string(nu) +
				             " bin=" + std::to_string(bin));
				const TTransitionGrid process(alpha, nu, grid);
				EXPECT_NEAR(total_probability(topology, process, TOneBin{bin}), 1.0, 1e-12);
			}
		}
	}
}

/// A tree may have more than one root, and each root's subtree is drawn independently from the
/// stationary distribution. The sum is still 1, which says the roots are scored once each.
TEST(NodeStateDensity, sums_to_one_over_a_forest_with_several_roots) {
	const TBinGrid grid(N_BINS);
	const TTransitionGrid process(0.3, 1.5, grid);
	const TPhylogeny topology = build_phylogeny(
	    {edge("a", "r1"), edge("b", "r1"), edge("c", "r2"), edge("d", "r3"), edge("e", "r3")});
	ASSERT_EQ(topology.n_roots(), 3u);

	EXPECT_NEAR(total_probability(topology, process, TOneBin{2}), 1.0, 1e-12);
}

/// The shapes a balanced fixture never produces: a chain, where every internal node has one child,
/// and a star, where one node has all of them.
TEST(NodeStateDensity, sums_to_one_over_a_chain_and_over_a_star) {
	const TBinGrid grid(N_BINS);
	const TTransitionGrid process(0.6, 0.7, grid);

	EXPECT_NEAR(total_probability(build_phylogeny(phylo::chain(7)), process, TOneBin{1}), 1.0,
	            1e-12);
	EXPECT_NEAR(total_probability(build_phylogeny(phylo::star(6)), process, TOneBin{3}), 1.0,
	            1e-12);
}

// -------------------------------------------------------------------------
// It is the density it claims to be
// -------------------------------------------------------------------------

/// Written out by hand for the smallest tree there is: one root and one leaf below it.
TEST(NodeStateDensity, is_the_stationary_term_of_the_root_and_the_branch_term_of_the_leaf) {
	const TBinGrid grid(N_BINS);
	const TTransitionGrid process(0.25, 1.1, grid);
	const TPhylogeny topology = build_phylogeny({edge("leaf", "root")});

	TStates states(topology.n_nodes());
	states.set(topology.index_of("root"), true);
	states.set(topology.index_of("leaf"), false);

	const double expected =
	    std::log(process.stationary(true)) + std::log(process.probability(2, true, false));
	EXPECT_NEAR(log_density_of_clique(topology, process, states, TOneBin{2}), expected, 1e-12);
}

/// A root has no branch, so nothing about it moves when the branch lengths do. A leaf's does.
TEST(NodeStateDensity, moves_with_the_bin_a_branch_sits_in) {
	const TBinGrid grid(N_BINS);
	const TTransitionGrid process(0.25, 1.1, grid);
	const TPhylogeny topology = build_phylogeny({edge("leaf", "root")});

	TStates states(topology.n_nodes());
	states.set(topology.index_of("root"), true);
	states.set(topology.index_of("leaf"), false);

	EXPECT_NE(log_density_of_clique(topology, process, states, TOneBin{0}),
	          log_density_of_clique(topology, process, states, TOneBin{N_BINS - 1}));
}

} // namespace
