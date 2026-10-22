//
// One clique's node state, reduced to its transition counts, and scored off them.
//
// The property that matters is that a score read off the counts is the clique's log density: the
// one tree/node_state_density.h computes term by term over the nodes. The alpha and nu moves
// decide on scores read off the counts, so a count in the wrong slot, a branch read from the wrong
// node or a root left out is a wrong Metropolis ratio. The suite asserts that equality over random
// forests, random configurations, a bin per branch and grids on either side of the stationary
// short-circuit -- and with one count serving every grid, the way the tree update uses it.
//

#include "clique_columns.h"
#include "phylogeny_generators.h"
#include "tree/TPhylogeny.h"
#include "tree/branch/TBinGrid.h"
#include "tree/branch/TTransitionGrid.h"
#include "tree/node_state_density.h"
#include "tree/transition_counts.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <span>
#include <vector>

namespace {

using clique::TColumn;
using node_state_density::log_density_of_clique;
using transition_counts::count_transitions;
using transition_counts::log_likelihood;
using transition_counts::n_counts;

/// A bin per branch, drawn at random, so that a term read from the wrong node's branch does not
/// land in the right bin by accident. A root has no branch and is never asked.
struct TRandomBins {
	std::vector<size_t> bin_of_node;
	size_t operator()(size_t node) const { return bin_of_node[node]; }
};

TRandomBins random_bins(std::mt19937_64 &rng, const TPhylogeny &topology, size_t n_bins) {
	std::uniform_int_distribution<size_t> pick(0, n_bins - 1);
	TRandomBins bins{std::vector<size_t>(topology.n_nodes())};
	for (const size_t node : topology.branches()) { bins.bin_of_node[node] = pick(rng); }
	return bins;
}

TColumn random_states(std::mt19937_64 &rng, const TPhylogeny &topology, double p_one) {
	std::bernoulli_distribution one(p_one);
	TColumn states(topology.n_nodes());
	for (size_t node = 0; node < topology.n_nodes(); ++node) { states.set_state(node, one(rng)); }
	return states;
}

/// Grids over one bin grid: alphas near both ends and in the middle, and rates from slow decay up
/// past the threshold above which the grid is the stationary distribution in every bin.
std::vector<TTransitionGrid> grids_over(std::mt19937_64 &rng, const TBinGrid &bin_grid) {
	std::uniform_real_distribution<double> alpha(0.01, 0.99);
	std::uniform_real_distribution<double> log_nu(std::log(0.01), std::log(60.0));
	std::vector<TTransitionGrid> grids;
	for (size_t g = 0; g < 6; ++g) {
		grids.emplace_back(alpha(rng), std::exp(log_nu(rng)), bin_grid);
	}
	grids.emplace_back(0.002, 0.05, bin_grid);
	grids.emplace_back(0.3, 40.0, bin_grid);
	return grids;
}

/// The two scores of one configuration under every grid, from one count.
void expect_scores_match(const TPhylogeny &topology, const TColumn &states, const TRandomBins &bins,
                         const std::vector<TTransitionGrid> &grids, size_t n_bins) {
	std::vector<uint32_t> counts(n_counts(n_bins));
	count_transitions(topology, bins, states, std::span<uint32_t>(counts));
	for (size_t g = 0; g < grids.size(); ++g) {
		const double expected = log_density_of_clique(topology, grids[g], states, bins);
		const double scored   = log_likelihood(std::span<const uint32_t>(counts), grids[g]);
		ASSERT_TRUE(std::isfinite(expected)) << "grid " << g;
		EXPECT_NEAR(scored, expected, 1e-10 * std::max(1.0, std::abs(expected))) << "grid " << g;
	}
}

// -------------------------------------------------------------------------
// The score off the counts is the clique's log density
// -------------------------------------------------------------------------

TEST(TransitionCounts, score_the_log_density_of_random_forests_under_random_grids) {
	std::mt19937_64 rng(20261022);
	std::uniform_int_distribution<size_t> pick_bins(1, 12);
	std::uniform_int_distribution<size_t> pick_roots(1, 4);
	std::uniform_int_distribution<size_t> pick_extra_nodes(0, 60);
	std::uniform_real_distribution<double> pick_p_one(0.0, 1.0);

	for (size_t trial = 0; trial < 200; ++trial) {
		SCOPED_TRACE("trial " + std::to_string(trial));
		const size_t n_roots      = pick_roots(rng);
		const TPhylogeny topology = build_phylogeny(
		    phylo::random_forest(rng, 2 * n_roots + pick_extra_nodes(rng), n_roots));
		const size_t n_bins = pick_bins(rng);
		const TBinGrid bin_grid(n_bins);

		const TColumn states   = random_states(rng, topology, pick_p_one(rng));
		const TRandomBins bins = random_bins(rng, topology, n_bins);
		expect_scores_match(topology, states, bins, grids_over(rng, bin_grid), n_bins);
	}
}

/// The shapes a random forest seldom produces: a chain, where every internal node has one child,
/// and a star, where one node has all of them. Each is scored with every node a zero, every node a
/// one, and a random mix.
TEST(TransitionCounts, score_the_log_density_of_a_chain_and_of_a_star) {
	std::mt19937_64 rng(20261023);
	constexpr size_t N_BINS = 7;
	const TBinGrid bin_grid(N_BINS);
	const std::array<TPhylogeny, 2> shapes{build_phylogeny(phylo::chain(40)),
	                                       build_phylogeny(phylo::star(150))};
	for (size_t s = 0; s < shapes.size(); ++s) {
		for (const double p_one : {0.0, 1.0, 0.4}) {
			SCOPED_TRACE("shape " + std::to_string(s) + " p_one=" + std::to_string(p_one));
			const TColumn states   = random_states(rng, shapes[s], p_one);
			const TRandomBins bins = random_bins(rng, shapes[s], N_BINS);
			expect_scores_match(shapes[s], states, bins, grids_over(rng, bin_grid), N_BINS);
		}
	}
}

// -------------------------------------------------------------------------
// A kind no node falls in
// -------------------------------------------------------------------------

/// A process that rules one transition out entirely: a parent in state 1 never has a child in
/// state 0, in any bin.
struct TNoLossProcess {
	[[nodiscard]] double probability(size_t /*bin*/, bool from, bool to) const {
		if (from) { return to ? 1.0 : 0.0; }
		return to ? 0.25 : 0.75;
	}
	[[nodiscard]] double stationary(bool state) const { return state ? 0.4 : 0.6; }
};

/// A clique none of whose nodes makes the ruled-out transition is a configuration the process
/// allows, so its score is finite, and the density's. Zero nodes times the log of a zero
/// probability would be NaN.
TEST(TransitionCounts, a_probability_of_zero_no_node_reaches_leaves_the_score_finite) {
	constexpr size_t N_BINS = 3;
	// root -> inner -> {a, b}, root -> c, every node a one: no parent in state 1 has a child in
	// state 0, so the ruled-out transition has a count of zero.
	const TPhylogeny topology =
	    build_phylogeny({phylo::edge("a", "inner"), phylo::edge("b", "inner"),
		                 phylo::edge("inner", "root"), phylo::edge("c", "root")});
	TColumn states(topology.n_nodes());
	for (size_t node = 0; node < topology.n_nodes(); ++node) { states.set_state(node, true); }
	const clique::TBinPerBranch bins{N_BINS};

	std::vector<uint32_t> counts(n_counts(N_BINS));
	count_transitions(topology, bins, states, std::span<uint32_t>(counts));
	const double scored = log_likelihood(std::span<const uint32_t>(counts), TNoLossProcess{});

	EXPECT_TRUE(std::isfinite(scored));
	EXPECT_NEAR(scored, log_density_of_clique(topology, TNoLossProcess{}, states, bins), 1e-12);
}

} // namespace
