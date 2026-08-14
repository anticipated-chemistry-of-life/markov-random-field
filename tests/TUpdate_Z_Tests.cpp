#include "TClique.h"
#include "coretools/Main/TRandomGenerator.h"
#include "gtest/gtest.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

/// Correctness tests for the exact block sampler of a clique's internal states
/// (sample_clique_states_exact in TClique.h).
///
/// The sampler claims to draw from p(Z | leaf states, alpha, nu) exactly. On a tree small
/// enough to enumerate, that claim is directly checkable: compute the target by brute force
/// over all 2^(#internal) configurations, draw many times, and compare frequencies.
///
/// The tree is described with plain vectors, so nothing here depends on TTree, TCurrentState,
/// TStorageZMatrix or the stattools parameter graph.
namespace {

/// A hand-built tree. Node 0 is the root.
///
///            0                 internal: 0, 1, 2
///          /   \               leaves:   3, 4, 5, 6
///         1     2
///        / \   / \
///       3   4 5   6
struct TinyTree {
	static constexpr size_t n_nodes    = 7;
	static constexpr size_t n_internal = 3;

	std::vector<std::vector<size_t>> children{{1, 2}, {3, 4}, {5, 6}, {}, {}, {}, {}};
	std::vector<size_t> parent{0, 0, 0, 1, 1, 2, 2}; // parent[0] unused
	std::vector<size_t> bin{0, 3, 7, 1, 5, 2, 9};    // branch above each node; bin[0] unused
	// internal nodes, children before parents
	std::vector<size_t> post_order{1, 2, 0};

	[[nodiscard]] bool is_leaf(size_t v) const { return children[v].empty(); }
	[[nodiscard]] bool is_root(size_t v) const { return v == 0; }
};

/// A three-node chain, to cover a topology with a single child per node and a deeper nesting.
///
///   0 -> 1 -> 2 -> 3(leaf)
struct ChainTree {
	static constexpr size_t n_nodes    = 4;
	static constexpr size_t n_internal = 3;

	std::vector<std::vector<size_t>> children{{1}, {2}, {3}, {}};
	std::vector<size_t> parent{0, 0, 1, 2};
	std::vector<size_t> bin{0, 4, 11, 6};
	std::vector<size_t> post_order{2, 1, 0};

	[[nodiscard]] bool is_leaf(size_t v) const { return children[v].empty(); }
	[[nodiscard]] bool is_root(size_t v) const { return v == 0; }
};

/// A root with a single leaf child, used to check the root's stationary marginal.
///
///   0 -> 1(leaf)
struct OneLeafTree {
	static constexpr size_t n_nodes    = 2;
	static constexpr size_t n_internal = 1;

	std::vector<std::vector<size_t>> children{{1}, {}};
	std::vector<size_t> parent{0, 0};
	std::vector<size_t> bin{0, 4};
	std::vector<size_t> post_order{0};

	[[nodiscard]] bool is_leaf(size_t v) const { return children[v].empty(); }
	[[nodiscard]] static bool is_root(size_t v) { return v == 0; }
};

TMatrices make_matrices(size_t n_bins, double alpha, double nu) {
	const double delta = 2.0 / (static_cast<double>(n_bins) + 1.0);
	TMatrices matrices(n_bins, delta);
	matrices.set_lambda(alpha, TypeNu(nu));
	return matrices;
}

/// Brute-force p(Z | leaf states) over every configuration of the internal nodes.
template<typename Tree>
std::vector<double> exact_posterior(const Tree &tree, const TMatrices &matrices, double alpha,
                                    const std::vector<bool> &leaf_states) {
	const size_t n_config = 1u << Tree::n_internal;
	std::vector<double> weight(n_config, 0.0);

	for (size_t config = 0; config < n_config; ++config) {
		// bit i of config is the state of post_order[i]; leaves take their observed state
		std::vector<bool> state(Tree::n_nodes, false);
		for (size_t i = 0; i < Tree::n_internal; ++i) {
			state[tree.post_order[i]] = ((config >> i) & 1u) != 0u;
		}
		for (size_t v = 0; v < Tree::n_nodes; ++v) {
			if (tree.is_leaf(v)) { state[v] = leaf_states[v]; }
		}

		double w = state[0] ? alpha : (1.0 - alpha); // root stationary
		for (size_t v = 0; v < Tree::n_nodes; ++v) {
			if (tree.is_root(v)) { continue; }
			w *= matrices[tree.bin[v]](state[tree.parent[v]], state[v]);
		}
		weight[config] = w;
	}

	double total = 0.0;
	for (double w : weight) { total += w; }
	for (double &w : weight) { w /= total; }
	return weight;
}

/// Draw n_draws times and return the empirical distribution over internal configurations,
/// encoded the same way as exact_posterior.
template<typename Tree>
std::vector<double> empirical_posterior(const Tree &tree, const TMatrices &matrices, double alpha,
                                        const std::vector<bool> &leaf_states, size_t n_draws) {
	const size_t n_config = 1u << Tree::n_internal;
	std::vector<double> count(n_config, 0.0);

	std::vector<bool> state = leaf_states;
	state.resize(Tree::n_nodes, false);
	std::vector<std::array<double, 2>> log_L;

	for (size_t draw = 0; draw < n_draws; ++draw) {
		sample_clique_states_exact(
		    tree.post_order, Tree::n_nodes, matrices, alpha,
		    [&](size_t v) -> const std::vector<size_t> & { return tree.children[v]; },
		    [&](size_t v) { return tree.is_leaf(v); }, [&](size_t v) { return tree.parent[v]; },
		    [&](size_t v) { return tree.is_root(v); }, [&](size_t v) { return tree.bin[v]; },
		    [&](size_t v) { return static_cast<bool>(state[v]); },
		    [&](size_t v, bool s) { state[v] = s; }, log_L);

		size_t config = 0;
		for (size_t i = 0; i < Tree::n_internal; ++i) {
			if (state[tree.post_order[i]]) { config |= (1u << i); }
		}
		count[config] += 1.0;
	}
	for (double &c : count) { c /= static_cast<double>(n_draws); }
	return count;
}

/// Largest absolute difference between the exact and empirical distributions, in units of the
/// Monte Carlo standard error of a proportion. Comparing on that scale means the tolerance
/// does not have to be retuned when n_draws changes.
double max_z_score(const std::vector<double> &exact, const std::vector<double> &empirical,
                   size_t n_draws) {
	double worst = 0.0;
	for (size_t i = 0; i < exact.size(); ++i) {
		const double se = std::sqrt(std::max(exact[i] * (1.0 - exact[i]), 1e-12) /
		                            static_cast<double>(n_draws));
		worst           = std::max(worst, std::abs(empirical[i] - exact[i]) / se);
	}
	return worst;
}

constexpr size_t N_DRAWS = 200000;

//----------------------------------------------------------------------------------------
// Joint (Y, Z) block: every node is sampled, leaves carry an emission likelihood.
//----------------------------------------------------------------------------------------

/// Emission table: log e_v(s) for every node and state. Internal nodes get 0 (no emission).
using Emissions = std::vector<std::array<double, 2>>;

/// Deterministic pseudo-random emissions for the leaves, so the tests are reproducible and do
/// not accidentally pick a degenerate table.
template<typename Tree> Emissions make_leaf_emissions(const Tree &tree, unsigned seed) {
	Emissions e(Tree::n_nodes, {0.0, 0.0});
	unsigned x = seed;
	auto next  = [&x]() {
		x = x * 1664525u + 1013904223u;
		// keep away from 0 and 1 so the logs stay finite
		return 0.05 + 0.9 * static_cast<double>((x >> 8) % 1000u) / 1000.0;
	};
	for (size_t v = 0; v < Tree::n_nodes; ++v) {
		if (!tree.is_leaf(v)) { continue; }
		e[v][0] = std::log(next());
		e[v][1] = std::log(next());
	}
	return e;
}

/// Brute-force p(all states | emissions) over every configuration of every node.
template<typename Tree>
std::vector<double> exact_posterior_with_emissions(const Tree &tree, const TMatrices &matrices,
                                                   double alpha, const Emissions &log_emission) {
	const size_t n_config = 1u << Tree::n_nodes;
	std::vector<double> weight(n_config, 0.0);

	for (size_t config = 0; config < n_config; ++config) {
		std::vector<bool> state(Tree::n_nodes, false);
		for (size_t v = 0; v < Tree::n_nodes; ++v) { state[v] = ((config >> v) & 1u) != 0u; }

		double log_w = std::log(state[0] ? alpha : (1.0 - alpha)); // root stationary
		for (size_t v = 0; v < Tree::n_nodes; ++v) {
			if (!tree.is_root(v)) {
				log_w += std::log(matrices[tree.bin[v]](state[tree.parent[v]], state[v]));
			}
			log_w += log_emission[v][state[v] ? 1 : 0];
		}
		weight[config] = std::exp(log_w);
	}

	double total = 0.0;
	for (double w : weight) { total += w; }
	for (double &w : weight) { w /= total; }
	return weight;
}

/// Empirical distribution over full configurations from the joint sampler.
template<typename Tree>
std::vector<double> empirical_posterior_with_emissions(const Tree &tree, const TMatrices &matrices,
                                                       double alpha, const Emissions &log_emission,
                                                       const std::vector<size_t> &post_order_all,
                                                       size_t n_draws) {
	const size_t n_config = 1u << Tree::n_nodes;
	std::vector<double> count(n_config, 0.0);

	std::vector<char> state(Tree::n_nodes, 0);
	std::vector<std::array<double, 2>> log_L;

	for (size_t draw = 0; draw < n_draws; ++draw) {
		sample_clique_states_with_emissions(
		    post_order_all, Tree::n_nodes, matrices, alpha,
		    [&](size_t v) -> const std::vector<size_t> & { return tree.children[v]; },
		    [&](size_t v) { return tree.parent[v]; }, [&](size_t v) { return tree.is_root(v); },
		    [&](size_t v) { return tree.bin[v]; },
		    [&](size_t v, bool s) { return log_emission[v][s ? 1 : 0]; },
		    [&](size_t v) { return state[v] != 0; },
		    [&](size_t v, bool s) { state[v] = static_cast<char>(s); }, log_L);

		size_t config = 0;
		for (size_t v = 0; v < Tree::n_nodes; ++v) {
			if (state[v] != 0) { config |= (1u << v); }
		}
		count[config] += 1.0;
	}
	for (double &c : count) { c /= static_cast<double>(n_draws); }
	return count;
}

} // namespace

TEST(Update_Z_block, matches_exact_posterior_on_a_balanced_tree) {
	coretools::instances::randomGenerator().setSeed(42, true);

	const TinyTree tree;
	const double alpha    = 0.35;
	const auto matrices   = make_matrices(20, alpha, 0.8);
	// leaves are nodes 3..6; entries for internal nodes are ignored
	const std::vector<bool> leaves{false, false, false, true, false, true, true};

	const auto exact     = exact_posterior(tree, matrices, alpha, leaves);
	const auto empirical = empirical_posterior(tree, matrices, alpha, leaves, N_DRAWS);

	// sanity: the enumeration really does cover a non-degenerate distribution
	double max_p = 0.0;
	for (double p : exact) { max_p = std::max(max_p, p); }
	EXPECT_LT(max_p, 0.99) << "test configuration is degenerate, it proves nothing";

	EXPECT_LT(max_z_score(exact, empirical, N_DRAWS), 5.0)
	    << "block sampler does not reproduce p(Z | Y)";
}

TEST(Update_Z_block, matches_exact_posterior_with_a_high_switching_rate) {
	coretools::instances::randomGenerator().setSeed(7, true);

	const TinyTree tree;
	const double alpha    = 0.65;
	// large nu: branches are close to independent, so the leaves say little about the internals
	const auto matrices   = make_matrices(20, alpha, 8.0);
	const std::vector<bool> leaves{false, false, false, true, true, false, true};

	const auto exact     = exact_posterior(tree, matrices, alpha, leaves);
	const auto empirical = empirical_posterior(tree, matrices, alpha, leaves, N_DRAWS);

	EXPECT_LT(max_z_score(exact, empirical, N_DRAWS), 5.0)
	    << "block sampler does not reproduce p(Z | Y)";
}

TEST(Update_Z_block, matches_exact_posterior_with_a_low_switching_rate) {
	coretools::instances::randomGenerator().setSeed(11, true);

	const TinyTree tree;
	// deliberately not 0.5: at alpha = 0.5 the transition matrix is symmetric, so the test would
	// be blind to a transposed lookup
	const double alpha    = 0.35;
	// small nu: strong parent-child coupling, which is the regime the real chain drifts into
	const auto matrices   = make_matrices(20, alpha, 0.05);
	const std::vector<bool> leaves{false, false, false, true, true, true, false};

	const auto exact     = exact_posterior(tree, matrices, alpha, leaves);
	const auto empirical = empirical_posterior(tree, matrices, alpha, leaves, N_DRAWS);

	EXPECT_LT(max_z_score(exact, empirical, N_DRAWS), 5.0)
	    << "block sampler does not reproduce p(Z | Y)";
}

TEST(Update_Z_block, matches_exact_posterior_on_a_chain) {
	coretools::instances::randomGenerator().setSeed(3, true);

	const ChainTree tree;
	const double alpha    = 0.4;
	const auto matrices   = make_matrices(20, alpha, 1.3);
	const std::vector<bool> leaves{false, false, false, true};

	const auto exact     = exact_posterior(tree, matrices, alpha, leaves);
	const auto empirical = empirical_posterior(tree, matrices, alpha, leaves, N_DRAWS);

	EXPECT_LT(max_z_score(exact, empirical, N_DRAWS), 5.0)
	    << "block sampler does not reproduce p(Z | Y)";
}

TEST(Update_Z_block, respects_the_stationary_distribution_at_the_root) {
	// With a single leaf hanging off the root and a switching rate high enough that the code
	// uses the stationary matrices, every node is an independent draw from (1-alpha, alpha).
	// The root marginal must then be exactly alpha, independent of the leaf state.
	coretools::instances::randomGenerator().setSeed(5, true);

	const OneLeafTree tree;
	const double alpha  = 0.3;
	const auto matrices = make_matrices(20, alpha, 30.0); // above the stationary threshold

	for (bool leaf_state : {false, true}) {
		const std::vector<bool> leaves{false, leaf_state};
		const auto exact     = exact_posterior(tree, matrices, alpha, leaves);
		const auto empirical = empirical_posterior(tree, matrices, alpha, leaves, N_DRAWS);
		EXPECT_NEAR(exact[1], alpha, 1e-12) << "enumeration disagrees with the stationary law";
		EXPECT_LT(max_z_score(exact, empirical, N_DRAWS), 5.0);
	}
}

//----------------------------------------------------------------------------------------
// Joint (Y, Z) block tests
//----------------------------------------------------------------------------------------

TEST(Update_YZ_block, matches_exact_joint_posterior_on_a_balanced_tree) {
	coretools::instances::randomGenerator().setSeed(101, true);

	const TinyTree tree;
	const double alpha  = 0.35;
	const auto matrices = make_matrices(20, alpha, 0.8);
	const auto emission = make_leaf_emissions(tree, 12345u);
	// all nodes, children before parents
	const std::vector<size_t> post_order_all{3, 4, 5, 6, 1, 2, 0};

	const auto exact = exact_posterior_with_emissions(tree, matrices, alpha, emission);
	const auto empirical =
	    empirical_posterior_with_emissions(tree, matrices, alpha, emission, post_order_all, N_DRAWS);

	double max_p = 0.0;
	for (double p : exact) { max_p = std::max(max_p, p); }
	EXPECT_LT(max_p, 0.5) << "test configuration is degenerate, it proves nothing";

	EXPECT_LT(max_z_score(exact, empirical, N_DRAWS), 5.0)
	    << "joint (Y, Z) block does not reproduce p(Y, Z | emissions)";
}

TEST(Update_YZ_block, matches_exact_joint_posterior_with_a_low_switching_rate) {
	// strong parent-child coupling: the emissions and the tree pull against each other
	coretools::instances::randomGenerator().setSeed(202, true);

	const TinyTree tree;
	const double alpha  = 0.6;
	const auto matrices = make_matrices(20, alpha, 0.05);
	const auto emission = make_leaf_emissions(tree, 999u);
	const std::vector<size_t> post_order_all{3, 4, 5, 6, 1, 2, 0};

	const auto exact = exact_posterior_with_emissions(tree, matrices, alpha, emission);
	const auto empirical =
	    empirical_posterior_with_emissions(tree, matrices, alpha, emission, post_order_all, N_DRAWS);

	EXPECT_LT(max_z_score(exact, empirical, N_DRAWS), 5.0)
	    << "joint (Y, Z) block does not reproduce p(Y, Z | emissions)";
}

TEST(Update_YZ_block, matches_exact_joint_posterior_on_a_chain) {
	coretools::instances::randomGenerator().setSeed(303, true);

	const ChainTree tree;
	const double alpha  = 0.45;
	const auto matrices = make_matrices(20, alpha, 1.3);
	const auto emission = make_leaf_emissions(tree, 777u);
	const std::vector<size_t> post_order_all{3, 2, 1, 0};

	const auto exact = exact_posterior_with_emissions(tree, matrices, alpha, emission);
	const auto empirical =
	    empirical_posterior_with_emissions(tree, matrices, alpha, emission, post_order_all, N_DRAWS);

	EXPECT_LT(max_z_score(exact, empirical, N_DRAWS), 5.0)
	    << "joint (Y, Z) block does not reproduce p(Y, Z | emissions)";
}

TEST(Update_YZ_block, reduces_to_the_Z_only_sampler_when_emissions_pin_the_leaves) {
	// A sharp emission that all but fixes a leaf must reproduce the observed-leaf sampler.
	// This is the bridge between the two functions: same target, different interface.
	coretools::instances::randomGenerator().setSeed(404, true);

	const TinyTree tree;
	const double alpha             = 0.4;
	const auto matrices            = make_matrices(20, alpha, 0.7);
	const std::vector<bool> leaves = {false, false, false, true, false, true, true};

	Emissions emission(TinyTree::n_nodes, {0.0, 0.0});
	for (size_t v = 0; v < TinyTree::n_nodes; ++v) {
		if (!tree.is_leaf(v)) { continue; }
		emission[v][leaves[v] ? 1 : 0] = std::log(1.0);
		emission[v][leaves[v] ? 0 : 1] = std::log(1e-12); // effectively impossible
	}
	const std::vector<size_t> post_order_all{3, 4, 5, 6, 1, 2, 0};

	const auto joint = empirical_posterior_with_emissions(tree, matrices, alpha, emission,
	                                                      post_order_all, N_DRAWS);
	// marginalise the joint draws down to the internal configuration, encoded as in
	// exact_posterior (bit i is post_order[i], i.e. nodes 1, 2, 0)
	std::vector<double> joint_internal(1u << TinyTree::n_internal, 0.0);
	for (size_t config = 0; config < joint.size(); ++config) {
		size_t internal = 0;
		for (size_t i = 0; i < TinyTree::n_internal; ++i) {
			if (((config >> tree.post_order[i]) & 1u) != 0u) { internal |= (1u << i); }
		}
		joint_internal[internal] += joint[config];
	}

	const auto exact_z = exact_posterior(tree, matrices, alpha, leaves);
	EXPECT_LT(max_z_score(exact_z, joint_internal, N_DRAWS), 5.0)
	    << "with the leaves pinned, the joint block must agree with the Z-only sampler";
}
