//
// The tree update, driven through its seam.
//
// The tests build a tree update from plain objects: a phylogeny, a clique space, a node state, a
// transition-grid table, a bin grid, cell uniforms and a scripted moves adapter. They need no tree
// file, no command line and no stattools. They assert only on what comes out: the node states, the
// ratios the adapter receives, and the grids the table holds after the decisions.
//

#include "backend_pairings.h"
#include "phylogeny_generators.h"
#include "random/TCellUniforms.h"
#include "storages/z_storage/TStorageZDense.h"
#include "storages/z_storage/TStorageZSparse.h"
#include "thread_count.h"
#include "tree/TPhylogeny.h"
#include "tree/branch/TBinGrid.h"
#include "tree/branch/TTransitionGrid.h"
#include "tree/branch/TTransitionGridTable.h"
#include "tree/clique/TCliqueSpace.h"
#include "tree/clique/TCliqueView.h"
#include "tree/node_state_density.h"
#include "tree/node_state_shape.h"
#include "tree/tree_update.h"
#include "gtest/gtest.h"
#include "gtest/gtest-spi.h"

#include <cmath>
#include <cstddef>
#include <random>
#include <utility>
#include <vector>

namespace {

using threads::TThreadCount;

/// The dimension the tree under test owns. The other dimension holds the cliques.
constexpr size_t TREE = 0;

/// The moves adapter a test scripts.
///
/// It returns the proposals and the decisions the test writes, and it records every log-likelihood
/// ratio it is handed, in call order. It fails the running test if the tree update asks for an
/// accepted value while that parameter has a proposal pending. So every test that drives the tree
/// update through it also checks that rule of the seam.
class TScriptedTreeMoves {
public:
	/// One call to `decide_alpha`: the clique, and the ratio the tree update handed over.
	struct TDecision {
		size_t clique;
		double log_ratio;
		bool operator==(const TDecision &) const = default;
	};

private:
	std::vector<double> _accepted_alpha;
	std::vector<double> _accepted_log_nu;
	std::vector<double> _alpha_proposals;
	std::vector<bool> _accept_alpha;

	std::vector<bool> _alpha_pending;
	size_t _n_alpha_pending = 0;

	size_t _n_alpha_proposals = 0;
	std::vector<TDecision> _alpha_decisions;

public:
	/// @param accepted_alpha   each clique's alpha before the update
	/// @param accepted_log_nu  each clique's log nu before the update
	/// @param alpha_proposals  the alpha `propose_alphas` puts up for each clique
	/// @param accept_alpha     whether `decide_alpha` accepts each clique's proposal
	TScriptedTreeMoves(std::vector<double> accepted_alpha, std::vector<double> accepted_log_nu,
	                   std::vector<double> alpha_proposals, std::vector<bool> accept_alpha)
	    : _accepted_alpha(std::move(accepted_alpha)), _accepted_log_nu(std::move(accepted_log_nu)),
	      _alpha_proposals(std::move(alpha_proposals)), _accept_alpha(std::move(accept_alpha)),
	      _alpha_pending(_accepted_alpha.size(), false) {}

	void propose_alphas() {
		if (_n_alpha_pending != 0) { ADD_FAILURE() << "alpha proposed while a proposal is pending"; }
		_alpha_pending.assign(_alpha_pending.size(), true);
		_n_alpha_pending = _alpha_pending.size();
		++_n_alpha_proposals;
	}

	[[nodiscard]] double proposed_alpha(size_t clique) const {
		if (!_alpha_pending.at(clique)) {
			ADD_FAILURE() << "proposed alpha of clique " << clique << " read with no proposal";
		}
		return _alpha_proposals.at(clique);
	}

	[[nodiscard]] double accepted_alpha(size_t clique) const {
		if (_n_alpha_pending != 0) {
			ADD_FAILURE() << "accepted alpha read while an alpha proposal is pending";
		}
		return _accepted_alpha.at(clique);
	}

	bool decide_alpha(size_t clique, double log_ratio) {
		if (!_alpha_pending.at(clique)) {
			ADD_FAILURE() << "alpha of clique " << clique << " decided with no proposal";
		}
		_alpha_pending[clique] = false;
		--_n_alpha_pending;
		_alpha_decisions.push_back({clique, log_ratio});
		const bool accepted = _accept_alpha.at(clique);
		if (accepted) { _accepted_alpha[clique] = _alpha_proposals[clique]; }
		return accepted;
	}

	[[nodiscard]] double accepted_log_nu(size_t clique) const {
		return _accepted_log_nu.at(clique);
	}

	[[nodiscard]] size_t n_alpha_proposals() const { return _n_alpha_proposals; }
	[[nodiscard]] const std::vector<TDecision> &alpha_decisions() const {
		return _alpha_decisions;
	}
};
static_assert(tree_update::TreeMoves<TScriptedTreeMoves>);

/// One tree, with everything a tree update reads, and a scripted adapter over it.
///
/// The forest has several roots and an odd branch count. The other tree has a leaf count unlike
/// this tree's node count, so the node state is not square.
template<typename NodeState> struct TTreeFixture {
	static constexpr size_t N_BINS = 5;

	TPhylogeny topology;
	IndexArray leaf_counts;
	TCliqueSpace<> cliques;
	NodeState Z;
	TBinGrid bin_grid{N_BINS};
	std::vector<size_t> bin_of_node;
	TTransitionGridTable grids;
	std::vector<double> alpha;
	std::vector<double> log_nu;
	std::vector<double> alpha_proposals;
	TCellUniforms uniforms{20261020, TCellStream::node_state, 3, TREE};

	static TPhylogeny forest() {
		std::mt19937_64 rng(20261020);
		return build_phylogeny(phylo::random_forest(rng, 14, 3));
	}

	explicit TTreeFixture(size_t n_cliques)
	    : topology(forest()), leaf_counts{topology.n_leaves(), n_cliques},
	      cliques(leaf_counts, TREE), Z(node_state_dimensions(leaf_counts, TREE, topology)) {
		backends::seed_ones(Z, 7);

		std::mt19937_64 rng(4242);
		std::uniform_int_distribution<size_t> bin(0, N_BINS - 1);
		std::uniform_real_distribution<double> unit(0.1, 0.9);
		std::uniform_real_distribution<double> log_rate(-1.0, 1.0);
		bin_of_node.resize(topology.n_nodes());
		for (size_t &b : bin_of_node) { b = bin(rng); }

		const size_t n = cliques.n_cliques();
		grids.resize(n, N_BINS);
		for (size_t c = 0; c < n; ++c) {
			alpha.push_back(unit(rng));
			log_nu.push_back(log_rate(rng));
			alpha_proposals.push_back(unit(rng));
			grids.set(c, TTransitionGrid(alpha[c], std::exp(log_nu[c]), bin_grid));
		}
	}

	[[nodiscard]] auto bins() const {
		return [this](size_t node) { return bin_of_node.at(node); };
	}

	[[nodiscard]] TScriptedTreeMoves moves(std::vector<bool> accept_alpha) const {
		return {alpha, log_nu, alpha_proposals, std::move(accept_alpha)};
	}

	template<bool FixZ> std::vector<double> run(TScriptedTreeMoves &moves) {
		return tree_update::run<FixZ>(topology, cliques, Z, grids, bin_grid, bins(), uniforms,
		                              moves);
	}
};

/// Every clique accepted where `c % 3 == 0`, and every other one rejected.
std::vector<bool> every_third_accepted(size_t n_cliques) {
	std::vector<bool> accept(n_cliques);
	for (size_t c = 0; c < n_cliques; ++c) { accept[c] = c % 3 == 0; }
	return accept;
}

/// The grid at clique `c` of `table` is `expected`, in every bin and in its stationary terms.
void expect_installed(const TTransitionGridTable &table, size_t c, const TTransitionGrid &expected) {
	const TTransitionGridView installed = table.at(c);
	ASSERT_EQ(installed.n_bins(), expected.n_bins());
	for (size_t bin = 0; bin < expected.n_bins(); ++bin) {
		for (const bool from : {false, true}) {
			for (const bool to : {false, true}) {
				EXPECT_EQ(installed.probability(bin, from, to), expected.probability(bin, from, to))
				    << "clique " << c << ", bin " << bin;
			}
		}
	}
	EXPECT_EQ(installed.stationary(true), expected.stationary(true)) << "clique " << c;
	EXPECT_EQ(installed.stationary(false), expected.stationary(false)) << "clique " << c;
}

template<typename NodeState> class TreeUpdate : public ::testing::Test {};
using NodeStateBackends = ::testing::Types<TStorageZDense, TStorageZSparse>;
TYPED_TEST_SUITE(TreeUpdate, NodeStateBackends);

// -------------------------------------------------------------------------
// Grid installation
// -------------------------------------------------------------------------

TYPED_TEST(TreeUpdate, an_accepted_alpha_installs_its_grid_and_a_rejected_one_leaves_the_grid) {
	TTreeFixture<TypeParam> fx(7);
	const size_t n = fx.cliques.n_cliques();
	auto moves     = fx.moves(every_third_accepted(n));

	const std::vector<double> LL_after_alpha = fx.template run<false>(moves);

	for (size_t c = 0; c < n; ++c) {
		const bool accepted = c % 3 == 0;
		const double alpha  = accepted ? fx.alpha_proposals[c] : fx.alpha[c];
		const TTransitionGrid expected(alpha, std::exp(fx.log_nu[c]), fx.bin_grid);
		expect_installed(fx.grids, c, expected);
		EXPECT_EQ(moves.accepted_alpha(c), alpha);

		// The likelihood the clique ends alpha on is the density of its new states under the grid
		// alpha left installed.
		TCliqueView<TypeParam> states(fx.Z, fx.topology, fx.cliques.index_of(c), TREE);
		const double density = node_state_density::log_density_of_clique(fx.topology, expected,
		                                                                  states, fx.bins());
		EXPECT_NEAR(LL_after_alpha.at(c), density, 1e-9 * std::abs(density)) << "clique " << c;
	}
}

// -------------------------------------------------------------------------
// Thread count
// -------------------------------------------------------------------------

TYPED_TEST(TreeUpdate, one_thread_and_four_threads_give_the_same_node_states_and_alpha_ratios) {
	// A prime clique count, so a four-way split never lands evenly.
	constexpr size_t N_CLIQUES = 37;

	TTreeFixture<TypeParam> one(N_CLIQUES);
	const std::vector<size_t> before = one.Z.get_full_Z_binary_vector();
	auto moves_one                   = one.moves(every_third_accepted(N_CLIQUES));
	{
		TThreadCount guard(1);
		(void)one.template run<false>(moves_one);
	}

	TTreeFixture<TypeParam> four(N_CLIQUES);
	auto moves_four = four.moves(every_third_accepted(N_CLIQUES));
	{
		TThreadCount guard(4);
		(void)four.template run<false>(moves_four);
	}

	// The walk moved some state, so equal node states are not two untouched copies.
	EXPECT_NE(one.Z.get_full_Z_binary_vector(), before);
	EXPECT_EQ(four.Z.get_full_Z_binary_vector(), one.Z.get_full_Z_binary_vector());
	ASSERT_EQ(moves_one.alpha_decisions().size(), N_CLIQUES);
	EXPECT_EQ(moves_four.alpha_decisions(), moves_one.alpha_decisions());
}

// -------------------------------------------------------------------------
// FixZ
// -------------------------------------------------------------------------

TYPED_TEST(TreeUpdate, with_fixed_node_state_the_state_is_unchanged_and_alpha_still_moves) {
	TTreeFixture<TypeParam> fx(7);
	const size_t n                   = fx.cliques.n_cliques();
	const std::vector<size_t> before = fx.Z.get_full_Z_binary_vector();
	auto moves                       = fx.moves(every_third_accepted(n));

	(void)fx.template run<true>(moves);

	EXPECT_EQ(fx.Z.get_full_Z_binary_vector(), before);
	EXPECT_EQ(moves.n_alpha_proposals(), 1u);
	ASSERT_EQ(moves.alpha_decisions().size(), n);
	for (size_t c = 0; c < n; ++c) {
		EXPECT_EQ(moves.alpha_decisions()[c].clique, c);
		const double alpha = c % 3 == 0 ? fx.alpha_proposals[c] : fx.alpha[c];
		expect_installed(fx.grids, c, TTransitionGrid(alpha, std::exp(fx.log_nu[c]), fx.bin_grid));
	}
}

// -------------------------------------------------------------------------
// The scripted adapter enforces the seam's ordering rule
// -------------------------------------------------------------------------

TEST(ScriptedTreeMoves, fails_the_test_on_an_accepted_alpha_read_while_a_proposal_is_pending) {
	TScriptedTreeMoves moves({0.3, 0.4}, {0.0, 0.1}, {0.5, 0.6}, {true, false});
	EXPECT_EQ(moves.accepted_alpha(0), 0.3);

	moves.propose_alphas();
	EXPECT_NONFATAL_FAILURE((void)moves.accepted_alpha(0), "pending");
	EXPECT_TRUE(moves.decide_alpha(0, 0.0));
	// One clique is still pending, so the read still fails.
	EXPECT_NONFATAL_FAILURE((void)moves.accepted_alpha(0), "pending");
	EXPECT_FALSE(moves.decide_alpha(1, 0.0));

	EXPECT_EQ(moves.accepted_alpha(0), 0.5);
	EXPECT_EQ(moves.accepted_alpha(1), 0.4);
}

} // namespace
