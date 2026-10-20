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

/// One parameter of the scripted adapter: one value per clique, proposed and decided as the test
/// writes.
///
/// It records every log-likelihood ratio it is handed, in call order. It fails the running test if
/// the tree update asks for its accepted value while a proposal is pending, or reads or decides a
/// proposal that is not pending.
class TScriptedParameter {
public:
	/// One decision: the clique, and the ratio the tree update handed over.
	struct TDecision {
		size_t clique;
		double log_ratio;
		bool operator==(const TDecision &) const = default;
	};

private:
	const char *_name;
	std::vector<double> _accepted;
	std::vector<double> _proposals;
	std::vector<bool> _accept;

	std::vector<bool> _pending;
	size_t _n_pending = 0;

	size_t _n_proposals = 0;
	std::vector<TDecision> _decisions;

public:
	/// @param name       the parameter's name in failure messages
	/// @param accepted   each clique's value before the update
	/// @param proposals  the value a proposal puts up for each clique
	/// @param accept     whether a decision accepts each clique's proposal
	TScriptedParameter(const char *name, std::vector<double> accepted,
	                   std::vector<double> proposals, std::vector<bool> accept)
	    : _name(name), _accepted(std::move(accepted)), _proposals(std::move(proposals)),
	      _accept(std::move(accept)), _pending(_accepted.size(), false) {}

	void propose() {
		if (_n_pending != 0) { ADD_FAILURE() << _name << " proposed while a proposal is pending"; }
		_pending.assign(_pending.size(), true);
		_n_pending = _pending.size();
		++_n_proposals;
	}

	[[nodiscard]] double proposed(size_t clique) const {
		if (!_pending.at(clique)) {
			ADD_FAILURE() << "proposed " << _name << " of clique " << clique
			              << " read with no proposal";
		}
		return _proposals.at(clique);
	}

	[[nodiscard]] double accepted(size_t clique) const {
		if (_n_pending != 0) {
			ADD_FAILURE() << "accepted " << _name << " read while a proposal is pending";
		}
		return _accepted.at(clique);
	}

	bool decide(size_t clique, double log_ratio) {
		if (!_pending.at(clique)) {
			ADD_FAILURE() << _name << " of clique " << clique << " decided with no proposal";
		}
		_pending[clique] = false;
		--_n_pending;
		_decisions.push_back({clique, log_ratio});
		const bool accepted = _accept.at(clique);
		if (accepted) { _accepted[clique] = _proposals[clique]; }
		return accepted;
	}

	[[nodiscard]] size_t n_proposals() const { return _n_proposals; }
	[[nodiscard]] const std::vector<TDecision> &decisions() const { return _decisions; }
};

/// The moves adapter a test scripts: alpha and log nu, each a `TScriptedParameter`. So every test
/// that drives the tree update through it also checks the seam's ordering rule.
class TScriptedTreeMoves {
public:
	using TDecision = TScriptedParameter::TDecision;

private:
	TScriptedParameter _alpha;
	TScriptedParameter _log_nu;

public:
	TScriptedTreeMoves(TScriptedParameter alpha, TScriptedParameter log_nu)
	    : _alpha(std::move(alpha)), _log_nu(std::move(log_nu)) {}

	void propose_alphas() { _alpha.propose(); }
	[[nodiscard]] double proposed_alpha(size_t clique) const { return _alpha.proposed(clique); }
	[[nodiscard]] double accepted_alpha(size_t clique) const { return _alpha.accepted(clique); }
	bool decide_alpha(size_t clique, double log_ratio) { return _alpha.decide(clique, log_ratio); }

	void propose_log_nus() { _log_nu.propose(); }
	[[nodiscard]] double proposed_log_nu(size_t clique) const { return _log_nu.proposed(clique); }
	[[nodiscard]] double accepted_log_nu(size_t clique) const { return _log_nu.accepted(clique); }
	bool decide_log_nu(size_t clique, double log_ratio) {
		return _log_nu.decide(clique, log_ratio);
	}

	[[nodiscard]] const TScriptedParameter &alpha() const { return _alpha; }
	[[nodiscard]] const TScriptedParameter &log_nu() const { return _log_nu; }
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
	std::vector<double> log_nu_proposals;
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
			log_nu_proposals.push_back(log_rate(rng));
			grids.set(c, TTransitionGrid(alpha[c], std::exp(log_nu[c]), bin_grid));
		}
	}

	[[nodiscard]] auto bins() const {
		return [this](size_t node) { return bin_of_node.at(node); };
	}

	[[nodiscard]] TScriptedTreeMoves moves(std::vector<bool> accept_alpha,
	                                       std::vector<bool> accept_log_nu) const {
		return {TScriptedParameter("alpha", alpha, alpha_proposals, std::move(accept_alpha)),
		        TScriptedParameter("log nu", log_nu, log_nu_proposals, std::move(accept_log_nu))};
	}

	template<bool FixZ> void run(TScriptedTreeMoves &moves) {
		tree_update::run<FixZ>(topology, cliques, Z, grids, bin_grid, bins(), uniforms, moves);
	}

	/// Clique `c`'s alpha after the update, when the alpha decisions follow `alpha_accepted`.
	[[nodiscard]] double alpha_after(size_t c) const {
		return alpha_accepted(c) ? alpha_proposals[c] : alpha[c];
	}

	/// Clique `c`'s log nu after the update, when the nu decisions follow `log_nu_accepted`.
	[[nodiscard]] double log_nu_after(size_t c) const {
		return log_nu_accepted(c) ? log_nu_proposals[c] : log_nu[c];
	}

	/// The scripts the tests use: alpha accepted where `c % 3 == 0`, nu where `c % 2 == 0`. So the
	/// cliques show every pairing of an accepted or rejected alpha with an accepted or rejected nu.
	static bool alpha_accepted(size_t c) { return c % 3 == 0; }
	static bool log_nu_accepted(size_t c) { return c % 2 == 0; }

	[[nodiscard]] TScriptedTreeMoves scripted_moves() const {
		const size_t n = cliques.n_cliques();
		std::vector<bool> accept_alpha(n);
		std::vector<bool> accept_log_nu(n);
		for (size_t c = 0; c < n; ++c) {
			accept_alpha[c]  = alpha_accepted(c);
			accept_log_nu[c] = log_nu_accepted(c);
		}
		return moves(std::move(accept_alpha), std::move(accept_log_nu));
	}
};

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

/// The density of clique `c`'s states under `grid`, from the pure module.
template<typename NodeState>
double density_of_clique(TTreeFixture<NodeState> &fx, size_t c, const TTransitionGrid &grid) {
	TCliqueView<NodeState> states(fx.Z, fx.topology, fx.cliques.index_of(c), TREE);
	return node_state_density::log_density_of_clique(fx.topology, grid, states, fx.bins());
}

// -------------------------------------------------------------------------
// Grid installation
// -------------------------------------------------------------------------

TYPED_TEST(TreeUpdate, an_accepted_alpha_installs_its_grid_and_a_rejected_one_leaves_the_grid) {
	TTreeFixture<TypeParam> fx(7);
	const size_t n = fx.cliques.n_cliques();
	// Every nu rejected, so the grid each clique holds after the update is the one alpha left.
	std::vector<bool> accept_alpha(n);
	for (size_t c = 0; c < n; ++c) { accept_alpha[c] = fx.alpha_accepted(c); }
	auto moves = fx.moves(std::move(accept_alpha), std::vector<bool>(n, false));

	fx.template run<false>(moves);

	for (size_t c = 0; c < n; ++c) {
		const double alpha = fx.alpha_after(c);
		expect_installed(fx.grids, c, TTransitionGrid(alpha, std::exp(fx.log_nu[c]), fx.bin_grid));
		EXPECT_EQ(moves.accepted_alpha(c), alpha);
	}
}

TYPED_TEST(TreeUpdate, an_accepted_nu_installs_its_grid_and_a_rejected_one_leaves_the_grid) {
	TTreeFixture<TypeParam> fx(7);
	const size_t n = fx.cliques.n_cliques();
	auto moves     = fx.scripted_moves();

	fx.template run<false>(moves);

	// An accepted nu installs the grid of the accepted alpha and the accepted nu. A rejected one
	// leaves the grid alpha left: the same two values.
	for (size_t c = 0; c < n; ++c) {
		const double log_nu = fx.log_nu_after(c);
		expect_installed(fx.grids, c,
		                 TTransitionGrid(fx.alpha_after(c), std::exp(log_nu), fx.bin_grid));
		EXPECT_EQ(moves.accepted_log_nu(c), log_nu);
	}
}

// -------------------------------------------------------------------------
// The alpha-to-nu handoff
// -------------------------------------------------------------------------

TYPED_TEST(TreeUpdate, the_nu_ratio_is_measured_against_the_grid_alpha_left) {
	TTreeFixture<TypeParam> fx(7);
	const size_t n = fx.cliques.n_cliques();
	auto moves     = fx.scripted_moves();

	fx.template run<false>(moves);

	// The walk does not run again after the alpha pass, so the states now are the states nu
	// scored. An accepted alpha leaves its candidate grid; a rejected one leaves the grid alpha
	// started from.
	const auto &decisions = moves.log_nu().decisions();
	ASSERT_EQ(decisions.size(), n);
	for (size_t c = 0; c < n; ++c) {
		EXPECT_EQ(decisions[c].clique, c);
		const double alpha       = fx.alpha_after(c);
		const double after_alpha = density_of_clique(
		    fx, c, TTransitionGrid(alpha, std::exp(fx.log_nu[c]), fx.bin_grid));
		const double candidate = density_of_clique(
		    fx, c, TTransitionGrid(alpha, std::exp(fx.log_nu_proposals[c]), fx.bin_grid));
		const double tolerance = 1e-9 * (std::abs(candidate) + std::abs(after_alpha));
		EXPECT_NEAR(decisions[c].log_ratio, candidate - after_alpha, tolerance)
		    << "clique " << c << (fx.alpha_accepted(c) ? ", alpha accepted" : ", alpha rejected");
	}
}

// -------------------------------------------------------------------------
// Thread count
// -------------------------------------------------------------------------

TYPED_TEST(TreeUpdate, one_thread_and_four_threads_give_the_same_node_states_and_ratios) {
	// A prime clique count, so a four-way split never lands evenly.
	constexpr size_t N_CLIQUES = 37;

	TTreeFixture<TypeParam> one(N_CLIQUES);
	const std::vector<size_t> before = one.Z.get_full_Z_binary_vector();
	auto moves_one                   = one.scripted_moves();
	{
		TThreadCount guard(1);
		one.template run<false>(moves_one);
	}

	TTreeFixture<TypeParam> four(N_CLIQUES);
	auto moves_four = four.scripted_moves();
	{
		TThreadCount guard(4);
		four.template run<false>(moves_four);
	}

	// The walk moved some state, so equal node states are not two untouched copies.
	EXPECT_NE(one.Z.get_full_Z_binary_vector(), before);
	EXPECT_EQ(four.Z.get_full_Z_binary_vector(), one.Z.get_full_Z_binary_vector());
	ASSERT_EQ(moves_one.alpha().decisions().size(), N_CLIQUES);
	EXPECT_EQ(moves_four.alpha().decisions(), moves_one.alpha().decisions());
	ASSERT_EQ(moves_one.log_nu().decisions().size(), N_CLIQUES);
	EXPECT_EQ(moves_four.log_nu().decisions(), moves_one.log_nu().decisions());
}

// -------------------------------------------------------------------------
// FixZ
// -------------------------------------------------------------------------

TYPED_TEST(TreeUpdate, with_fixed_node_state_the_state_is_unchanged_and_alpha_and_nu_still_move) {
	TTreeFixture<TypeParam> fx(7);
	const size_t n                   = fx.cliques.n_cliques();
	const std::vector<size_t> before = fx.Z.get_full_Z_binary_vector();
	auto moves                       = fx.scripted_moves();

	fx.template run<true>(moves);

	EXPECT_EQ(fx.Z.get_full_Z_binary_vector(), before);
	EXPECT_EQ(moves.alpha().n_proposals(), 1u);
	EXPECT_EQ(moves.log_nu().n_proposals(), 1u);
	ASSERT_EQ(moves.alpha().decisions().size(), n);
	ASSERT_EQ(moves.log_nu().decisions().size(), n);
	for (size_t c = 0; c < n; ++c) {
		EXPECT_EQ(moves.alpha().decisions()[c].clique, c);
		EXPECT_EQ(moves.log_nu().decisions()[c].clique, c);
		expect_installed(fx.grids, c,
		                 TTransitionGrid(fx.alpha_after(c), std::exp(fx.log_nu_after(c)),
		                                 fx.bin_grid));
	}
}

// -------------------------------------------------------------------------
// The scripted adapter enforces the seam's ordering rule
// -------------------------------------------------------------------------

TEST(ScriptedTreeMoves, fails_the_test_on_an_accepted_value_read_while_a_proposal_is_pending) {
	TScriptedTreeMoves moves(TScriptedParameter("alpha", {0.3, 0.4}, {0.5, 0.6}, {true, false}),
	                         TScriptedParameter("log nu", {0.0, 0.1}, {0.2, 0.3}, {false, true}));
	EXPECT_EQ(moves.accepted_alpha(0), 0.3);

	moves.propose_alphas();
	EXPECT_NONFATAL_FAILURE((void)moves.accepted_alpha(0), "pending");
	EXPECT_TRUE(moves.decide_alpha(0, 0.0));
	// One clique is still pending, so the read still fails.
	EXPECT_NONFATAL_FAILURE((void)moves.accepted_alpha(0), "pending");
	EXPECT_FALSE(moves.decide_alpha(1, 0.0));

	EXPECT_EQ(moves.accepted_alpha(0), 0.5);
	EXPECT_EQ(moves.accepted_alpha(1), 0.4);

	moves.propose_log_nus();
	EXPECT_NONFATAL_FAILURE((void)moves.accepted_log_nu(1), "pending");
	EXPECT_FALSE(moves.decide_log_nu(0, 0.0));
	EXPECT_TRUE(moves.decide_log_nu(1, 0.0));

	EXPECT_EQ(moves.accepted_log_nu(0), 0.0);
	EXPECT_EQ(moves.accepted_log_nu(1), 0.3);
}

} // namespace
