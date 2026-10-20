//
// The tree update: the part of an iteration that belongs to one tree (GLOSSARY.md, "Tree update").
//
// It runs the node-state walk, the Metropolis move on alpha, the Metropolis move on nu and the
// Metropolis move on the branch lengths, in this order:
//   1. Read every branch's current bin, and propose every clique's alpha.
//   2. Parallel pass: the node-state walk (unless `FixZ`), then alpha's two scores per clique. One
//      score is under the installed grid. The other is under the candidate grid, built from the
//      proposed alpha and the accepted nu.
//   3. Serial alpha decisions in clique order. Each accepted candidate grid is installed.
//   4. Propose every clique's log nu.
//   5. Parallel pass: nu's candidate score per clique, under the grid built from the accepted alpha
//      and the proposed nu. The installed score is the one the clique ended alpha on.
//   6. Serial nu decisions in clique order. Each accepted candidate grid is installed.
//   7. Draw this iteration's branch pairs, then each pair's coin in pair order, only where the bin
//      grid leaves the direction free. Each pair's proposed bins follow the bin grid's step rule.
//   8. Parallel pass: each pair's log-likelihood ratio, summed over every clique in fixed slots
//      and reduced in slot order.
//   9. Serial pair decisions in pair order.
// Every move reads the bins the branches hold at the start. The proposed bins stay in this module
// until their pair's decision, so no pass reads a proposed bin as a current one.
//
// It reads plain data and reaches the parameters only through a `TreeMoves` adapter. So a test
// builds it from a phylogeny, a clique space, a node state, a transition-grid table, a bin grid and
// cell uniforms, and scripts the adapter.
//

#pragma once

#include "cli.h"
#include "coretools/Main/TError.h"
#include "coretools/Math/TSumLog.h"
#include "random/TCellUniforms.h"
#include "tree/TPhylogeny.h"
#include "tree/branch/TBinGrid.h"
#include "tree/branch/TTransitionGrid.h"
#include "tree/branch/TTransitionGridTable.h"
#include "tree/clique/TCliqueSpace.h"
#include "tree/clique/TCliqueView.h"
#include "tree/clique/clique_traversal.h"
#include "tree/node_state_density.h"
#include "tree/node_state_walk.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <utility>
#include <vector>

namespace tree_update {

using node_state_density::BranchBins;
using node_state_density::CliqueStates;

/// What the tree update asks of the parameters it moves. It speaks in Metropolis decisions.
///
/// alpha, per clique:
/// - `propose_alphas()` proposes every clique's alpha.
/// - `proposed_alpha(c)` is the proposed alpha of clique `c`, while that proposal is pending.
/// - `accepted_alpha(c)` is the accepted alpha of clique `c`, only while no alpha proposal is
///   pending.
/// - `decide_alpha(c, log_ratio)` accepts or rejects the proposal of clique `c`, from the
///   log-likelihood ratio of the candidate against the installed grid. It answers true on accept.
///   The adapter adds the prior and Hastings terms.
///
/// nu, per clique:
/// - `propose_log_nus()` proposes every clique's log nu.
/// - `proposed_log_nu(c)` is the proposed log nu of clique `c`, while that proposal is pending.
/// - `accepted_log_nu(c)` is the accepted log nu of clique `c`, only while no nu proposal is
///   pending.
/// - `decide_log_nu(c, log_ratio)` accepts or rejects the proposal of clique `c`, from the
///   log-likelihood ratio of the candidate against the grid installed after alpha. It answers true
///   on accept. The adapter adds the prior and Hastings terms.
///
/// branch lengths, per pair of branches:
/// - `sample_branch_pairs()` draws this iteration's pairs. A pair holds two branch indices, and no
///   branch is in two pairs. With an odd branch count one branch is in no pair.
/// - `n_branch_pairs()` is the number of pairs drawn, and `branch_pair(p)` is pair `p`.
/// - `branch_coin()` draws one coin. True asks the first branch of a pair to step down.
/// - `current_bin(branch)` is the bin branch `branch` holds now.
/// - `decide_branch_pair(p, first_bin, second_bin, log_ratio)` sets the two proposed bins of pair
///   `p` and accepts or rejects them, from their log-likelihood ratio against the current bins.
///   It answers true on accept. The adapter adds the prior and Hastings terms.
///
/// Two rules hold for every adapter:
/// - The tree update calls the adapter only from serial code, never from inside a parallel pass. An
///   adapter therefore needs no locks, and the generator's call order stays fixed (ADR-0007).
/// - The tree update asks for the accepted value of a parameter only while that parameter has no
///   proposal pending. Alpha's candidate grids read the accepted nu. Nu's candidate grids read the
///   accepted alpha, after every alpha decision. A branch has a proposal only inside its pair's
///   decision.
///
/// The seam has no `oldValue`. It has no query for whether a parameter is fixed either. A parameter
/// fixed with `--<name>.update false` makes its proposal a no-op and its decision a rejection, so
/// every move always runs.
template<typename T>
concept TreeMoves = requires(T &moves, const T &held, size_t clique, double log_ratio) {
	{ moves.propose_alphas() } -> std::same_as<void>;
	{ held.proposed_alpha(clique) } -> std::convertible_to<double>;
	{ held.accepted_alpha(clique) } -> std::convertible_to<double>;
	{ moves.decide_alpha(clique, log_ratio) } -> std::same_as<bool>;
	{ held.accepted_log_nu(clique) } -> std::convertible_to<double>;
	{ moves.propose_log_nus() } -> std::same_as<void>;
	{ held.proposed_log_nu(clique) } -> std::convertible_to<double>;
	{ moves.decide_log_nu(clique, log_ratio) } -> std::same_as<bool>;
	{ moves.sample_branch_pairs() } -> std::same_as<void>;
	{ held.n_branch_pairs() } -> std::convertible_to<size_t>;
	{ held.branch_pair(clique) } -> std::convertible_to<std::pair<size_t, size_t>>;
	{ moves.branch_coin() } -> std::same_as<bool>;
	{ held.current_bin(clique) } -> std::convertible_to<size_t>;
	{ moves.decide_branch_pair(clique, clique, clique, log_ratio) } -> std::same_as<bool>;
};

namespace detail {

/// One clique's log-likelihood under each of two processes, in one pass over its nodes.
///
/// A root is scored against the stationary distribution, and every other node against its parent.
/// The terms are added per node in node order. One pass reads each state once for both scores.
template<CliqueStates States, BranchBins Bins, TransitionGridLike ProcessA,
         TransitionGridLike ProcessB>
[[nodiscard]] std::array<double, 2> clique_LL(const TPhylogeny &topology, const Bins &bin_of,
                                              const States &states, const ProcessA &a,
                                              const ProcessB &b) {
	coretools::TSumLogProbability LL_a;
	coretools::TSumLogProbability LL_b;
	for (size_t node = 0; node < topology.n_nodes(); ++node) {
		const bool state = states.is_one(node);
		if (topology.is_root(node)) {
			LL_a.add(a.stationary(state));
			LL_b.add(b.stationary(state));
		} else {
			const size_t bin           = bin_of(node);
			const bool state_of_parent = states.is_one(topology.parent_of(node));
			LL_a.add(a.probability(bin, state_of_parent, state));
			LL_b.add(b.probability(bin, state_of_parent, state));
		}
	}
	return {LL_a.getSum(), LL_b.getSum()};
}

/// One clique's log-likelihood under one process. The terms are added in the order `clique_LL`
/// above adds them, so a score of the same grid is the same double.
template<CliqueStates States, BranchBins Bins, TransitionGridLike Process>
[[nodiscard]] double clique_LL(const TPhylogeny &topology, const Bins &bin_of,
                               const States &states, const Process &process) {
	coretools::TSumLogProbability LL;
	for (size_t node = 0; node < topology.n_nodes(); ++node) {
		const bool state = states.is_one(node);
		if (topology.is_root(node)) {
			LL.add(process.stationary(state));
		} else {
			const bool state_of_parent = states.is_one(topology.parent_of(node));
			LL.add(process.probability(bin_of(node), state_of_parent, state));
		}
	}
	return LL.getSum();
}

/// How many slots the branch-length ratio is summed in before it is reduced.
///
/// A fixed number, and not one per thread. The sum runs over every clique, and a thread-indexed
/// slot holds whichever cliques the schedule handed that thread. So the last bits of the total,
/// and with them the branch-length decision, moved with the thread count. A clique belongs to slot
/// `c / cliques_per_slot` whatever thread runs it, and the slots are reduced in slot order, so the
/// total is a function of the seed. Large enough to keep a big team busy, small enough that the
/// slots are a few tens of MB.
inline constexpr size_t N_BRANCH_LL_SLOTS = 64;

/// One pair of the branch-length move: its two branches, their current bins and their proposed
/// bins. The first branch moves by the step and the second by its negation.
struct TBranchPair {
	size_t first_node;
	size_t second_node;
	size_t first_bin;
	size_t second_bin;
	size_t proposed_first_bin;
	size_t proposed_second_bin;
};

/// Adds one clique's ratio of every pair to `log_sum`: the first branch's, then the second's.
///
/// Each branch is scored at its proposed and its current bin under the clique's installed grid.
/// The move proposes bins and not a grid, so one grid stands on both sides of the ratio.
template<CliqueStates States, TransitionGridLike Process>
void add_branch_ratios(const TPhylogeny &topology, const Process &process, const States &states,
                       const std::vector<TBranchPair> &pairs,
                       std::vector<coretools::TSumLogProbability> &log_sum) {
	auto ratio = [&](size_t node, size_t current, size_t proposed) {
		const bool state_of_parent = states.is_one(topology.parent_of(node));
		const bool state           = states.is_one(node);
		const double old_p         = process.probability(current, state_of_parent, state);
		const double new_p         = process.probability(proposed, state_of_parent, state);
		return new_p / old_p;
	};
	for (size_t p = 0; p < pairs.size(); ++p) {
		const TBranchPair &pair = pairs[p];
		log_sum[p].add(ratio(pair.first_node, pair.first_bin, pair.proposed_first_bin));
		log_sum[p].add(ratio(pair.second_node, pair.second_bin, pair.proposed_second_bin));
	}
}

} // namespace detail

/// Runs the node-state walk, the alpha move, the nu move and the branch-length move of one tree.
///
/// @param topology  the tree's phylogeny
/// @param cliques   the tree's clique space. Its dimension is the tree's own.
/// @param Z         the tree's node state, of either backend. Unchanged when `FixZ` is set.
/// @param grids     one installed transition grid per clique. An accepted alpha or nu installs its
///                  grid.
/// @param bin_grid  the bin grid every candidate grid is built on, and whose step rule proposes
///                  the branch-length pairs
/// @param uniforms  this iteration's cell uniforms for the node-state walk
/// @param moves     the adapter to the parameters
///
/// The first parallel pass writes through clique views. Their deferred inserts go into `Z` after
/// the pass ends (ADR-0009). The other passes only read.
template<bool FixZ, typename NodeState, CellUniforms Uniforms, TreeMoves Moves>
void run(const TPhylogeny &topology, const TCliqueSpace<> &cliques, NodeState &Z,
         TTransitionGridTable &grids, const TBinGrid &bin_grid, const Uniforms &uniforms,
         Moves &moves) {
	const size_t n = cliques.n_cliques();

	// The adapter is asked here, in serial code, and the passes below read these copies.
	std::vector<size_t> current_bin(topology.n_branches());
	for (size_t b = 0; b < current_bin.size(); ++b) { current_bin[b] = moves.current_bin(b); }
	auto bin_of = [&](size_t node) { return current_bin[topology.branch_index(node)]; };

	moves.propose_alphas();
	std::vector<double> proposed_alpha(n);
	std::vector<double> accepted_nu(n);
	for (size_t c = 0; c < n; ++c) {
		proposed_alpha[c] = moves.proposed_alpha(c);
		accepted_nu[c]    = std::exp(moves.accepted_log_nu(c));
	}

	// --- the node-state walk, and alpha's two likelihoods ---
	auto view_factory = [&](size_t c) {
		return TCliqueView<NodeState>(Z, topology, cliques.index_of(c), cliques.dimension());
	};
	auto kernel = [&](size_t c, TCliqueView<NodeState> &states) {
		// The view lives across the scores below, so they read the states the walk assigns.
		const auto installed = grids.at(c);
		if constexpr (!FixZ) {
			node_state_walk::update_clique(topology, installed, bin_of, uniforms, states);
		}
		return detail::clique_LL(topology, bin_of, states, installed,
		                         TTransitionGrid(proposed_alpha[c], accepted_nu[c], bin_grid));
	};
	auto [alpha_LL, inserts] = clique_traversal::run(n, view_factory, kernel);
	if constexpr (!FixZ) {
		Z.insert_in_Z(inserts);
	} else {
		// With `FixZ` the pass writes nothing, so no view defers an insert.
		DEBUG_ASSERT(std::all_of(inserts.begin(), inserts.end(),
		                         [](const auto &clique_inserts) { return clique_inserts.empty(); }));
	}

	// --- alpha's decisions, and the likelihood each clique ends alpha on ---
	std::vector<double> LL_after_alpha(n);
	for (size_t c = 0; c < n; ++c) {
		if (moves.decide_alpha(c, alpha_LL[c][1] - alpha_LL[c][0])) {
			grids.set(c, TTransitionGrid(proposed_alpha[c], accepted_nu[c], bin_grid));
			LL_after_alpha[c] = alpha_LL[c][1];
		} else {
			LL_after_alpha[c] = alpha_LL[c][0];
		}
	}

	// --- nu's proposals, read in serial code after every alpha decision ---
	moves.propose_log_nus();
	std::vector<double> accepted_alpha(n);
	std::vector<double> proposed_nu(n);
	for (size_t c = 0; c < n; ++c) {
		accepted_alpha[c] = moves.accepted_alpha(c);
		proposed_nu[c]    = std::exp(moves.proposed_log_nu(c));
	}

	// --- nu's candidate likelihood ---
	// The installed score is the one each clique ended alpha on, so this pass scores one grid.
	auto nu_kernel = [&](size_t c, TCliqueView<NodeState> &states) {
		return detail::clique_LL(topology, bin_of, states,
		                         TTransitionGrid(accepted_alpha[c], proposed_nu[c], bin_grid));
	};
	auto [nu_LL, nu_inserts] = clique_traversal::run(n, view_factory, nu_kernel);
	// This pass only reads, so no view defers an insert.
	DEBUG_ASSERT(std::all_of(nu_inserts.begin(), nu_inserts.end(),
	                         [](const auto &clique_inserts) { return clique_inserts.empty(); }));

	// --- nu's decisions ---
	for (size_t c = 0; c < n; ++c) {
		if (moves.decide_log_nu(c, nu_LL[c] - LL_after_alpha[c])) {
			grids.set(c, TTransitionGrid(accepted_alpha[c], proposed_nu[c], bin_grid));
		}
	}

	// --- the branch-length pairs and their proposed bins, in serial code ---
	// The pairs are drawn first, then the coins in pair order. A coin is drawn only where the bin
	// grid leaves the direction free: at a boundary the direction is forced, and a coin drawn
	// anyway would shift the generator's stream.
	moves.sample_branch_pairs();
	const size_t n_pairs = moves.n_branch_pairs();
	std::vector<detail::TBranchPair> pairs(n_pairs);
	for (size_t p = 0; p < n_pairs; ++p) {
		const auto [first, second] = static_cast<std::pair<size_t, size_t>>(moves.branch_pair(p));
		const size_t first_bin     = current_bin[first];
		const size_t second_bin    = current_bin[second];
		const bool decrease = bin_grid.step_is_free(first_bin, second_bin) && moves.branch_coin();
		const int step      = bin_grid.step_direction(first_bin, second_bin, decrease);
		pairs[p]            = {topology.branches()[first],
		                       topology.branches()[second],
		                       first_bin,
		                       second_bin,
		                       static_cast<size_t>(static_cast<int>(first_bin) + step),
		                       static_cast<size_t>(static_cast<int>(second_bin) - step)};
	}

	// --- the branch-length ratio, summed in fixed slots ---
	// A clique belongs to slot `c / cliques_per_slot`, and `schedule(static, cliques_per_slot)`
	// hands each of those blocks to one thread. So a slot has one writer and its contents do not
	// depend on the team size.
	const size_t cliques_per_slot =
	    (n + detail::N_BRANCH_LL_SLOTS - 1) / detail::N_BRANCH_LL_SLOTS;
	std::vector<std::vector<coretools::TSumLogProbability>> slots(
	    detail::N_BRANCH_LL_SLOTS, std::vector<coretools::TSumLogProbability>(n_pairs));
#pragma omp parallel for num_threads(ProgramOptions::NUMBER_OF_THREADS) default(none)              \
    schedule(static, cliques_per_slot)                                                             \
    shared(n, cliques_per_slot, view_factory, topology, grids, pairs, slots)
	for (size_t c = 0; c < n; ++c) {
		auto states = view_factory(c);
		detail::add_branch_ratios(topology, grids.at(c), states, pairs,
		                          slots[c / cliques_per_slot]);
		DEBUG_ASSERT(states.take_deferred_inserts().empty());
	}
	std::vector<coretools::TSumLogProbability> log_ratio = slots[0];
	for (size_t slot = 1; slot < slots.size(); ++slot) {
		for (size_t p = 0; p < n_pairs; ++p) { log_ratio[p] = log_ratio[p] + slots[slot][p]; }
	}

	// --- the branch-length decisions, in pair order ---
	for (size_t p = 0; p < n_pairs; ++p) {
		const detail::TBranchPair &pair = pairs[p];
		(void)moves.decide_branch_pair(p, pair.proposed_first_bin, pair.proposed_second_bin,
		                               log_ratio[p].getSum());
	}
}

} // namespace tree_update
