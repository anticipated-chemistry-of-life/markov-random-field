//
// The tree update: the part of an iteration that belongs to one tree (GLOSSARY.md, "Tree update").
//
// It runs the node-state walk and the Metropolis move on alpha, in this order:
//   1. Propose every clique's alpha.
//   2. Parallel pass: the node-state walk (unless `FixZ`), then alpha's two scores per clique. One
//      score is under the installed grid. The other is under the candidate grid, built from the
//      proposed alpha and the accepted nu.
//   3. Serial alpha decisions in clique order. Each accepted candidate grid is installed.
// The nu move and the branch-length move still run in TTree, after this.
//
// It reads plain data and reaches the parameters only through a `TreeMoves` adapter. So a test
// builds it from a phylogeny, a clique space, a node state, a transition-grid table, a bin grid and
// cell uniforms, and scripts the adapter.
//

#pragma once

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
/// - `accepted_log_nu(c)` is the accepted log nu of clique `c`.
///
/// Two rules hold for every adapter:
/// - The tree update calls the adapter only from serial code, never from inside a parallel pass. An
///   adapter therefore needs no locks, and the generator's call order stays fixed (ADR-0007).
/// - The tree update asks for the accepted value of a parameter only while that parameter has no
///   proposal pending. Alpha's candidate grids read the accepted nu.
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

} // namespace detail

/// Runs the node-state walk and the alpha move of one tree.
///
/// @param topology  the tree's phylogeny
/// @param cliques   the tree's clique space. Its dimension is the tree's own.
/// @param Z         the tree's node state, of either backend. Unchanged when `FixZ` is set.
/// @param grids     one installed transition grid per clique. An accepted alpha installs its grid.
/// @param bin_grid  the bin grid every candidate grid is built on
/// @param bin_of    the bin each branch holds now, addressed by node index
/// @param uniforms  this iteration's cell uniforms for the node-state walk
/// @param moves     the adapter to the parameters
/// @return per clique, the log-likelihood of the clique under the grid alpha left installed. The
///         nu move divides by it.
///
/// The parallel pass writes through clique views. Their deferred inserts go into `Z` after the
/// pass ends (ADR-0009).
template<bool FixZ, typename NodeState, BranchBins Bins, CellUniforms Uniforms,
         TreeMoves Moves>
[[nodiscard]] std::vector<double> run(const TPhylogeny &topology, const TCliqueSpace<> &cliques,
                                      NodeState &Z, TTransitionGridTable &grids,
                                      const TBinGrid &bin_grid, const Bins &bin_of,
                                      const Uniforms &uniforms, Moves &moves) {
	const size_t n = cliques.n_cliques();

	// The adapter is asked here, in serial code, and the pass below reads these copies.
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
	return LL_after_alpha;
}

} // namespace tree_update
