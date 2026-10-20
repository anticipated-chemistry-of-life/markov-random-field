//
// Created by Marco Visani on 26.06.23.
//

#ifndef METABOLITE_INFERENCE_TREE_H
#define METABOLITE_INFERENCE_TREE_H

#include "Types.h"
#include "cli.h"
#include "constants.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Files/TOutputFile.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Math/TSumLog.h"
#include "coretools/algorithms.h"
#include "omp.h"
#include "random/TCellUniforms.h"
#include "stattools/ParametersObservations/TParameter.h"
#include "storages/storage_backend.h"
#include "tree/TPhylogeny.h"
#include "tree/TStattoolsTreeMoves.h"
#include "tree/branch/TBinGrid.h"
#include "tree/branch/TTransitionGrid.h"
#include "tree/branch/TTransitionGridTable.h"
#include "tree/clique/TCliqueSpace.h"
#include "tree/clique/TCliqueView.h"
#include "tree/clique/clique_traversal.h"
#include "tree/node_state_density.h"
#include "tree/tree_update.h"
#include <algorithm>
#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

/// Note: All indices are within the tree itself
class TTree : public stattools::prior::TStochasticBase<stattools::TParameterBase, TypeMarkovField,
                                                       NumDimMarkovField> {
public:
	// some type aliases, for better readability
	using BoxType = TTree;
	using Base    = stattools::prior::TStochasticBase<stattools::TParameterBase, TypeMarkovField,
	                                                  NumDimMarkovField>;
	using typename Base::Storage;
	using typename Base::UpdatedStorage;

	using TypeParamAlpha       = stattools::TParameter<SpecAlpha, TTree>;
	using TypeParamLogNu       = stattools::TParameter<SpecLogNu, TTree>;
	using TypeParamBinBranches = stattools::TParameter<SpecBinnedBranches, TTree>;

private:
	std::string _tree_name;

	// The topology, and everything derived from it. Set once by _load_from_file; every use goes
	// through _topology(), which throws rather than reading a half-built tree.
	std::optional<TPhylogeny> _phylogeny;

	[[nodiscard]] const TPhylogeny &_topology() const { return _phylogeny.value(); }

	// dimension of the tree
	size_t _dimension;

	// For binning branch lengths. The grid is set once by _initialize_grid_branch_lengths, during
	// _load_from_file; every use goes through _grid(), which throws rather than reading a
	// half-built grid if that ever stops being true.
	std::optional<TBinGrid> _bin_grid;
	std::vector<size_t> _binned_branch_lengths_from_tree;
	TypeParamBinBranches *_binned_branch_lengths = nullptr;

	[[nodiscard]] const TBinGrid &_grid() const { return _bin_grid.value(); }

	/// One transition grid per clique, in clique order. A clique is its grid. Nothing else about a
	/// clique is stored: its multidimensional index follows from its position here, and its cells
	/// follow from that index (tree/clique/TCliqueView.h).
	///
	/// A grid is installed once the parameters exist (guessInitialValues) and replaced wholesale
	/// whenever a proposal on alpha or nu is accepted. There is no mutable "try" copy. A grid is
	/// empty until then, so asking for one before the parameters are drawn throws instead of
	/// reading a grid of zeros.
	TTransitionGridTable _transition_grids;

	// The cliques of this tree, numbered. Set once by _initialize_cliques; every use goes through
	// _cliques(), which throws rather than reading extents nothing has filled in yet.
	std::optional<TCliqueSpace<>> _clique_space;
	std::vector<std::string> _clique_names;

	[[nodiscard]] const TCliqueSpace<> &_cliques() const { return _clique_space.value(); }

	// Nus
	TypeParamLogNu *_log_nu_c = nullptr;

	// Alphas
	TypeParamAlpha *_alpha_c = nullptr;

	// Set Z
	TNodeStateStorage _Z;

	// private functions
	void _set_initial_branch_lengths(bool is_simulation);
	[[nodiscard]] std::vector<size_t>
	_bin_branch_lengths(const std::vector<double> &branch_lengths) const;
	void _bin_branch_lengths_from_tree();
	void _initialize_grid_branch_lengths();
	void _initialize_Z(IndexArray num_leaves_per_tree,
	                   const std::vector<std::unique_ptr<TTree>> &all_trees);
	void _initialize_cliques(const IndexArray &num_leaves_per_tree,
	                         const std::vector<std::unique_ptr<TTree>> &all_trees);
	/// @brief Load tree from file
	void _load_from_file(const std::string &filename, const std::string &tree_name);

	// updating branch lengths
	[[nodiscard]] stattools::TPairIndexSampler _build_pairs_branch_lengths() const;
	void _propose_new_branch_lengths(const stattools::TPairIndexSampler &pairs);
	void _propose_new_branch_lengths(size_t p1, size_t p2, int val);
	void _add_to_LL_branch_lengths(size_t c, const TNodeStateCliqueView &states,
	                               std::vector<coretools::TSumLogProbability> &log_sum,
	                               const stattools::TPairIndexSampler &pairs) const;
	template<TransitionGridLike Process>
	[[nodiscard]] double
	_calculate_likelihood_ratio_branch_length(size_t index_in_binned_branch_length,
	                                          const Process &process,
	                                          const TNodeStateCliqueView &states) const;

	/// The multidimensional index of clique `c`: a leaf in every dimension but this tree's own,
	/// which carries a 0. Clique space owns the convention (ADR-0011).
	[[nodiscard]] IndexArray _clique_index(size_t c) const;

	/// Installs the process of clique `c`. Called once the parameters exist, and again whenever a
	/// proposal on alpha or nu is accepted. Nothing outside the tree installs a grid.
	void _set_transition_grid(size_t c, const TTransitionGrid &grid) {
		_transition_grids.set(c, grid);
	}

	/// P(node | parent) under an explicitly given process, so a Metropolis proposal can ask the
	/// same question of the clique's current grid and of its candidate.
	template<TransitionGridLike Process>
	[[nodiscard]] double
	_prob_to_parent(size_t index_in_tree, TypeBinnedBranchLengths binned_branch_length,
	                const TNodeStateCliqueView &states, const Process &process) const {
		const size_t parent_index = _topology().parent_of(index_in_tree);
		const bool parent_state   = states.is_one(parent_index);
		const bool child_state    = states.is_one(index_in_tree);
		return process.probability(binned_branch_length, parent_state, child_state);
	}

	/// The bottom-up start of one clique. It shares its child terms with the node-state walk.
	void _initialize_clique_from_children(size_t c, TNodeStateCliqueView &states) const;
	template<TransitionGridLike Process>
	void _initialize_node_from_children(size_t node_index, const Process &process,
	                                    TNodeStateCliqueView &states) const;

	/// The bin every branch of this tree holds now. The node-state walk, the alpha and nu moves,
	/// the chain start and the node-state density are handed this, so none of them reads a
	/// parameter.
	[[nodiscard]] auto _bins() const {
		return [this](size_t node) { return get_binned_branch_length(node); };
	}

	void _simulateUnderPrior(Storage *) override;

	/// One node's contribution to a clique's log-likelihood under `process`: scored against its
	/// parent, or against the stationary distribution if it is a root.
	template<TransitionGridLike Process>
	void
	_add_node_to_clique_LL(size_t index_in_tree, bool state_of_node,
	                       coretools::TSumLogProbability &LL, const TNodeStateCliqueView &states,
	                       std::optional<size_t> branch_len_bin, const Process &process) const {
		if (_topology().is_root(index_in_tree)) {
			LL.add(process.stationary(state_of_node));
		} else {
			LL.add(_prob_to_parent(index_in_tree, branch_len_bin.value(), states, process));
		}
	}

	/// One clique's log-likelihood under one process. The nu move needs only this. Its ratio divides
	/// by the likelihood `tree_update::run` returns for the clique. That sum is over the same nodes
	/// in the same order, so it is the same double and not a near one.
	template<TransitionGridLike Process>
	[[nodiscard]] double _clique_LL(const TNodeStateCliqueView &states,
	                                const Process &process) const {
		coretools::TSumLogProbability LL;
		const auto &topology = _topology();
		for (size_t i = 0; i < topology.n_nodes(); ++i) {
			const bool state_of_node = states.is_one(i);
			std::optional<size_t> branch_len_bin;
			if (!topology.is_root(i)) { branch_len_bin = get_binned_branch_length(i); }
			_add_node_to_clique_LL(i, state_of_node, LL, states, branch_len_bin, process);
		}
		return LL.getSum();
	}

	/// The grid a proposed nu would give clique `c`, with the accepted alpha. Ask for it before
	/// `acceptOrReject` puts a rejected value back.
	[[nodiscard]] TTransitionGrid _nu_candidate_grid(size_t c, double proposed) const {
		return TTransitionGrid(_alpha_c->value(c), std::exp(proposed), _grid());
	}

	/// The cells of clique `c`, addressed by node index. One place holds the three things a view
	/// is built from, so the loops below cannot drift apart.
	[[nodiscard]] TNodeStateCliqueView _clique_view(size_t c) {
		return {_Z, _topology(), _clique_index(c), _dimension};
	}

	void _evalute_update_branch_length(std::vector<coretools::TSumLogProbability> &log_sum,
	                                   const stattools::TPairIndexSampler &pairs);

	/// How many slots the branch-length likelihood is summed in before it is reduced.
	///
	/// A fixed number, and not one per thread. The sum runs over every clique, and a
	/// thread-indexed slot holds whichever cliques the schedule handed that thread -- so the last
	/// bits of the total, and with them the branch-length move's accept or reject, moved with the
	/// thread count and with `schedule(dynamic)`'s run-to-run assignment. A clique belongs to slot
	/// `c * N / n_cliques` whatever runs it, and the slots are reduced in slot order, so the total
	/// is a function of the seed. Large enough to keep a big team busy, small enough that the
	/// slots are a few tens of MB.
	static constexpr size_t N_BRANCH_LL_SLOTS = 64;

	/// The slots above, reduced in slot order.
	[[nodiscard]] static std::vector<coretools::TSumLogProbability>
	_reduce_branch_LL_slots(std::vector<std::vector<coretools::TSumLogProbability>> &slots,
	                        size_t n_pairs) {
		auto reduced = slots[0];
		for (size_t slot = 1; slot < slots.size(); ++slot) {
			for (size_t p = 0; p < n_pairs; ++p) { reduced[p] = reduced[p] + slots[slot][p]; }
		}
		return reduced;
	}

public:
	TTree(size_t dimension, const std::string &filename, const std::string &tree_name,
	      TypeParamAlpha *Alpha, TypeParamLogNu *LogNu,
	      TypeParamBinBranches *Binned_Branch_Lenghts);
	~TTree() override;

	/// The topology this tree is built on. Everything that does not need the parameters -- the
	/// current state, the sheet, the clique's walks -- should ask this rather than the tree.
	[[nodiscard]] const TPhylogeny &phylogeny() const { return _topology(); }

	/** @return the number of leaves in the tree
	 */
	[[nodiscard]] size_t get_number_of_leaves() const { return _topology().n_leaves(); }

	/** @param node_index: the index of the node within the tree
	 * @return The index of the node in leaf space. Meaningless if the node is not a leaf; the node
	 * index alone says whether it is one.
	 */
	[[nodiscard]] size_t get_index_within_leaves(size_t node_index) const {
		return _topology().leaf_index(node_index);
	}
	[[nodiscard]] size_t get_index_within_leaves(const std::string &node_name) const {
		return _topology().leaf_index(_topology().index_of(node_name));
	}
	[[nodiscard]] size_t get_node_index_from_leaf_index(size_t leaf_index) const {
		return _topology().leaves()[leaf_index];
	}

	// stattools stuff
	[[nodiscard]] std::string name() const override;
	void initialize() override;
	[[nodiscard]] double getSumLogPriorDensity(const Storage &) const override;
	void guessInitialValues() override;
	[[nodiscard]] double getDensity(const Storage &, size_t) const override;
	[[nodiscard]] double getLogDensityRatio(const UpdatedStorage &, size_t) const override;

	void initialize_cliques_and_Z(const std::vector<std::unique_ptr<TTree>> &all_trees);

	/// The number of cliques of this tree: the product of every other tree's leaf count.
	///
	/// Clique space answers it, and the grid vector is sized from the same answer. The count is a
	/// property of clique space and not of how many grids happen to be held.
	[[nodiscard]] size_t n_cliques() const { return _cliques().n_cliques(); }

	/// The process of clique `c`.
	[[nodiscard]] TTransitionGridView transition_grid(size_t c) const {
		return _transition_grids.at(c);
	}

	/// The process of the clique a cell belongs to. A clique of this tree is named by a leaf of
	/// every other tree, so the cell's own dimension is dropped on the way in.
	[[nodiscard]] TTransitionGridView
	transition_grid_of_cell(const IndexArray &index_in_leaves_space) const {
		return transition_grid(_cliques().clique_of(index_in_leaves_space));
	};

	/// The clique a cell belongs to. A clique of this tree is named by a leaf of every other tree,
	/// so the cell's own dimension is dropped on the way in -- which is what makes the clique of
	/// every cell of one field row the same, in the tree that does not own the row's dimension.
	/// The block update leans on that (field/TBlockUpdate.h) and checks it in a debug build.
	[[nodiscard]] size_t clique_of_cell(const IndexArray &index_in_leaves_space) const {
		return _cliques().clique_of(index_in_leaves_space);
	}

	[[nodiscard]] const TNodeStateStorage &get_Z() const;
	TNodeStateStorage &get_Z();

	[[nodiscard]] std::string get_node_id(size_t index) const { return _topology().id_of(index); }

	/// The tree update (GLOSSARY.md): the node-state update, then the moves on alpha, nu and the
	/// branch lengths. `tree_update::run` runs the node-state walk and the alpha move. The nu move
	/// and the branch-length move follow here.
	///
	/// The field is not named here at all. Each tree owns a leaf-level view of it: its tree field,
	/// the leaf block of this tree's node state. So the walk and the moves read one tree and
	/// nothing else. See ADR-0005. The leaves themselves are drawn with the field, as one
	/// eight-state block, before this runs.
	///
	/// No random number is drawn inside a parallel region. Every proposal and every decision is
	/// taken in serial code, in clique order. `coretools::instances::randomGenerator()` is
	/// `thread_local`, and a worker thread seeds its copy from the wall clock. So a draw inside a
	/// region would make the chain depend on the schedule and not on the seed. The node-state walk
	/// draws from cell uniforms instead (ADR-0007).
	template<bool FixZ> void update_Z_and_nus_and_alphas_and_branch_lengths(size_t iteration) {
		const size_t n = n_cliques();

		// The stream this tree's node state draws from this iteration. Each tree names its own
		// dimension, so the two trees never share a uniform.
		const TCellUniforms node_state_uniforms(run_seed(), TCellStream::node_state, iteration,
		                                        _dimension);

		// --- the node-state walk and the alpha move ---
		TStattoolsTreeMoves moves(*_alpha_c, *_log_nu_c, n);
		const std::vector<double> LL_after_alpha =
		    tree_update::run<FixZ>(_topology(), _cliques(), _Z, _transition_grids, _grid(), _bins(),
		                           node_state_uniforms, moves);

		// --- the nu move: the candidate likelihood ---
		// The ratio divides by the likelihood each clique ended alpha on, so this scores one grid.
		auto view_factory = [this](size_t c) { return _clique_view(c); };
		_log_nu_c->propose(coretools::TRange(0, n, 1));
		auto nu_kernel = [this](size_t c, TNodeStateCliqueView &states) {
			return _clique_LL(states, _nu_candidate_grid(c, _log_nu_c->value(c)));
		};
		auto [nu_LL_new, nu_indices_to_insert] = clique_traversal::run(n, view_factory, nu_kernel);

		// This region is read-only: it opens a clique view only to score the candidate grid, and
		// writes nothing. Checked once here, against the whole collected result, rather than once
		// per clique inside the loop.
		DEBUG_ASSERT(
		    std::all_of(nu_indices_to_insert.begin(), nu_indices_to_insert.end(),
			            [](const auto &clique_inserts) { return clique_inserts.empty(); }));

		for (size_t i = 0; i < n; ++i) {
			const double proposed = _log_nu_c->value(i);
			const double logH = nu_LL_new[i] - LL_after_alpha[i] + _log_nu_c->getLogDensityRatio(i);
			if (_log_nu_c->acceptOrReject(logH, coretools::TRange(i))) {
				_set_transition_grid(i, _nu_candidate_grid(i, proposed));
			}
		}

		// --- the branch-length move: proposed only now ---
		// Every move above reads the bin each branch holds now. A proposal sets `value` to the
		// candidate, so the branch lengths are proposed after those moves and not before. Only the
		// ratio below reads `oldValue`, and only for the branches it pairs.
		auto pairs         = _build_pairs_branch_lengths();
		const auto n_pairs = pairs.length();
		_propose_new_branch_lengths(pairs);

		// --- region 3: the branch-length likelihood, summed in fixed slots ---
		// A clique belongs to slot `i / cliques_per_slot`, and `schedule(static, cliques_per_slot)`
		// hands each of those blocks to one thread, so a slot has one writer and its contents do
		// not depend on the team size.
		const size_t cliques_per_slot = (n + N_BRANCH_LL_SLOTS - 1) / N_BRANCH_LL_SLOTS;
		std::vector<std::vector<coretools::TSumLogProbability>> branch_LL(
		    N_BRANCH_LL_SLOTS, std::vector<coretools::TSumLogProbability>(n_pairs));
#pragma omp parallel for num_threads(ProgramOptions::NUMBER_OF_THREADS) default(none)              \
    schedule(static, cliques_per_slot) shared(pairs, branch_LL, n, cliques_per_slot)
		for (size_t i = 0; i < n; ++i) {
			auto states = _clique_view(i);
			_add_to_LL_branch_lengths(i, states, branch_LL[i / cliques_per_slot], pairs);
			DEBUG_ASSERT(states.take_deferred_inserts().empty());
		}

		auto log_sum = TTree::_reduce_branch_LL_slots(branch_LL, n_pairs);
		_evalute_update_branch_length(log_sum, pairs);
	}

	[[nodiscard]] TypeBinnedBranchLengths get_binned_branch_length(size_t index_in_tree) const {
		return _binned_branch_lengths->value(_topology().branch_index(index_in_tree));
	}

	[[nodiscard]] const std::string &get_tree_name() const { return _tree_name; }

	/// The branch length each bin stands for. A stored branch is only ever a bin index, so this is
	/// what turns one back into a length.
	[[nodiscard]] const std::vector<double> &grid_branch_lengths() const {
		return _grid().grid_branch_lengths();
	}

	void simulate_Z();

	/// `log p(Z | theta)` for this tree: every clique's own column of the node state, scored under
	/// that clique's process.
	///
	/// One term per node -- a root against the stationary distribution, everything else against its
	/// parent -- so each branch is counted exactly once and the result is a density. That is what
	/// makes it the no-drift instrument ADR-0002 says the old sum of tree likelihoods was not.
	///
	/// It reads the states the node state holds now and the bins the branches sit in now, so it is
	/// the density of the configuration as the iteration leaves it. It is a pass over the whole
	/// node state, which is why writing the joint density is behind a command-line flag.
	[[nodiscard]] double log_node_state_density() {
		// The kernel returns one double per clique, and the helper gathers them into a vector in
		// clique order (clique_traversal::run). A thread-indexed accumulator would add its
		// cliques in whatever order the schedule handed them out, so the sum's rounding would
		// move with the thread count and the trace would stop being reproducible from its seed.
		auto view_factory = [this](size_t c) { return _clique_view(c); };
		auto kernel       = [this](size_t c, TNodeStateCliqueView &states) {
			return node_state_density::log_density_of_clique(_topology(), transition_grid(c),
			                                                 states, _bins());
		};
		auto [per_clique, inserts] = clique_traversal::run(n_cliques(), view_factory, kernel);

		// This pass is read-only, so every view defers nothing. Checked once here, against the
		// whole collected result, rather than once per clique inside the loop.
		DEBUG_ASSERT(std::all_of(inserts.begin(), inserts.end(), [](const auto &clique_inserts) {
			return clique_inserts.empty();
		}));

		// In clique order, so the answer does not depend on how the cliques were shared out.
		return coretools::containerSum(per_clique);
	}

	/// Starts every internal node of this tree at the state its children make most likely. One
	/// forward pass over each clique. This is initialisation and not a sampler move. It runs once,
	/// before the chain's first update, and reads the tree field the chain start wrote
	/// (TMarkovField::_start_the_chain). A node state the run supplied stays as it is.
	void initialize_Z_from_children() {
		std::string set_Z_cli_command = "set_" + get_tree_name() + "_Z";
		if (coretools::instances::parameters().exists(set_Z_cli_command)) { return; }

		auto view_factory = [this](size_t c) { return _clique_view(c); };
		auto kernel       = [this](size_t c, TNodeStateCliqueView &states) {
			_initialize_clique_from_children(c, states);
		};
		const auto indices_to_insert = clique_traversal::run(n_cliques(), view_factory, kernel);

		_Z.insert_in_Z(indices_to_insert);
	};
};
#endif // METABOLITE_INFERENCE_TREE_H
