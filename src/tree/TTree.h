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
	std::vector<double> _nu_c;

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

	/// The bin every branch of this tree sat in before the current round of proposals. Branch
	/// lengths are proposed before the loop over cliques starts, so `value` inside that loop is
	/// already the candidate. The node-state walk and the bottom-up start are handed this, so
	/// neither reads a parameter.
	[[nodiscard]] auto _previous_bins() const {
		return [this](size_t node) { return get_previous_binned_branch_length(node); };
	}

	void _simulateUnderPrior(Storage *) override;

	/// One node's contribution to a clique's log-likelihood under `process`: scored against its
	/// parent, or against the stationary distribution if it is a root.
	///
	/// `branch_len_bin` is the bin the branch sat in *before* this iteration's proposals, because
	/// branch lengths are proposed before the cliques are walked and `value` is already the
	/// candidate by the time this runs.
	template<TransitionGridLike Process>
	void _add_node_to_clique_LL(size_t index_in_tree, bool state_of_node,
	                            coretools::TSumLogProbability &LL,
	                            const TNodeStateCliqueView &states,
	                            std::optional<size_t> branch_len_bin,
	                            const Process &process) const {
		if (_topology().is_root(index_in_tree)) {
			LL.add(process.stationary(state_of_node));
		} else {
			LL.add(_prob_to_parent(index_in_tree, branch_len_bin.value(), states, process));
		}
	}

	/// One clique's log-likelihood under each of two processes, in one pass over its nodes.
	///
	/// Two processes rather than one because a Metropolis move on alpha or nu needs the current
	/// grid and the candidate scored over the same states, and a node's state is the expensive
	/// thing to read: one pass reads it once and the compiler shares the load between both scores.
	template<TransitionGridLike ProcessA, TransitionGridLike ProcessB>
	[[nodiscard]] std::array<double, 2> _clique_LL(const TNodeStateCliqueView &states,
	                                              const ProcessA &a, const ProcessB &b) const {
		coretools::TSumLogProbability LL_a;
		coretools::TSumLogProbability LL_b;
		const auto &topology = _topology();
		for (size_t i = 0; i < topology.n_nodes(); ++i) {
			const bool state_of_node = states.is_one(i);
			std::optional<size_t> branch_len_bin;
			if (!topology.is_root(i)) { branch_len_bin = get_previous_binned_branch_length(i); }
			_add_node_to_clique_LL(i, state_of_node, LL_a, states, branch_len_bin, a);
			_add_node_to_clique_LL(i, state_of_node, LL_b, states, branch_len_bin, b);
		}
		return {LL_a.getSum(), LL_b.getSum()};
	}

	/// One clique's log-likelihood under one process. The nu move needs only this: the likelihood
	/// its ratio divides by is the one the alpha move already computed -- after the alpha move a
	/// clique's grid is either the one alpha started from or the candidate it accepted, and both
	/// of those sums were summed over these same nodes in this same order, so reading them back is
	/// the same double rather than a near one.
	template<TransitionGridLike Process>
	[[nodiscard]] double _clique_LL(const TNodeStateCliqueView &states,
	                                const Process &process) const {
		coretools::TSumLogProbability LL;
		const auto &topology = _topology();
		for (size_t i = 0; i < topology.n_nodes(); ++i) {
			const bool state_of_node = states.is_one(i);
			std::optional<size_t> branch_len_bin;
			if (!topology.is_root(i)) { branch_len_bin = get_previous_binned_branch_length(i); }
			_add_node_to_clique_LL(i, state_of_node, LL, states, branch_len_bin, process);
		}
		return LL.getSum();
	}

	/// The grid a proposed alpha or nu would give clique `c`. Built from the value the parameter
	/// currently holds, so it has to be asked before `acceptOrReject` puts a rejected value back.
	template<bool IsAlpha> [[nodiscard]] TTransitionGrid _candidate_grid(size_t c,
	                                                                     double proposed) const {
		if constexpr (IsAlpha) {
			return TTransitionGrid(proposed, _nu_c[c], _grid());
		} else {
			return TTransitionGrid(_alpha_c->value(c), std::exp(proposed), _grid());
		}
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

	/// The field is not named here at all. Each tree owns a leaf-level view of it: its tree field,
	/// the leaf block of this tree's node state. So the walk below, and the alpha and nu moves
	/// after it, read one tree and nothing else. See ADR-0005. The leaves themselves are drawn
	/// with the field, as one eight-state block, before this runs.
	///
	/// Every read and write of a cell goes through the clique's view, which is the one place a
	/// node index becomes a cell of the node state.
	///
	/// Not one random number is drawn inside a parallel region. Every proposal and every accept
	/// or reject is taken here, single-threaded, walking the cliques in order; the regions only
	/// read states and add up likelihoods. That is what makes the chain a function of the seed:
	/// `coretools::instances::randomGenerator()` is `thread_local`, and a worker thread's copy is
	/// constructed inside the region, where `setSeed` seeds it from the wall clock -- so
	/// `--fixedSeed` never reached the alpha and nu moves at all, and which clique got which draw
	/// depended on how `schedule(dynamic)` happened to hand the cliques out that run. ADR-0007
	/// bought the node state this property by hashing a cell's uniform from its position; these
	/// moves get it by drawing where there is only one thread to draw on.
	///
	/// The cost is three regions where there was one, and two serial passes over the cliques that
	/// do no likelihood work. A proposal and an accept are a few tens of nanoseconds against the
	/// tens of milliseconds a clique's likelihood costs.
	template<bool FixZ> void update_Z_and_nus_and_alphas_and_branch_lengths(size_t iteration) {
		const size_t n = n_cliques();
		std::vector<std::vector<size_t>> indices_to_insert(n);

		// The stream this tree's node state draws from this iteration, built before the parallel
		// region (see run_seed). Each tree names its own dimension, so the two never share a
		// uniform.
		const TCellUniforms node_state_uniforms(run_seed(), TCellStream::node_state, iteration,
		                                        _dimension);

		// build pairs of branch lengths to update, and propose them: both draw, both serial
		auto pairs         = _build_pairs_branch_lengths();
		const auto n_pairs = pairs.length();
		_propose_new_branch_lengths(pairs);

		// Every clique's alpha, proposed here rather than inside the region.
		_alpha_c->propose(coretools::TRange(0, n, 1));

		// --- region 1: the node-state walk, and alpha's two likelihoods ---
		std::vector<std::array<double, 2>> alpha_LL(n);
#pragma omp parallel for num_threads(ProgramOptions::NUMBER_OF_THREADS) default(none)              \
    schedule(dynamic) shared(indices_to_insert, node_state_uniforms, alpha_LL, n)
		for (size_t i = 0; i < n; ++i) {
			// The cells this clique reads and writes. The view lives across the moves below,
			// because those moves read the states the walk assigns.
			auto states = _clique_view(i);
			if constexpr (!FixZ) {
				node_state_walk::update_clique(_topology(), transition_grid(i), _previous_bins(),
				                               node_state_uniforms, states);
			}
			alpha_LL[i] = _clique_LL(states, transition_grid(i),
			                         _candidate_grid<true>(i, _alpha_c->value(i)));
			// The view ends here, inside the parallel region, so it hands its inserts out rather
			// than making them. The list is taken whether or not the walk ran.
			indices_to_insert[i] = states.take_deferred_inserts();
		}
		if constexpr (!FixZ) { _Z.insert_in_Z(indices_to_insert); }

		// --- alpha's decision, and the likelihood the nu move inherits ---
		// Whichever of the two alpha ends on is the likelihood of the clique as nu finds it, so
		// the nu move below scores one grid instead of two.
		std::vector<double> LL_after_alpha(n);
		for (size_t i = 0; i < n; ++i) {
			const double proposed = _alpha_c->value(i);
			const double logH     = alpha_LL[i][1] - alpha_LL[i][0] + _alpha_c->getLogDensityRatio(i);
			if (_alpha_c->acceptOrReject(logH, coretools::TRange(i))) {
				_set_transition_grid(i, _candidate_grid<true>(i, proposed));
				LL_after_alpha[i] = alpha_LL[i][1];
			} else {
				LL_after_alpha[i] = alpha_LL[i][0];
			}
		}

		// --- region 2: nu's candidate likelihood ---
		_log_nu_c->propose(coretools::TRange(0, n, 1));
		std::vector<double> nu_LL_new(n);
#pragma omp parallel for num_threads(ProgramOptions::NUMBER_OF_THREADS) default(none)              \
    schedule(dynamic) shared(nu_LL_new, n)
		for (size_t i = 0; i < n; ++i) {
			auto states  = _clique_view(i);
			nu_LL_new[i] = _clique_LL(states, _candidate_grid<false>(i, _log_nu_c->value(i)));
			DEBUG_ASSERT(states.take_deferred_inserts().empty());
		}

		for (size_t i = 0; i < n; ++i) {
			const double proposed = _log_nu_c->value(i);
			const double logH = nu_LL_new[i] - LL_after_alpha[i] + _log_nu_c->getLogDensityRatio(i);
			if (_log_nu_c->acceptOrReject(logH, coretools::TRange(i))) {
				_set_transition_grid(i, _candidate_grid<false>(i, proposed));
				_nu_c[i] = std::exp(proposed);
			}
		}

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

	/// The bin this branch had before the current round of proposals. The clique update needs it
	/// because branch lengths are proposed before the loop over cliques starts, so `value` inside
	/// the loop is already the candidate. Exposed here so that a clique can ask the tree rather
	/// than be handed a pointer to the parameter itself.
	[[nodiscard]] TypeBinnedBranchLengths
	get_previous_binned_branch_length(size_t index_in_tree) const {
		return _binned_branch_lengths->oldValue(_topology().branch_index(index_in_tree));
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
		auto kernel = [this](size_t c, TNodeStateCliqueView &states) {
			return node_state_density::log_density_of_clique(
			    _topology(), transition_grid(c), states,
			    [this](size_t node) { return get_binned_branch_length(node); });
		};
		auto [per_clique, inserts] = clique_traversal::run(n_cliques(), view_factory, kernel);

		// This pass is read-only, so every view defers nothing. Checked once here, against the
		// whole collected result, rather than once per clique inside the loop.
		DEBUG_ASSERT(std::all_of(inserts.begin(), inserts.end(),
		                         [](const auto &clique_inserts) { return clique_inserts.empty(); }));

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
		auto kernel = [this](size_t c, TNodeStateCliqueView &states) {
			_initialize_clique_from_children(c, states);
		};
		const auto indices_to_insert = clique_traversal::run(n_cliques(), view_factory, kernel);

		_Z.insert_in_Z(indices_to_insert);
	};
};
#endif // METABOLITE_INFERENCE_TREE_H
