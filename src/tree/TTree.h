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
#include "tree/clique/TCliqueView.h"
#include "tree/node_state_density.h"
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
	std::vector<std::optional<TTransitionGrid>> _transition_grids;
	IndexArray _dimension_cliques;
	std::vector<std::string> _clique_names;

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
	[[nodiscard]] double
	_calculate_likelihood_ratio_branch_length(size_t index_in_binned_branch_length,
	                                          const TTransitionGrid &process,
	                                          const TNodeStateCliqueView &states) const;

	/// The multidimensional index of clique `c`: a leaf in every dimension but this tree's own,
	/// which carries a 0.
	///
	/// It reads `c` as a row-major subscript, and transition_grid_of_cell walks a column-major
	/// stride. The two agree because clique space has one dimension above 1, this tree's own
	/// carrying a 1 and there being two trees. A third tree would need one convention here.
	[[nodiscard]] IndexArray _clique_index(size_t c) const;

	/// Installs the process of clique `c`. Called once the parameters exist, and again whenever a
	/// proposal on alpha or nu is accepted. Nothing outside the tree installs a grid.
	void _set_transition_grid(size_t c, TTransitionGrid grid) {
		_transition_grids[c] = std::move(grid);
	}

	/// P(node | parent) under an explicitly given process, so a Metropolis proposal can ask the
	/// same question of the clique's current grid and of its candidate.
	[[nodiscard]] double _prob_to_parent(size_t index_in_tree,
	                                     TypeBinnedBranchLengths binned_branch_length,
	                                     const TNodeStateCliqueView &states,
	                                     const TTransitionGrid &process) const {
		const size_t parent_index = _topology().parent_of(index_in_tree);
		const bool parent_state   = states.is_one(parent_index);
		const bool child_state    = states.is_one(index_in_tree);
		return process.probability(binned_branch_length, parent_state, child_state);
	}

	/// The node-state walk over one clique, and the bottom-up start it shares its arithmetic with.
	void _update_Z_of_clique(size_t c, TNodeStateCliqueView &states,
	                         const TCellUniforms &uniforms) const;
	void _initialize_clique_from_children(size_t c, TNodeStateCliqueView &states) const;
	void _initialize_node_from_children(size_t node_index, const TTransitionGrid &process,
	                                    TNodeStateCliqueView &states) const;
	static void _log_prob_root(double stationary_0,
	                           std::array<coretools::TSumLogProbability, 2> &sum_log);
	void _log_prob_node_to_children(size_t index_in_tree, const TTransitionGrid &process,
	                                const TNodeStateCliqueView &states,
	                                std::array<coretools::TSumLogProbability, 2> &sum_log) const;

	void _simulateUnderPrior(Storage *) override;

	/// One node's contribution to a clique's log-likelihood under `process`. Called twice per
	/// node, once with the clique's current grid and once with the proposal's candidate.
	void _compute_LL_old_and_new_nu_or_alpha(size_t index_in_tree, bool state_of_node,
	                                         coretools::TSumLogProbability &LL,
	                                         const TNodeStateCliqueView &states,
	                                         std::optional<size_t> branch_len_bin,
	                                         const TTransitionGrid &process) const {
		if (_topology().is_root(index_in_tree)) {
			LL.add(process.stationary(state_of_node));
		} else {
			LL.add(_prob_to_parent(index_in_tree, branch_len_bin.value(), states, process));
		}
	}

	template<bool IsAlpha, typename TypeParam>
	void _update_nu_or_alpha(const TNodeStateCliqueView &states, size_t c, TypeParam *param) {
		// propose a new value
		param->propose(coretools::TRange(c));

		double new_value;
		if constexpr (IsAlpha) {
			new_value = param->value(c);
		} else {
			new_value = std::exp(param->value(c));
		}

		// No need to mutate anything: the candidate is a second grid built from the proposed value,
		// and the clique keeps whichever of the two is accepted. The old value is not read back
		// from the parameter either -- the clique's current grid still carries it.
		const TTransitionGrid &current  = transition_grid(c);
		const TTransitionGrid candidate = [&] {
			if constexpr (IsAlpha) {
				return TTransitionGrid(new_value, _nu_c[c], _grid());
			} else {
				return TTransitionGrid(_alpha_c->value(c), new_value, _grid());
			}
		}();

		coretools::TSumLogProbability LL_old;
		coretools::TSumLogProbability LL_new;
		const auto &topology = _topology();
		for (size_t i = 0; i < topology.n_nodes(); ++i) {
			bool state_of_node = states.is_one(i);

			// Note: need to take oldValue because we update _binned_branch_length before
			// starting the loop!!!
			std::optional<size_t> branch_len_bin;
			if (!topology.is_root(i)) { branch_len_bin = get_previous_binned_branch_length(i); }

			_compute_LL_old_and_new_nu_or_alpha(i, state_of_node, LL_old, states, branch_len_bin,
			                                    current);
			_compute_LL_old_and_new_nu_or_alpha(i, state_of_node, LL_new, states, branch_len_bin,
			                                    candidate);
		}

		// calculate Hastings ratio
		const double LLRatio       = LL_new.getSum() - LL_old.getSum();
		const double logPriorRatio = param->getLogDensityRatio(c);
		const double logH          = LLRatio + logPriorRatio;

		// accept or reject
		bool accepted = param->acceptOrReject(logH, coretools::TRange(c));
		if (accepted) {
			_set_transition_grid(c, candidate);
			if constexpr (!IsAlpha) { _nu_c[c] = new_value; }
		}
	}

	/// The cells of clique `c`, addressed by node index. One place holds the three things a view
	/// is built from, so the loops below cannot drift apart.
	[[nodiscard]] TNodeStateCliqueView _clique_view(size_t c) {
		return {_Z, _topology(), _clique_index(c), _dimension};
	}

	void _evalute_update_branch_length(std::vector<coretools::TSumLogProbability> &log_sum,
	                                   const stattools::TPairIndexSampler &pairs);

	/// @brief Helper function to reduce the parallelized log_sum into a log_sum
	static std::vector<coretools::TSumLogProbability> _reduce_log_sum_per_thread(
	    std::vector<std::vector<coretools::TSumLogProbability>> &log_sum_per_thread,
	    size_t n_pairs) {
		auto &log_sum_b = log_sum_per_thread[0];
		for (size_t t = 1; t < ProgramOptions::NUMBER_OF_THREADS; ++t) {
			for (size_t p = 0; p < n_pairs; ++p) {
				log_sum_b[p] = log_sum_b[p] + log_sum_per_thread[t][p];
			}
		}
		return log_sum_b;
	};

public:
	TTree(size_t dimension, const std::string &filename, const std::string &tree_name,
	      TypeParamAlpha *Alpha, TypeParamLogNu *LogNu,
	      TypeParamBinBranches *Binned_Branch_Lenghts);
	~TTree() override;

	[[nodiscard]] size_t size() const { return _topology().n_nodes(); };

	/// The topology this tree is built on. Everything that does not need the parameters -- the
	/// current state, the sheet, the clique's walks -- should ask this rather than the tree.
	[[nodiscard]] const TPhylogeny &phylogeny() const { return _topology(); }

	[[nodiscard]] size_t parent_of(size_t index) const { return _topology().parent_of(index); }
	[[nodiscard]] std::span<const size_t> children_of(size_t index) const {
		return _topology().children_of(index);
	}
	[[nodiscard]] bool is_root(size_t index) const { return _topology().is_root(index); }
	[[nodiscard]] bool isLeaf(size_t index) const { return _topology().is_leaf(index); }

	/** Get the index of a node by its id
	 * @param Id: the id of the node
	 * @return the index of the node with the given id
	 */
	[[nodiscard]] size_t get_node_index(const std::string &Id) const {
		return _topology().index_of(Id);
	}

	/** @return the number of leaves in the tree
	 */
	[[nodiscard]] size_t get_number_of_leaves() const { return _topology().n_leaves(); }
	[[nodiscard]] size_t get_number_of_nodes() const { return _topology().n_nodes(); }
	[[nodiscard]] size_t get_number_of_roots() const { return _topology().n_roots(); }

	/** @param node_index: the index of the node within the tree
	 * @return The index of the node in leaf space. Meaningless if the node is not a leaf; the node
	 * index alone says whether it is one.
	 */
	[[nodiscard]] size_t get_index_within_leaves(size_t node_index) const {
		return _topology().leaf_index(node_index);
	}
	[[nodiscard]] size_t get_index_within_leaves(const std::string &node_name) const {
		return _topology().leaf_index(get_node_index(node_name));
	}
	[[nodiscard]] size_t get_node_index_from_leaf_index(size_t leaf_index) const {
		return _topology().leaves()[leaf_index];
	}

	/** @return The root nodes of the tree, as the range of node indices they occupy
	 */
	[[nodiscard]] auto get_root_nodes() const { return _topology().roots(); }
	[[nodiscard]] auto get_internal_nodes() const { return _topology().internal_nodes(); }

	/** Checks whether a node is in the tree
	 * @param node_id: the id of the node
	 * @return true if the node is in the tree, false otherwise
	 */
	[[nodiscard]] bool in_tree(const std::string &node_id) const {
		return _topology().contains(node_id);
	};

	// stattools stuff
	[[nodiscard]] std::string name() const override;
	void initialize() override;
	[[nodiscard]] double getSumLogPriorDensity(const Storage &) const override;
	void guessInitialValues() override;
	[[nodiscard]] double getDensity(const Storage &, size_t) const override;
	[[nodiscard]] double getLogDensityRatio(const UpdatedStorage &, size_t) const override;

	void initialize_cliques_and_Z(const std::vector<std::unique_ptr<TTree>> &all_trees);

	/// The number of cliques of this tree, which is the number of transition grids it holds.
	[[nodiscard]] size_t n_cliques() const { return _transition_grids.size(); }

	/// The process of clique `c`.
	[[nodiscard]] const TTransitionGrid &transition_grid(size_t c) const {
		return _transition_grids[c].value();
	}

	/// The process of the clique a cell belongs to. A clique of this tree is named by a leaf of
	/// every other tree, so the cell's own dimension is dropped on the way in.
	[[nodiscard]] const TTransitionGrid &
	transition_grid_of_cell(const IndexArray &index_in_leaves_space) const;
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
	template<bool FixZ> void update_Z_and_nus_and_alphas_and_branch_lengths(size_t iteration) {
		std::vector<std::vector<size_t>> indices_to_insert(n_cliques());

		// The stream this tree's node state draws from this iteration, built before the parallel
		// region (see run_seed). Each tree names its own dimension, so the two never share a
		// uniform.
		const TCellUniforms node_state_uniforms(run_seed(), TCellStream::node_state, iteration,
		                                        _dimension);

		// build pairs of branch lengths to update
		auto pairs         = _build_pairs_branch_lengths();
		const auto n_pairs = pairs.length();
		std::vector<std::vector<coretools::TSumLogProbability>> log_sum_per_thread(
		    ProgramOptions::NUMBER_OF_THREADS, std::vector<coretools::TSumLogProbability>(n_pairs));

		// propose new branch lengths
		_propose_new_branch_lengths(pairs);

#pragma omp parallel for num_threads(ProgramOptions::NUMBER_OF_THREADS) default(none)              \
    schedule(dynamic) shared(pairs, log_sum_per_thread, indices_to_insert, node_state_uniforms)
		for (size_t i = 0; i < n_cliques(); ++i) {
			auto &log_sum_local = log_sum_per_thread[omp_get_thread_num()];
			// The cells this clique reads and writes. The view lives across the moves below,
			// because those moves read the states the walk assigns.
			auto states         = _clique_view(i);
			// update Z
			if constexpr (!FixZ) { _update_Z_of_clique(i, states, node_state_uniforms); }

			// update nu and alpha
			_update_nu_or_alpha<true>(states, i, _alpha_c);
			_update_nu_or_alpha<false>(states, i, _log_nu_c);

			// add to likelihood ratio for branch length
			_add_to_LL_branch_lengths(i, states, log_sum_local, pairs);

			// The view ends here, inside the parallel region, so it hands its inserts out rather
			// than making them. The list is taken whether or not the walk ran.
			indices_to_insert[i] = states.take_deferred_inserts();
		}

		// update branch lengths
		auto log_sum_b = TTree::_reduce_log_sum_per_thread(log_sum_per_thread, n_pairs);
		_evalute_update_branch_length(log_sum_b, pairs);

		if constexpr (!FixZ) { _Z.insert_in_Z(indices_to_insert); }
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
		// One slot per clique, and not one per thread. A thread-indexed accumulator adds its
		// cliques in whatever order the schedule handed them out, so the sum's rounding would move
		// with the thread count and the trace would stop being reproducible from its seed.
		std::vector<double> per_clique(n_cliques(), 0.0);

#pragma omp parallel for num_threads(ProgramOptions::NUMBER_OF_THREADS)                            \
    schedule(dynamic) default(none) shared(per_clique)
		for (size_t i = 0; i < n_cliques(); ++i) {
			auto states   = _clique_view(i);
			per_clique[i] = node_state_density::log_density_of_clique(
			    _topology(), transition_grid(i), states,
			    [this](size_t node) { return get_binned_branch_length(node); });
			// The view is read and never written, so it defers nothing. Its list is taken here
			// all the same, so that no view reaches its storage from inside the parallel region.
			const std::vector<size_t> inserts = states.take_deferred_inserts();
			DEBUG_ASSERT(inserts.empty());
		}
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

		// Each clique is independent of each other so we should be able to parallelize this
		std::vector<std::vector<size_t>> indices_to_insert(n_cliques());

#pragma omp parallel for num_threads(ProgramOptions::NUMBER_OF_THREADS)                            \
    schedule(dynamic) default(none) shared(indices_to_insert)
		for (size_t i = 0; i < n_cliques(); ++i) {
			auto states = _clique_view(i);
			_initialize_clique_from_children(i, states);
			indices_to_insert[i] = states.take_deferred_inserts();
		}

		_Z.insert_in_Z(indices_to_insert);
	};
};
#endif // METABOLITE_INFERENCE_TREE_H
