//
// Created by Marco Visani on 26.06.23.
//

#pragma once

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
#include "stattools/ParametersObservations/TParameter.h"
#include "storages/storage_backend.h"
#include "tree/TPhylogeny.h"
#include "tree/branch/TBinGrid.h"
#include "tree/branch/TTransitionGrid.h"
#include <cstddef>
#include <optional>
#include <span>
#include <string>
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

	// cliques
	std::vector<TTransitionGrid> _transition_grid_per_clique;
	size_t _n_cliques;
	IndexArray _dimension_cliques;
	std::vector<std::string> _clique_names;

	// Nus
	TypeParamLogNu *_log_nu_c = nullptr;
	std::vector<double> _nu_c;

	// Alphas
	TypeParamAlpha *_alpha_c = nullptr;

	// Set Z
	TTreeStateStorage _Z;

	// Joint probability density
	std::vector<double> _joint_log_prob_density;

	// private functions
	void _reset_joint_log_prob_density() {
		_joint_log_prob_density.clear();
		_joint_log_prob_density.resize(ProgramOptions::NUMBER_OF_THREADS);
	}
	void _set_initial_branch_lengths(bool is_simulation);
	[[nodiscard]] std::vector<size_t>
	_bin_branch_lengths(const std::vector<double> &branch_lengths) const;
	void _bin_branch_lengths_from_tree();
	void _initialize_grid_branch_lengths();
	void _initialize_Z(IndexArray num_leaves_per_tree,
	                   [[maybe_unused]] const std::vector<std::unique_ptr<TTree>> &all_trees);
	void _initialize_cliques(const IndexArray &num_leaves_per_tree,
	                         const std::vector<std::unique_ptr<TTree>> &all_trees);
	/// @brief Load tree from file
	void _load_from_file(const std::string &filename, const std::string &tree_name);

	// updating branch lengths
	[[nodiscard]] stattools::TPairIndexSampler _build_pairs_branch_lengths() const;
	void _propose_new_branch_lengths(const stattools::TPairIndexSampler &pairs);
	void _propose_new_branch_lengths(size_t p1, size_t p2, int val);
	void _add_to_LL_branch_lengths(size_t clique_index,
	                               std::vector<coretools::TSumLogProbability> &log_sum,
	                               const stattools::TPairIndexSampler &pairs) const;
	[[nodiscard]] double
	_calculate_likelihood_ratio_branch_length(size_t index_in_binned_branch_length,
	                                          size_t clique_index) const;

	void _simulateUnderPrior(Storage *) override;

	/// One node's contribution to a clique's log-likelihood under `process`. Called twice per
	/// node, once with the clique's current grid and once with the proposal's candidate.
	void _compute_LL_old_and_new_nu_or_alpha(size_t index_in_tree, IndexArray index,
	                                         bool state_of_node, coretools::TSumLogProbability &LL,
	                                         std::optional<size_t> branch_len_bin,
	                                         const TTransitionGrid &process) const {
		if (_topology().is_root(index_in_tree)) {
			LL.add(process.stationary(state_of_node));
		} else {
			double prob = _calculate_prob_to_parent(index_in_tree, index, state_of_node,
			                                        branch_len_bin.value(), process);
			LL.add(prob);
		}
	}

	[[nodiscard]] double _calculate_prob_to_parent(size_t index_in_tree, IndexArray index,
	                                               bool node_state,
	                                               TypeBinnedBranchLengths binned_branch_length,
	                                               const TTransitionGrid &process) const {
		size_t parent_index              = parent_of(index_in_tree);
		IndexArray parent_index_in_tree  = index;
		parent_index_in_tree[_dimension] = parent_index;
		bool parent_state                = _Z.is_one(parent_index_in_tree).is_one;
		return process.probability(binned_branch_length, parent_state, node_state);
	}

	template<bool IsAlpha, typename TypeParam>
	void _update_nu_or_alpha(size_t clique_index, TypeParam *param) {
		// propose a new value
		param->propose(coretools::TRange(clique_index));

		double new_value;
		if constexpr (IsAlpha) {
			new_value = param->value(clique_index);
		} else {
			new_value = std::exp(param->value(clique_index));
		}

		// No need to mutate anything: the candidate is a second grid built from the proposed value,
		// and the clique keeps whichever of the two is accepted. The old value is not read back
		// from the parameter either -- the clique's current grid still carries it.
		const auto &current             = _transition_grid_per_clique[clique_index];
		const TTransitionGrid candidate = [&] {
			if constexpr (IsAlpha) {
				return TTransitionGrid(new_value, _nu_c[clique_index], _grid());
			} else {
				return TTransitionGrid(_alpha_c->value(clique_index), new_value, _grid());
			}
		}();

		coretools::TSumLogProbability LL_old;
		coretools::TSumLogProbability LL_new;
		const auto &topology = _topology();
		IndexArray multi_dim_index =
		    coretools::getSubscriptsAsArray(clique_index, _dimension_cliques);
		for (size_t i = 0; i < topology.n_nodes(); ++i) {
			// the _dimension tells along which axis we are going so multi_dim_index[_dimension] =
			// i; the other dimension is the clique_index
			multi_dim_index[_dimension] = i;
			bool state_of_node          = _Z.is_one(multi_dim_index).is_one;

			// Note: need to take oldValue because we update _binned_branch_length before
			// starting the loop!!!
			std::optional<size_t> branch_len_bin;
			if (!topology.is_root(i)) { branch_len_bin = get_previous_binned_branch_length(i); }

			_compute_LL_old_and_new_nu_or_alpha(i, multi_dim_index, state_of_node, LL_old,
			                                    branch_len_bin, current);
			_compute_LL_old_and_new_nu_or_alpha(i, multi_dim_index, state_of_node, LL_new,
			                                    branch_len_bin, candidate);
		}

		// calculate Hastings ratio
		const double LLRatio       = LL_new.getSum() - LL_old.getSum();
		const double logPriorRatio = param->getLogDensityRatio(clique_index);
		const double logH          = LLRatio + logPriorRatio;

		// accept or reject
		bool accepted = param->acceptOrReject(logH, coretools::TRange(clique_index));
		if (accepted) {
			_transition_grid_per_clique[clique_index] = candidate;
			if constexpr (!IsAlpha) { _nu_c[clique_index] = new_value; }
		}
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

	/// @brief Calculates the log probability of the root given the stationary distribution
	static void _calculate_log_prob_root(double stationary_0,
	                                     std::array<coretools::TSumLogProbability, 2> &sum_log);

	/// @brief Calculates the log probability of a node to its children
	void _calculate_log_prob_node_to_children(
	    size_t index_in_tree, size_t clique_index,
	    std::array<coretools::TSumLogProbability, 2> &sum_log) const;

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
	[[nodiscard]] size_t get_number_of_internal_nodes() const {
		return _topology().n_internal_nodes();
	}
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
	[[nodiscard]] size_t get_node_index_from_internal_nodes_index(size_t internal_index) const {
		return _topology().internal_nodes()[internal_index];
	}

	/** @param node_index: the index of the node within the tree
	 * @return The index of the node in internal-node space. Meaningless if the node is a leaf; the
	 * node index alone says whether it is one.
	 */
	[[nodiscard]] size_t get_index_within_internal_nodes(size_t node_index) const {
		return _topology().internal_index(node_index);
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

	[[nodiscard]] const TTreeStateStorage &get_Z() const;
	TTreeStateStorage &get_Z();

	[[nodiscard]] std::string get_node_id(size_t index) const { return _topology().id_of(index); }

	void update_Z_clique(size_t clique_index, std::vector<double> &joint_prob_density,
	                     const TFieldStorage &Y);
	void _simulation_prepare_cliques(size_t c);

	template<bool IsSimulation, bool FixZ>
	void update_Z_and_nus_and_alphas_and_branch_lengths(const TFieldStorage &Y) {
		_reset_joint_log_prob_density();

		// build pairs of branch lengths to update
		auto pairs         = _build_pairs_branch_lengths();
		const auto n_pairs = pairs.length();
		std::vector<std::vector<coretools::TSumLogProbability>> log_sum_per_thread(
		    ProgramOptions::NUMBER_OF_THREADS, std::vector<coretools::TSumLogProbability>(n_pairs));

		// propose new branch lengths
		if constexpr (!IsSimulation) { _propose_new_branch_lengths(pairs); }

#pragma omp parallel for num_threads(ProgramOptions::NUMBER_OF_THREADS) default(none)              \
    schedule(dynamic) shared(pairs, log_sum_per_thread, Y)
		for (size_t i = 0; i < _transition_grid_per_clique.size(); ++i) {
			auto &log_sum_local = log_sum_per_thread[omp_get_thread_num()];
			// update Z
			if constexpr (!FixZ) { update_Z_clique(i, _joint_log_prob_density, Y); }

			// update nu and alpha
			if constexpr (!IsSimulation) {
				_update_nu_or_alpha<true>(i, _alpha_c);
				_update_nu_or_alpha<false>(i, _log_nu_c);

				// add to likelihood ratio for branch length
				_add_to_LL_branch_lengths(i, log_sum_local, pairs);
			}
		}

		// update branch lengths
		if constexpr (!IsSimulation) {
			auto log_sum_b = TTree::_reduce_log_sum_per_thread(log_sum_per_thread, n_pairs);
			_evalute_update_branch_length(log_sum_b, pairs);
		}
	}

	[[nodiscard]] TypeBinnedBranchLengths get_binned_branch_length(size_t index_in_tree) const {
		return _binned_branch_lengths->value(_topology().branch_index(index_in_tree));
	}

	/// The bin this branch had before the current round of proposals. The clique sweep needs it
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

	[[nodiscard]] double get_complete_joint_density() const {
		return coretools::containerSum(_joint_log_prob_density);
	}

	void initialize_Z_from_children(const TFieldStorage &Y);
	void
	calculate_log_prob_parent_to_node(size_t index_in_tree, size_t clique_index,
	                                  TypeBinnedBranchLengths binned_branch_length,
	                                  std::array<coretools::TSumLogProbability, 2> &sum_log) const;

	[[nodiscard]] size_t clique_index(const IndexArray &index_in_leaves_space) const;
};

inline bool sample(std::array<coretools::TSumLogProbability, 2> &sum_log) {
	const double log_Q = sum_log[1].getSum() - sum_log[0].getSum();
	return coretools::TAcceptOddsRatio::accept(log_Q);
}

inline bool sample(double log_prob_0, double log_prob_1) {
	const double log_Q = log_prob_1 - log_prob_0;
	return coretools::TAcceptOddsRatio::accept(log_Q);
}
