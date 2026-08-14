//
// Created by Marco Visani on 26.06.23.
//

#ifndef METABOLITE_INFERENCE_TREE_H
#define METABOLITE_INFERENCE_TREE_H

#include "TClique.h"
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
#include "storages/y_storage/TStorageYMatrix.h"
#include "storages/z_storage/TStorageZMatrix.h"
#include "tree/node.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
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
	std::vector<TNode> _nodes;                         // a map to store nodes with their ids
	std::unordered_map<std::string, size_t> _node_map; // for fast access to nodes
	std::vector<size_t> _leaves;
	std::vector<size_t> _roots;
	std::vector<size_t> _internal_nodes;
	std::vector<size_t> _internal_nodes_without_roots;
	std::vector<size_t> _leaves_and_internal_nodes_without_roots;

	/// Internal nodes (roots included) ordered so that every node appears after all of its
	/// children. Built once when the tree is read. The block sampler walks it forwards for the
	/// pruning pass and backwards (parents before children) for the sampling pass.
	std::vector<size_t> _internal_nodes_post_order;

	/// 0, 1, ... n_nodes - 1. Lets the pseudo-likelihood loop switch between "internal nodes
	/// only" and "every node" by picking a vector rather than branching per node.
	std::vector<size_t> _all_nodes_indices;

	// The four vectors below have size _nodes.size()
	std::vector<size_t> _leafIndices;
	std::vector<size_t> _rootIndices;
	std::vector<size_t> _internalIndices;
	std::vector<size_t> _leaves_and_internal_nodes_without_roots_indices;
	std::vector<size_t> _internalIndicesWithoutRoots;

	// dimension of the tree
	size_t _dimension;

	// For binning branch lengths
	double _delta          = 0.0;
	size_t _number_of_bins = 0;
	std::vector<double> _grid_branch_lengths;
	std::vector<size_t> _binned_branch_lengths_from_tree;
	TypeParamBinBranches *_binned_branch_lengths = nullptr;

	// cliques
	std::vector<TClique> _cliques;
	IndexArray _dimension_cliques;
	std::vector<std::string> _clique_names;

	// Nus
	TypeParamLogNu *_log_nu_c = nullptr;
	std::vector<double> _nu_c;

	// Alphas
	TypeParamAlpha *_alpha_c = nullptr;

	// All trees, needed by the leaf terms of the pseudo-likelihood: the full conditional of a Y
	// cell couples both trees. Set in initialize_cliques_and_Z, which already receives them.
	// Owned by TModel, which outlives every TTree.
	const std::vector<std::unique_ptr<TTree>> *_all_trees = nullptr;

	/// State of the *other* tree's parent for every Y cell of this tree, laid out for this
	/// tree's access pattern: index [clique * number_of_leaves + index_within_leaves].
	///
	/// For the species tree the underlying quantity is Z_molecules[i, parent(j)] with i the
	/// species leaf and j the molecule leaf; reading that straight out of the other tree's Z
	/// walks it with a stride, so it is transposed into this cache once per sweep instead of
	/// being chased per node. One bit per cell: 3 MB at 5000 x 5000.
	std::vector<bool> _other_parent_state;
	void _refresh_other_parent_state();

	/// P_other(parent state -> s) for s = 0, 1, at the Y cell (clique, index_within_leaves).
	/// This is the factor that stops the leaf conditional from collapsing to this tree's own
	/// transition probability; it does not depend on this tree's parameters but does not cancel,
	/// because it sits inside the normaliser.
	[[nodiscard]] std::array<double, 2> _prob_other_tree_at_leaf(size_t clique,
	                                                             size_t index_within_leaves) const;

	// Set Z
	TStorageZMatrix _Z;

	// Joint probability density
	std::vector<double> _joint_log_prob_density;

	// private functions
	void _reset_joint_log_prob_density() {
		_joint_log_prob_density.clear();
		_joint_log_prob_density.resize(ProgramOptions::NUMBER_OF_THREADS);
	}
	void _set_initial_branch_lengths(bool is_simulation);
	[[nodiscard]] std::vector<size_t> _bin_branch_lengths(const std::vector<double> &branch_lengths,
	                                                      bool exclude_root) const;
	void _bin_branch_lengths_from_tree(std::vector<double> &branch_lengths);
	void _initialize_grid_branch_lengths();
	void _initialize_Z(IndexArray num_leaves_per_tree);
	void _initialize_cliques(const IndexArray &num_leaves_per_tree,
	                         const std::vector<std::unique_ptr<TTree>> &all_trees);
	/// @brief Load tree from file
	void _load_from_file(const std::string &filename, const std::string &tree_name);
	void _add_parent(const std::string &parent_id, std::vector<double> &branch_lengths);
	void _add_child(const std::string &child_id, size_t parent_index, bool is_root,
	                std::vector<double> &branch_lengths, double branch_length_of_child);
	void _simulation_prepare_cliques(size_t c, TClique &clique) const;
	void _simulate_one(const TClique &clique, TCurrentState &current_state, size_t tree_index,
	                   size_t node_index_in_tree);

	// updating branch lengths
	[[nodiscard]] stattools::TPairIndexSampler _build_pairs_branch_lengths() const;
	void _propose_new_branch_lengths(const stattools::TPairIndexSampler &pairs);
	void _propose_new_branch_lengths(size_t p1, size_t p2, int val);
	void _add_to_LL_branch_lengths(size_t c, const TCurrentState &current_state,
	                               std::vector<coretools::TSumLogProbability> &log_sum,
	                               const stattools::TPairIndexSampler &pairs) const;
	[[nodiscard]] double
	_calculate_likelihood_ratio_branch_length(size_t index_in_binned_branch_length,
	                                          const TClique &clique,
	                                          const TCurrentState &current_state) const;

	void _simulateUnderPrior(Storage *) override;

	/// @brief Adds one node's contribution to the pseudo-likelihood of nu / alpha.
	///
	/// The Markov field is a product of the per-tree factors that share the leaf layer Y:
	///   p(Y, Z_species, Z_molecules | theta) = f_species * f_molecules / C(theta).
	/// Scoring a proposal with the raw factor f alone (what this function used to do) drops
	/// C(theta), which depends on nu, alpha and the branch lengths of *both* trees. That makes the
	/// estimator biased, and it leaves the nu -> 0 direction unbounded: as nu shrinks every
	/// P(child | parent) tends to 1, so log f grows without the compensating -log C.
	///
	/// Instead we score with the node's *normalised full conditional*
	///   p(x_i | x_-i) = P(x_i | pa) * prod_children P(child | x_i)
	///                   / sum_{x' in {0,1}} P(x' | pa) * prod_children P(child | x'),
	/// i.e. maximum pseudo-likelihood. Every term is a genuine probability, C(theta) cancels
	/// between numerator and denominator, and no normalising constant has to be evaluated.
	///
	/// Only internal nodes are scored (see _update_nu_or_alpha). Their conditionals involve this
	/// tree's transition matrices only, so no cross-tree information is needed. Leaves are shared
	/// with the other tree and their conditional would need that tree's parent term; leaving them
	/// out keeps the estimator consistent (a composite likelihood over any subset of correct
	/// conditionals is), at the cost of some efficiency, and it also removes the double counting of
	/// the leaf layer that the old code performed (each Y cell was scored once per tree).
	template<bool UseTryMatrix>
	void _add_pseudo_likelihood_term(const TNode &node, size_t index_in_tree, const TClique &clique,
	                                 bool state_of_node, coretools::TSumLogProbability &LL,
	                                 const TCurrentState &current_state,
	                                 std::optional<size_t> branch_len_bin, double alpha,
	                                 size_t clique_index) const {
		std::array<coretools::TSumLogProbability, 2> sum_log;

		// P(node = 0 | parent) and P(node = 1 | parent), or the stationary probabilities at a root
		if (node.is_root()) {
			for (size_t i = 0; i < 2; ++i) {
				sum_log[i].add(TClique::get_stationary_probability(i == 1, alpha));
			}
		} else {
			clique.calculate_log_prob_parent_to_node<TCurrentState, UseTryMatrix>(
			    index_in_tree, branch_len_bin.value(), this, 0, current_state, sum_log);
		}

		if (node.is_leaf()) {
			// A leaf is a cell of Y, which both trees claim. Its full conditional under the field
			// is proportional to the product of the two trees' parent terms, so the other tree's
			// factor has to be carried: it is constant in this tree's parameters, but it sits
			// inside the normaliser below and therefore does not cancel from the ratio.
			//
			// Note the data likelihood does *not* appear. The pseudo-likelihood being maximised
			// is that of the field p(Y, Z | theta); P(data | Y) carries no theta and drops out of
			// the Hastings ratio entirely.
			const auto prob_other =
			    _prob_other_tree_at_leaf(clique_index, _leafIndices[index_in_tree]);
			for (size_t i = 0; i < 2; ++i) { sum_log[i].add(prob_other[i]); }
		} else {
			// P(child | node = 0) and P(child | node = 1) for every child
			clique.calculate_log_prob_node_to_children(
			    index_in_tree, this, current_state, sum_log, _binned_branch_lengths,
			    _leaves_and_internal_nodes_without_roots_indices, UseTryMatrix);
		}

		// normalise over the two states of this node. Done as a logistic of the log difference so
		// that the two products never have to be exponentiated on their own; the clamp keeps a
		// vanishing conditional at ~1e-305 instead of an exact 0, which would make LL -infinity and
		// poison the Hastings ratio.
		const double log_p_state = sum_log[static_cast<size_t>(state_of_node)].getSum();
		const double log_p_other = sum_log[static_cast<size_t>(!state_of_node)].getSum();
		const double diff        = std::min(log_p_other - log_p_state, 700.0);
		LL.add(1.0 / (1.0 + std::exp(diff)));
	}

	template<bool IsAlpha, typename TypeParam>
	void _update_nu_or_alpha(const TCurrentState &current_state, size_t c, TypeParam *param) {
		// propose a new value
		param->propose(coretools::TRange(c));

		// calculate LL for old mu
		// No need to change Lambda (rate matrix), just go over entire tree and calculate
		// probabilities
		double old_value;
		double new_value;
		if constexpr (IsAlpha) {
			old_value = param->oldValue(c);
			new_value = param->value(c);
		} else {
			old_value = _nu_c[c];
			new_value = std::exp(param->value(c));
		}

		coretools::TSumLogProbability LL_old;
		coretools::TSumLogProbability LL_new;
		if constexpr (IsAlpha) {
			_cliques[c].update_lambda(new_value, _nu_c[c]);
		} else {
			_cliques[c].update_lambda(_alpha_c->value(c), new_value);
		}

		const auto &clique = _cliques[c];
		std::optional<size_t> branch_len_bin;
		// Every node, internal and leaf. The internal conditionals involve this tree alone; the
		// leaf conditionals additionally carry the other tree's parent term (Phase 2). Dropping
		// the leaves is still available via --no_leaf_pseudo_likelihood, because leaving theta's
		// only information in the latent Z layer is what closes the feedback loop that biases nu.
		const auto &nodes_to_score =
		    ProgramOptions::LEAF_PSEUDO_LIKELIHOOD ? _all_nodes_indices : _internal_nodes;
		for (const size_t i : nodes_to_score) {
			const auto &node   = _nodes[i];
			bool state_of_node = current_state.get(i);

			// Note: need to take oldValue because we update _binned_branch_length before
			// starting the loop!!!
			if (!node.is_root()) {
				branch_len_bin = _binned_branch_lengths->oldValue(
				    _leaves_and_internal_nodes_without_roots_indices[i]);
			} else {
				branch_len_bin = std::nullopt;
			}

			if constexpr (IsAlpha) {
				_add_pseudo_likelihood_term<false>(node, i, clique, state_of_node, LL_old,
				                                   current_state, branch_len_bin, old_value, c);
				_add_pseudo_likelihood_term<true>(node, i, clique, state_of_node, LL_new,
				                                  current_state, branch_len_bin, new_value, c);
			} else {
				_add_pseudo_likelihood_term<false>(node, i, clique, state_of_node, LL_old,
				                                   current_state, branch_len_bin,
				                                   _alpha_c->value(c), c);
				_add_pseudo_likelihood_term<true>(node, i, clique, state_of_node, LL_new,
				                                  current_state, branch_len_bin,
				                                  _alpha_c->value(c), c);
			}
		}

		// calculate Hastings ratio
		const double LLRatio       = LL_new.getSum() - LL_old.getSum();
		const double logPriorRatio = param->getLogDensityRatio(c);
		const double logH          = LLRatio + logPriorRatio;

		// accept or reject
		bool accepted = param->acceptOrReject(logH, coretools::TRange(c));
		if (accepted) {
			_cliques[c].accept_update_mu();
			if constexpr (!IsAlpha) { _nu_c[c] = new_value; }
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

public:
	TTree(size_t dimension, const std::string &filename, const std::string &tree_name,
	      TypeParamAlpha *Alpha, TypeParamLogNu *LogNu,
	      TypeParamBinBranches *Binned_Branch_Lenghts);
	~TTree() override;

	[[nodiscard]] size_t size() const { return _nodes.size(); };

	/** Get node by its id
	 * @param Id: the id of the node
	 * @return a reference to the node with the given id
	 */
	[[nodiscard]] const TNode &get_node(const std::string &Id) const;

	/** Get node by its index
	 * @param Id: the index of the node wihtin the tree
	 * @return a reference to the Node with the given id
	 */
	[[nodiscard]] const TNode &get_node(size_t index) const;
	[[nodiscard]] bool isLeaf(size_t index) const;

	/** Get the index of a node by its id
	 * @param Id: the id of the node
	 * @return the index of the node with the given id
	 */
	[[nodiscard]] size_t get_node_index(const std::string &Id) const;

	/** Gives the number of roots within the tree
	 * @return the number of roots
	 */
	[[nodiscard]] size_t number_of_roots() const { return _roots.size(); }

	/** Method to get all the leaves of the tree.
	 * @return Returns a vector of length equal to the number of leaves in the tree. Each
	 * element of the vector is the index of the leaf node within the tree.
	 */
	[[nodiscard]] const std::vector<size_t> &get_leaf_nodes() const { return _leaves; }

	/** @return the number of leaves in the tree
	 */
	[[nodiscard]] size_t get_number_of_leaves() const { return _leaves.size(); }
	[[nodiscard]] size_t get_number_of_nodes() const { return _nodes.size(); }
	[[nodiscard]] size_t get_number_of_internal_nodes() const { return _internal_nodes.size(); }
	[[nodiscard]] size_t get_number_of_roots() const { return _roots.size(); }

	/** @param node_index: the index of the node within the tree
	 * @return The index of the node within the leaves vector (which is smaller than the total
	 * number of nodes in the tree). If the node is not a leaf, the function will return -1.
	 */
	[[nodiscard]] size_t get_index_within_leaves(size_t node_index) const {
		return _leafIndices[node_index];
	}
	[[nodiscard]] size_t get_index_within_leaves(const std::string &node_name) const {
		return _leafIndices[get_node_index(node_name)];
	}
	[[nodiscard]] size_t get_node_index_from_leaf_index(size_t leaf_index) const {
		return _leaves[leaf_index];
	}
	[[nodiscard]] size_t get_node_index_from_internal_nodes_index(size_t internal_index) const {
		return _internal_nodes[internal_index];
	}

	/** @param node_index: the index of the node within the tree
	 * @return The index of the node within the internal nodes vector (which is smaller than the
	 * total number of nodes in the tree). If the node is not an internal node, the function
	 * will return -1.
	 */
	[[nodiscard]] size_t get_index_within_internal_nodes(size_t node_index) const {
		return _internalIndices[node_index];
	}

	/** @return The root nodes of the tree
	 */
	[[nodiscard]] const std::vector<size_t> &get_root_nodes() const { return _roots; }
	[[nodiscard]] const std::vector<size_t> &get_internal_nodes() const { return _internal_nodes; }
	[[nodiscard]] const std::vector<size_t> &get_internal_indicies() const {
		return _internalIndices;
	}
	[[nodiscard]] const std::vector<size_t> &get_internal_nodes_without_roots() const {
		return _internal_nodes_without_roots;
	}

	/// @brief Cross-check the cross-tree lookup used by the leaf pseudo-likelihood terms.
	///
	/// Recomputes P_other(parent -> s) for every Y cell by a route that does not assume the
	/// clique/leaf role swap, and compares. Returns the number of disagreeing cells; 0 is the
	/// only acceptable answer. Run with --check_leaf_lookup.
	size_t check_leaf_lookup_consistency();

	/// @return Internal nodes ordered children-before-parents (see _internal_nodes_post_order).
	[[nodiscard]] const std::vector<size_t> &get_internal_nodes_post_order() const {
		return _internal_nodes_post_order;
	}

	/** Checks whether a node is in the tree
	 * @param node_id: the id of the node
	 * @return true if the node is in the tree, false otherwise
	 */
	[[nodiscard]] bool in_tree(const std::string &node_id) const {
		return _node_map.find(node_id) != _node_map.end();
	};

	[[nodiscard]] std::vector<size_t> get_all_binned_branch_lengths_from_tree() const {
		return _binned_branch_lengths_from_tree;
	}

	// stattools stuff
	[[nodiscard]] std::string name() const override;
	void initialize() override;
	[[nodiscard]] double getSumLogPriorDensity(const Storage &) const override;
	void guessInitialValues() override;
	[[nodiscard]] double getDensity(const Storage &, size_t) const override;
	[[nodiscard]] double getLogDensityRatio(const UpdatedStorage &, size_t) const override;

	void initialize_cliques_and_Z(const std::vector<std::unique_ptr<TTree>> &all_trees);

	[[nodiscard]] double get_delta() const { return _delta; }
	[[nodiscard]] size_t get_number_of_bins() const { return _number_of_bins; }
	std::vector<TClique> &get_cliques();
	[[nodiscard]] const TClique &get_clique_by_index(size_t clique_index) const {
		return _cliques[clique_index];
	}
	[[nodiscard]] const TClique &get_clique(const IndexArray &index_in_leaves_space) const;
	TClique &get_clique(const IndexArray &index_in_leaves_space);
	[[nodiscard]] const TStorageZMatrix &get_Z() const;
	TStorageZMatrix &get_Z();

	[[nodiscard]] std::string get_node_id(size_t index) const { return _nodes[index].get_id(); }

	template<bool IsSimulation, bool FixZ>
	void update_Z_and_nus_and_alphas_and_branch_lengths(const TStorageYMatrix &Y) {
		_reset_joint_log_prob_density();
		// The leaf terms of the pseudo-likelihood read the other tree's parent state for every
		// Y cell. Transpose those bits once here rather than chasing them per node inside the
		// parallel region (see _other_parent_state).
		if constexpr (!IsSimulation) {
			if (ProgramOptions::LEAF_PSEUDO_LIKELIHOOD) { _refresh_other_parent_state(); }
		}
		std::vector<std::vector<size_t>> indices_to_insert(this->_cliques.size());

		// build pairs of branch lengths to update
		auto pairs         = _build_pairs_branch_lengths();
		const auto n_pairs = pairs.length();
		std::vector<std::vector<coretools::TSumLogProbability>> log_sum_per_thread(
		    ProgramOptions::NUMBER_OF_THREADS, std::vector<coretools::TSumLogProbability>(n_pairs));

		// propose new branch lengths
		if constexpr (!IsSimulation) { _propose_new_branch_lengths(pairs); }

#pragma omp parallel for num_threads(ProgramOptions::NUMBER_OF_THREADS) default(none)              \
    schedule(dynamic) shared(pairs, log_sum_per_thread, Y, indices_to_insert)
		for (size_t i = 0; i < _cliques.size(); ++i) {
			auto &log_sum_local = log_sum_per_thread[omp_get_thread_num()];
			// fill the current state for this clique
			auto current_state  = _cliques[i].create_current_state(Y, _Z, *this);
			// update Z
			if constexpr (!FixZ) {
				indices_to_insert[i] = _cliques[i].update_Z(
				    _joint_log_prob_density, current_state, _Z, this, _alpha_c->value(i),
				    _binned_branch_lengths, _leaves_and_internal_nodes_without_roots_indices);
			}

			// update nu and alpha
			if constexpr (!IsSimulation) {
				_update_nu_or_alpha<true>(current_state, i, _alpha_c);
				_update_nu_or_alpha<false>(current_state, i, _log_nu_c);

				// add to likelihood ratio for branch length
				_add_to_LL_branch_lengths(i, current_state, log_sum_local, pairs);
			}
		}

		// update branch lengths
		if constexpr (!IsSimulation) {
			auto log_sum_b = TTree::_reduce_log_sum_per_thread(log_sum_per_thread, n_pairs);
			_evalute_update_branch_length(log_sum_b, pairs);
		}

		if constexpr (!FixZ) { _Z.insert_in_Z(indices_to_insert); }
	}

	[[nodiscard]] TypeBinnedBranchLengths get_binned_branch_length(size_t index_in_tree) const {
		return _binned_branch_lengths->value(
		    _leaves_and_internal_nodes_without_roots_indices[index_in_tree]);
	}

	[[nodiscard]] const std::string &get_tree_name() const { return _tree_name; }

	void simulate_Z(size_t tree_index);

	template<bool WriteFullZ>
	void write_Z_to_file(const std::string &filename, std::vector<std::unique_ptr<TTree>> &trees,
	                     size_t dimension_number_of_tree) const {
		std::vector<std::string> header;
		header.reserve(trees.size());
		for (const auto &tree : trees) { header.push_back(tree->get_tree_name()); }
		header.emplace_back("position");
		header.emplace_back("Z_state");

		if constexpr (WriteFullZ) {
			std::array<size_t, 2> line{};
			coretools::TOutputFile file(filename, header, "\t");
			for (size_t i = 0; i < _Z.total_size_of_container_space(); ++i) {
				// a missing cell reads as state 0, so a direct point lookup covers both cases
				line                = {i, _Z.is_one(i)};
				auto multidim_index = _Z.get_multi_dimensional_index(i);
				std::vector<std::string> node_names;
				for (size_t idx = 0; idx < multidim_index.size(); ++idx) {
					if (idx == dimension_number_of_tree) {
						size_t node_idx = trees[idx]->get_node_index_from_internal_nodes_index(
						    multidim_index[idx]);
						node_names.push_back(trees[idx]->get_node_id(node_idx));
					} else {
						size_t node_idx =
						    trees[idx]->get_node_index_from_leaf_index(multidim_index[idx]);
						node_names.push_back(trees[idx]->get_node_id(node_idx));
					};
				};
				file.writeln(node_names, line);
			}
		} else {
			std::array<size_t, 2> line{};
			coretools::TOutputFile file(filename, header, "\t");
			// iterate only the stored (non-default) cells, in ascending linear-index order
			for (const auto &[linear_index_in_Z_space, storage] : _Z.get_stored_entries()) {
				const auto state    = storage.is_one();
				line                = {linear_index_in_Z_space, state};
				auto multidim_index = _Z.get_multi_dimensional_index(linear_index_in_Z_space);
				std::vector<std::string> node_names;
				for (size_t idx = 0; idx < multidim_index.size(); ++idx) {
					if (idx == dimension_number_of_tree) {
						size_t node_idx = trees[idx]->get_node_index_from_internal_nodes_index(
						    multidim_index[idx]);
						node_names.push_back(trees[idx]->get_node_id(node_idx));
					} else {
						size_t node_idx =
						    trees[idx]->get_node_index_from_leaf_index(multidim_index[idx]);
						node_names.push_back(trees[idx]->get_node_id(node_idx));
					};
				};
				file.writeln(node_names, line);
			}
		}

		if (ProgramOptions::WRITE_BRANCH_LENGTHS) {
			std::vector<std::string> header_branch_len = {"grid_position", "branch_length"};
			coretools::TOutputFile branch_len_file("acol_simulated_" + get_tree_name() +
			                                           "_branch_length_grid.txt",
			                                       header_branch_len, "\t");
			for (size_t i = 0; i < _grid_branch_lengths.size(); ++i) {
				branch_len_file.writeln(i, _grid_branch_lengths[i]);
			}
		}
	}

	[[nodiscard]] double get_complete_joint_density() const {
		return coretools::containerSum(_joint_log_prob_density);
	}

	void initialize_Z_from_children(const TStorageYMatrix &Y) {
		std::string set_Z_cli_command = "set_" + get_tree_name() + "_Z";
		if (coretools::instances::parameters().exists(set_Z_cli_command)) { return; }

		// Each clique is independent of each other so we should be able to parallelize this
		std::vector<std::vector<size_t>> indices_to_insert(this->_cliques.size());

#pragma omp parallel for num_threads(ProgramOptions::NUMBER_OF_THREADS)                            \
    schedule(dynamic) default(none) shared(indices_to_insert, Y)
		for (size_t i = 0; i < _cliques.size(); ++i) {
			auto current_state   = _cliques[i].create_current_state(Y, _Z, *this);
			indices_to_insert[i] = _cliques[i].initialize_Z_from_children(
			    current_state, _Z, this, _binned_branch_lengths,
			    _leaves_and_internal_nodes_without_roots_indices);
		}

		_Z.insert_in_Z(indices_to_insert);
	};

	[[nodiscard]] std::vector<double> get_paper_counts() const {
		std::string parameter_name = get_tree_name() + "_paper_counts";
		if (!coretools::instances::parameters().exists(parameter_name)) {
			throw coretools::TUserError("Parameter '", parameter_name,
			                            "' not found. Please provide it.");
		}

		const auto filename = coretools::instances::parameters().get<std::string>(parameter_name);
		coretools::TInputFile file(filename, coretools::FileType::Header);

		if (file.numCols() != 2) {
			throw coretools::TUserError("File '", filename,
			                            "' is expected to have 2 columns, but has ", file.numCols(),
			                            " !");
		}

		// now we initilise the vector of paper counts. The entries should only be leaves
		std::vector<double> paper_counts(get_number_of_leaves(), 0.0);

		for (; !file.empty(); file.popFront()) {
			const std::string leaf_name = std::string(file.get(0));
			const auto count            = file.get<size_t>(1);

			// get the node index from the leaf name
			const size_t node_index = this->get_node_index(leaf_name);
			if (!isLeaf(node_index)) { throw coretools::TUserError("All nodes should be leaves."); }
			const size_t leaf_index = this->get_index_within_leaves(node_index);

			if (leaf_index >= paper_counts.size()) {
				throw coretools::TUserError("Leaf index ", leaf_index,
				                            " is out of bounds for paper counts vector of size ",
				                            paper_counts.size(), ".");
			}

			paper_counts[leaf_index] = std::log(count + 1);
		}

		return paper_counts;
	}
};
#endif // METABOLITE_INFERENCE_TREE_H
