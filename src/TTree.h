//
// Created by Marco Visani on 26.06.23.
//

#ifndef METABOLITE_INFERENCE_TREE_H
#define METABOLITE_INFERENCE_TREE_H

#include "TClique.h"
#include "TStorageYVector.h"
#include "TStorageZVector.h"
#include "Types.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Files/TOutputFile.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Math/TSumLog.h"
#include "coretools/Types/commonWeakTypes.h"
#include "coretools/Types/probability.h"
#include "coretools/algorithms.h"
#include "coretools/devtools.h"
#include "stattools/ParametersObservations/TParameter.h"
#include <cstddef>
#include <list>
#include <string>
#include <vector>

class TNode {
private:
	std::string _id;                               // unique identifier for the node
	size_t _parentIndex_in_tree;                   // pointer to parent node
	std::vector<size_t> _children_indices_in_tree; // vector to child nodes indices
	bool _is_root;

public:
	// Constructor, initializes a TNode with id, _branchLengthToParent and parent node
	TNode(std::string IdString, size_t Parent_index_in_tree, bool is_root)
	    : _id(std::move(IdString)), _parentIndex_in_tree(Parent_index_in_tree), _is_root(is_root) {};

	// Function to add child
	void addChild_index_in_tree(size_t ChildIndex) { _children_indices_in_tree.push_back(ChildIndex); };

	// Method to check if the node is a leaf (has no _children)
	[[nodiscard]] const std::string &id() const { return _id; };
	[[nodiscard]] size_t parentIndex_in_tree() const { return _parentIndex_in_tree; };
	[[nodiscard]] size_t numChildren() const { return _children_indices_in_tree.size(); };
	[[nodiscard]] const std::vector<size_t> &children_indices_in_tree() const { return _children_indices_in_tree; };
	[[nodiscard]] inline bool isLeaf() const { return _children_indices_in_tree.empty(); };
	[[nodiscard]] bool isRoot() const { return _is_root; };
	[[nodiscard]] bool isInternalNode() const { return !_is_root && !_children_indices_in_tree.empty(); };
	void set_is_root(bool is_root) { _is_root = is_root; }
	void set_parent_index_in_tree(size_t parent_index_in_tree) { _parentIndex_in_tree = parent_index_in_tree; }
	std::string get_id() const { return _id; };

	bool operator==(const std::string &Id) const { return _id == Id; };
};

/// Note: All indices are within the tree itself
class TTree : public stattools::prior::TStochasticBase<TypeMarkovField, NumDimMarkovField> {
public:
	// some type aliases, for better readability
	using BoxType = TTree;
	using Base    = stattools::prior::TStochasticBase<TypeMarkovField, NumDimMarkovField>;
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
	std::vector<size_t> _dimension_cliques;

	// Nus
	TypeParamLogNu *_log_nu_c = nullptr;
	std::vector<double> _nu_c;

	// Alphas
	TypeParamAlpha *_alpha_c = nullptr;

	// Set Z
	TStorageZVector _Z;

	// Joint probability density
	std::vector<double> _joint_log_prob_density;

	// private functions
	void _reset_joint_log_prob_density() {
		_joint_log_prob_density.clear();
		_joint_log_prob_density.resize(NUMBER_OF_THREADS);
	}
	void _set_initial_branch_lengths(bool is_simulation);
	std::vector<size_t> _bin_branch_lengths(const std::vector<double> &branch_lengths, bool exclude_root) const;
	void _bin_branch_lengths_from_tree(std::vector<double> &branch_lengths);
	void _initialize_grid_branch_lengths(size_t number_of_branches);
	void _initialize_Z(std::vector<size_t> num_leaves_per_tree);
	void _initialize_cliques(const std::vector<size_t> &num_leaves_per_tree,
	                         const std::vector<std::unique_ptr<TTree>> &all_trees);
	void _load_from_file(const std::string &filename, const std::string &tree_name);
	void _simulation_prepare_cliques(size_t c, TClique &clique) const;
	void _simulate_one(const TClique &clique, TCurrentState &current_state, size_t tree_index,
	                   size_t node_index_in_tree);

	// updating branch lengths
	stattools::TPairIndexSampler _build_pairs_branch_lengths() const;
	void _propose_new_branch_lengths(const stattools::TPairIndexSampler &pairs);
	void _propose_new_branch_lengths(size_t p1, size_t p2, int val);
	void _add_to_LL_branch_lengths(size_t c, const TCurrentState &current_state,
	                               std::vector<coretools::TSumLogProbability> &log_sum,
	                               const stattools::TPairIndexSampler &pairs) const;
	double _calculate_likelihood_ratio_branch_length(size_t index_in_binned_branch_length, const TClique &clique,
	                                                 const TCurrentState &current_state) const;

	void _simulateUnderPrior(Storage *) override;

	template<bool UseTryMatrix>
	void _compute_LL_old_and_new_mu(size_t index_in_tree, const TClique &clique, bool state_of_node,
	                                coretools::TSumLogProbability &LL, const TCurrentState &current_state,
	                                size_t branch_len_bin, double alpha) {
		if (_nodes[index_in_tree].isRoot()) {
			LL.add(TClique::get_stationary_probability(state_of_node, alpha));
		} else {
			double prob =
			    clique.calculate_prob_to_parent<UseTryMatrix>(index_in_tree, this, branch_len_bin, current_state);
			LL.add(prob);
		}
	}

	template<bool IsAlpha, typename TypeParam>
	void _update_mu(const TCurrentState &current_state, size_t c, TypeParam *param) {
		// propose a new value
		param->propose(coretools::TRange(c));

		// calculate LL for old mu
		// No need to change Lambda (rate matrix), just go over entire tree and calculate probabilities
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
		for (size_t i = 0; i < _nodes.size(); ++i) {
			bool state_of_node = current_state.get(i);
			// Note: need to take oldValue because we update _binned_branch_length before starting the loop!!!
			const auto branch_len_bin =
			    _binned_branch_lengths->oldValue(_leaves_and_internal_nodes_without_roots_indices[i]);

			if constexpr (IsAlpha) {
				_compute_LL_old_and_new_mu<false>(i, clique, state_of_node, LL_old, current_state, branch_len_bin,
				                                  old_value);
				_compute_LL_old_and_new_mu<true>(i, clique, state_of_node, LL_new, current_state, branch_len_bin,
				                                 new_value);
			} else {
				_compute_LL_old_and_new_mu<false>(i, clique, state_of_node, LL_old, current_state, branch_len_bin,
				                                  _alpha_c->value(c));
				_compute_LL_old_and_new_mu<true>(i, clique, state_of_node, LL_new, current_state, branch_len_bin,
				                                 _alpha_c->value(c));
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

public:
	TTree(size_t dimension, const std::string &filename, const std::string &tree_name, TypeParamAlpha *Alpha,
	      TypeParamLogNu *LogNu, TypeParamBinBranches *Binned_Branch_Lenghts);
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
	const TNode &get_node(size_t index) const;
	bool isLeaf(size_t index) const;

	/** Get the index of a node by its id
	 * @param Id: the id of the node
	 * @return the index of the node with the given id
	 */
	size_t get_node_index(const std::string &Id) const;

	/** Gives the number of roots within the tree
	 * @return the number of roots
	 */
	[[nodiscard]] size_t number_of_roots() const { return _roots.size(); }

	/** Method to get all the leaves of the tree.
	 * @return Returns a vector of length equal to the number of leaves in the tree. Each element of the vector is the
	 * index of the leaf node within the tree.
	 */
	[[nodiscard]] const std::vector<size_t> &get_leaf_nodes() const { return _leaves; }

	/** @return the number of leaves in the tree
	 */
	size_t get_number_of_leaves() const { return _leaves.size(); }
	size_t get_number_of_nodes() const { return _nodes.size(); }
	size_t get_number_of_internal_nodes() const { return _internal_nodes.size(); }
	size_t get_number_of_roots() const { return _roots.size(); }

	/** @param node_index: the index of the node within the tree
	 * @return The index of the node within the leaves vector (which is smaller than the total number of nodes in the
	 * tree). If the node is not a leaf, the function will return -1.
	 */
	size_t get_index_within_leaves(size_t node_index) const { return _leafIndices[node_index]; }
	size_t get_index_within_leaves(const std::string &node_name) const {
		return _leafIndices[get_node_index(node_name)];
	}
	size_t get_node_index_from_leaf_index(size_t leaf_index) const { return _leaves[leaf_index]; }
	size_t get_node_index_from_internal_nodes_index(size_t internal_index) const {
		return _internal_nodes[internal_index];
	}

	/** @param node_index: the index of the node within the tree
	 * @return The index of the node within the internal nodes vector (which is smaller than the total number of nodes
	 * in the tree). If the node is not an internal node, the function will return -1.
	 */
	size_t get_index_within_internal_nodes(size_t node_index) const { return _internalIndices[node_index]; }

	/** @return The root nodes of the tree
	 */
	[[nodiscard]] const std::vector<size_t> &get_root_nodes() const { return _roots; }
	const std::vector<size_t> &get_internal_nodes() const { return _internal_nodes; }
	const std::vector<size_t> &get_internal_indicies() const { return _internalIndices; }
	const std::vector<size_t> &get_internal_nodes_without_roots() const { return _internal_nodes_without_roots; }

	/** Checks whether a node is in the tree
	 * @param node_id: the id of the node
	 * @return true if the node is in the tree, false otherwise
	 */
	bool in_tree(const std::string &node_id) const { return _node_map.find(node_id) != _node_map.end(); };

	std::vector<size_t> get_all_binned_branch_lengths_from_tree() const { return _binned_branch_lengths_from_tree; }

	// stattools stuff
	[[nodiscard]] std::string name() const override;
	void initialize() override;
	double getSumLogPriorDensity(const Storage &) const override;
	void guessInitialValues() override;
	double getDensity(const Storage &, size_t) const override;
	double getLogDensityRatio(const UpdatedStorage &, size_t) const override;

	void initialize_cliques_and_Z(const std::vector<std::unique_ptr<TTree>> &all_trees);

	double get_delta() const { return _delta; }
	size_t get_number_of_bins() const { return _number_of_bins; }
	std::vector<TClique> &get_cliques();
	const TClique &get_clique(const std::vector<size_t> &index_in_leaves_space) const;
	TClique &get_clique(const std::vector<size_t> &index_in_leaves_space);
	const TStorageZVector &get_Z() const;
	TStorageZVector &get_Z();

	std::string get_node_id(size_t index) const { return _nodes[index].get_id(); }

	template<bool IsSimulation, bool FixZ> void update_Z_and_mus_and_branch_lengths(const TStorageYVector &Y) {
		_reset_joint_log_prob_density();
		std::vector<std::vector<TStorageZ>> indices_to_insert(this->_cliques.size());

		// build pairs of branch lengths to update
		auto pairs = _build_pairs_branch_lengths();
		std::vector<coretools::TSumLogProbability> log_sum_b(pairs.length());

		// propose new branch lengths
		if constexpr (!IsSimulation) { _propose_new_branch_lengths(pairs); }

#pragma omp parallel for num_threads(NUMBER_OF_THREADS) schedule(static)
		for (size_t i = 0; i < _cliques.size(); ++i) {
			// fill the current state for this clique
			auto current_state = _cliques[i].create_current_state(Y, _Z, *this);
			// update Z
			if constexpr (!FixZ) {
				indices_to_insert[i] =
				    _cliques[i].update_Z(_joint_log_prob_density, current_state, _Z, this, _alpha_c->value(i),
				                         _binned_branch_lengths, _leaves_and_internal_nodes_without_roots_indices);
			}

			// update mu
			if constexpr (!IsSimulation) {
				_update_mu<true>(current_state, i, _alpha_c);
				_update_mu<false>(current_state, i, _log_nu_c);

				// add to likelihood ratio for branch length
				_add_to_LL_branch_lengths(i, current_state, log_sum_b, pairs);
			}
		}
		if constexpr (!FixZ) { _Z.insert_in_Z(indices_to_insert); }

		// update branch lengths
		if constexpr (!IsSimulation) { _evalute_update_branch_length(log_sum_b, pairs); }
	}

	TypeBinnedBranchLengths get_binned_branch_length(size_t index_in_tree) const {
		return _binned_branch_lengths->value(_leaves_and_internal_nodes_without_roots_indices[index_in_tree]);
	}

	const std::string &get_tree_name() const { return _tree_name; }

	void simulate_Z(size_t tree_index);

	template<bool WriteFullZ>
	void write_Z_to_file(const std::string &filename, std::vector<std::unique_ptr<TTree>> &trees,
	                     size_t dimension_number_of_tree) const {
		std::vector<std::string> header;
		for (const auto &tree : trees) { header.push_back(tree->get_tree_name()); }
		header.emplace_back("position");
		header.emplace_back("Z_state");

		if constexpr (WriteFullZ) {
			std::array<size_t, 2> line{};
			coretools::TOutputFile file(filename, header, "\t");
			for (size_t i = 0; i < _Z.total_size_of_container_space(); ++i) {
				auto [found, position] = _Z.binary_search(i);
				if (found) {
					line = {i, _Z.is_one(position)};
				} else {
					line = {i, 0};
				}
				std::vector<size_t> multidim_index = _Z.get_multi_dimensional_index(i);
				std::vector<std::string> node_names;
				for (size_t idx = 0; idx < multidim_index.size(); ++idx) {
					if (idx == dimension_number_of_tree) {
						size_t node_idx = trees[idx]->get_node_index_from_internal_nodes_index(multidim_index[idx]);
						node_names.push_back(trees[idx]->get_node_id(node_idx));
					} else {
						size_t node_idx = trees[idx]->get_node_index_from_leaf_index(multidim_index[idx]);
						node_names.push_back(trees[idx]->get_node_id(node_idx));
					};
				};
				file.writeln(node_names, line);
			}
		} else {
			std::array<size_t, 2> line{};
			coretools::TOutputFile file(filename, header, "\t");
			for (size_t i = 0; i < _Z.size(); ++i) {
				line                               = {_Z[i].get_linear_index_in_Z_space(), _Z[i].is_one()};
				std::vector<size_t> multidim_index = _Z.get_multi_dimensional_index(i);
				std::vector<std::string> node_names;
				for (size_t idx = 0; idx < multidim_index.size(); ++idx) {
					if (idx == dimension_number_of_tree) {
						size_t node_idx = trees[idx]->get_node_index_from_internal_nodes_index(multidim_index[idx]);
						node_names.push_back(trees[idx]->get_node_id(node_idx));
					} else {
						size_t node_idx = trees[idx]->get_node_index_from_leaf_index(multidim_index[idx]);
						node_names.push_back(trees[idx]->get_node_id(node_idx));
					};
				};
				file.writeln(node_names, line);
			}
		}

		if (WRITE_BRANCH_LENGTHS) {
			std::vector<std::string> header_branch_len = {"grid_position", "branch_length"};
			coretools::TOutputFile branch_len_file("acol_simulated_" + get_tree_name() + "_branch_length_grid.txt",
			                                       header_branch_len, "\t");
			for (size_t i = 0; i < _grid_branch_lengths.size(); ++i) {
				branch_len_file.writeln(i, _grid_branch_lengths[i]);
			}
		}
	}

	double get_complete_joint_density() const { return coretools::containerSum(_joint_log_prob_density); }

	void initialize_Z_from_children(const TStorageYVector &Y) {
		std::string set_Z_cli_command = "set_" + get_tree_name() + "_Z";
		if (coretools::instances::parameters().exists(set_Z_cli_command)) { return; }

		// Each clique is independent of each other so we should be able to parallelize this
		std::vector<std::vector<TStorageZ>> indices_to_insert(this->_cliques.size());

#pragma omp parallel for num_threads(NUMBER_OF_THREADS) schedule(static)
		for (size_t i = 0; i < _cliques.size(); ++i) {
			auto current_state   = _cliques[i].create_current_state(Y, _Z, *this);
			indices_to_insert[i] = _cliques[i].initialize_Z_from_children(
			    current_state, _Z, this, _binned_branch_lengths, _leaves_and_internal_nodes_without_roots_indices);
		}

		_Z.insert_in_Z(indices_to_insert);
	};

	std::vector<size_t> get_paper_counts() const {
		std::string parameter_name = get_tree_name() + "_paper_counts";
		if (!coretools::instances::parameters().exists(parameter_name)) {
			UERROR("Parameter '", parameter_name, "' not found. Please provide it.");
		}

		const auto filename = coretools::instances::parameters().get<std::string>(parameter_name);
		coretools::TInputFile file(filename, coretools::FileType::Header);

		if (file.numCols() != 2) {
			UERROR("File '", filename, "' is expected to have 2 columns, but has ", file.numCols(), " !");
		}

		// now we initilise the vector of paper counts. The entries should only be leaves
		std::vector<size_t> paper_counts(get_number_of_leaves(), 0);

		for (; !file.empty(); file.popFront()) {
			const std::string leaf_name = std::string(file.get(0));
			const auto count            = file.get<size_t>(1);

			// get the node index from the leaf name
			const size_t node_index = this->get_node_index(leaf_name);
			if (!isLeaf(node_index)) { UERROR("All nodes should be leaves."); }
			const size_t leaf_index = this->get_index_within_leaves(node_index);

			if (leaf_index >= paper_counts.size()) {
				UERROR("Leaf index ", leaf_index, " is out of bounds for paper counts vector of size ",
				       paper_counts.size(), ".");
			}

			paper_counts[leaf_index] = count;
		}

		return paper_counts;
	}
};
#endif // METABOLITE_INFERENCE_TREE_H
