//
// Created by Marco Visani on 26.06.23.
//

#include "TTree.h"
#include "TClique.h"
#include "cli.h"
#include "constants.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Math/TSumLog.h"
#include "coretools/Types/probability.h"

#include <cstddef>
#include <cstdlib>
#include <queue>
#include <string>
#include <vector>

TTree::TTree(size_t dimension, const std::string &filename, const std::string &tree_name,
             TypeParamAlpha *Alpha, TypeParamLogNu *LogNu,
             TypeParamBinBranches *Binned_Branch_Lenghts)
    : _dimension(dimension), _binned_branch_lengths(Binned_Branch_Lenghts), _log_nu_c(LogNu),
      _alpha_c(Alpha) {

	// tell stattools that these parameters belong to a prior distribution
	this->addPriorParameter({_binned_branch_lengths, _alpha_c, _log_nu_c});

	_load_from_file(filename, tree_name);
}

TTree::~TTree() = default;

const TNode &TTree::get_node(const std::string &Id) const {
	auto it = _node_map.find(Id);
	if (it == _node_map.end()) { throw coretools::TUserError("Node '", Id, "' does not exist!"); }
	return _nodes[it->second]; // Retrieve node from vector using the index
}

const TNode &TTree::get_node(size_t index) const { return _nodes[index]; }
bool TTree::isLeaf(size_t index) const { return _nodes.at(index).is_leaf(); }

size_t TTree::get_node_index(const std::string &Id) const {
	auto it = _node_map.find(Id);
	if (it == _node_map.end()) {
		throw coretools::TUserError("Node '", Id, "' does not exist in the tree !");
	}
	return it->second; // Return the index from the map
}

void TTree::initialize_cliques_and_Z(const std::vector<std::unique_ptr<TTree>> &all_trees) {
	// kept for the leaf terms of the pseudo-likelihood, which need the other tree
	_all_trees = &all_trees;

	// we initialize the number of leaves we have in each tree
	IndexArray num_leaves_per_tree;
	for (size_t i = 0; i < all_trees.size(); ++i) {
		num_leaves_per_tree[i] = all_trees[i]->get_number_of_leaves();
	}

	_initialize_Z(num_leaves_per_tree);
	_initialize_cliques(num_leaves_per_tree, all_trees);
}

[[nodiscard]] std::string TTree::name() const { return "tree"; }

void TTree::initialize() {
	// stattools initialization function
	_alpha_c->initStorage(this, {_cliques.size()},
	                      {std::make_shared<coretools::TNamesStrings>(_clique_names)});

	// now we initialize the mu_c_1
	_log_nu_c->initStorage(this, {_cliques.size()},
	                       {std::make_shared<coretools::TNamesStrings>(_clique_names)});
	_nu_c.resize(_cliques.size());
	for (size_t c = 0; c < _cliques.size(); ++c) { _nu_c[c] = std::exp(_log_nu_c->value(c)); }

	// number of branches = number of leaves + number of internal nodes without roots
	std::vector<std::string> branch_names;
	branch_names.reserve(_leaves_and_internal_nodes_without_roots.size());
	for (size_t node_idx : _leaves_and_internal_nodes_without_roots) {
		branch_names.push_back(get_node_id(node_idx));
	}
	_binned_branch_lengths->initStorage(this, {get_number_of_nodes() - get_number_of_roots()},
	                                    {std::make_shared<coretools::TNamesStrings>(branch_names)});
}

void TTree::_refresh_other_parent_state() {
	// Only meaningful for two trees: with more, a Y cell would have one parent term per tree and
	// the leaf conditional would need all of them.
	static_assert(NUMBER_OF_TREES == 2,
	              "the leaf pseudo-likelihood term assumes exactly two trees");
	const TTree &other = *(*_all_trees)[1 - _dimension];

	const size_t n_cliques = _cliques.size();
	const size_t n_leaves  = get_number_of_leaves();
	_other_parent_state.assign(n_cliques * n_leaves, false);

	const auto &other_Z = other.get_Z();
	for (size_t clique = 0; clique < n_cliques; ++clique) {
		// our clique index is a leaf of the other tree; we need that leaf's parent there
		const size_t node_in_other  = other.get_node_index_from_leaf_index(clique);
		const size_t parent_ix      = other.get_node(node_in_other).parent_index_in_tree();
		const size_t parent_interna = other.get_index_within_internal_nodes(parent_ix);

		IndexArray index_in_other_Z{};
		index_in_other_Z[other._dimension] = parent_interna;
		for (size_t leaf = 0; leaf < n_leaves; ++leaf) {
			// symmetrically, our leaf indexes the other tree's clique dimension
			index_in_other_Z[_dimension] = leaf;
			_other_parent_state[clique * n_leaves + leaf] =
			    other_Z.is_one(other_Z.get_linear_index_in_container_space(index_in_other_Z));
		}
	}
}

std::array<double, 2> TTree::_prob_other_tree_at_leaf(size_t clique,
                                                      size_t index_within_leaves) const {
	const TTree &other = *(*_all_trees)[1 - _dimension];

	// The roles swap across the two trees: our clique index is a leaf over there, and our leaf
	// index selects which of its cliques (hence which alpha and nu) applies to this cell.
	const size_t node_in_other = other.get_node_index_from_leaf_index(clique);
	const auto &other_clique   = other.get_clique_by_index(index_within_leaves);
	// value(), not oldValue(): the other tree's branch lengths are not mid-proposal here, only
	// our own are (they are proposed at the top of our own update).
	const auto bin          = other.get_binned_branch_length(node_in_other);
	const bool parent_state = _other_parent_state[clique * get_number_of_leaves() +
	                                              index_within_leaves];
	const auto &matrix      = other_clique.get_matrix<false>(bin);
	return {matrix(parent_state, 0), matrix(parent_state, 1)};
}

size_t TTree::check_leaf_lookup_consistency() {
	// _prob_other_tree_at_leaf reaches the other tree by swapping roles: our clique index is a
	// leaf over there, our leaf index selects which of its cliques applies. That swap is the one
	// piece of Phase 2 that the neutral-tree control cannot exercise -- when the other tree is
	// neutral its factor is {0.5, 0.5} whatever cell is read, so a wrong index is invisible.
	//
	// So recompute the same quantity by a route that does not assume the swap: build the cell's
	// index in Y space from (clique, leaf), and from that derive the other tree's clique with
	// get_clique(IndexArray) -- the same call the Y sweep makes in _calculate_log_prob_field --
	// its leaf node, and its parent's state. The two must agree cell for cell.
	_refresh_other_parent_state();

	const TTree &other       = *(*_all_trees)[1 - _dimension];
	const size_t n_cliques   = _cliques.size();
	const size_t n_leaves    = get_number_of_leaves();
	const auto &other_Z      = other.get_Z();
	size_t n_mismatch        = 0;
	size_t first_bad_clique  = 0;
	size_t first_bad_leaf    = 0;

	for (size_t clique = 0; clique < n_cliques; ++clique) {
		for (size_t leaf = 0; leaf < n_leaves; ++leaf) {
			// the cell, in Y space: our leaf runs along our own dimension, our clique index
			// along the other one
			IndexArray index_in_Y{};
			index_in_Y[_dimension]     = leaf;
			index_in_Y[1 - _dimension] = clique;

			// --- independent route ---
			const size_t other_leaf_node =
			    other.get_node_index_from_leaf_index(index_in_Y[1 - _dimension]);
			const size_t other_parent = other.get_node(other_leaf_node).parent_index_in_tree();
			IndexArray index_in_other_Z    = index_in_Y;
			index_in_other_Z[1 - _dimension] = other.get_index_within_internal_nodes(other_parent);
			const bool parent_state =
			    other_Z.is_one(other_Z.get_linear_index_in_container_space(index_in_other_Z));
			const auto &matrix = other.get_clique(index_in_Y).get_matrix<false>(
			    other.get_binned_branch_length(other_leaf_node));
			const std::array<double, 2> expected{matrix(parent_state, 0), matrix(parent_state, 1)};

			// --- the route Phase 2 actually uses ---
			const auto got = _prob_other_tree_at_leaf(clique, leaf);

			if (std::abs(got[0] - expected[0]) > 1e-12 ||
			    std::abs(got[1] - expected[1]) > 1e-12) {
				if (n_mismatch == 0) {
					first_bad_clique = clique;
					first_bad_leaf   = leaf;
				}
				++n_mismatch;
			}
		}
	}

	if (n_mismatch == 0) {
		coretools::instances::logfile().list("Leaf lookup check (", get_tree_name(), "): all ",
		                                     n_cliques * n_leaves, " cells agree.");
	} else {
		coretools::instances::logfile().list(
		    "Leaf lookup check (", get_tree_name(), "): ", n_mismatch, " of ",
		    n_cliques * n_leaves, " cells DISAGREE, first at clique ", first_bad_clique,
		    " leaf ", first_bad_leaf, ".");
	}
	return n_mismatch;
}

void TTree::guessInitialValues() {
	for (size_t c = 0; c < _cliques.size(); ++c) {
		// Draw log_nu[c] ~ Normal(LOG_NU_C, LOG_NU_C_INIT_SD^2) instead of setting every clique to
		// the same constant. Identical initial values would make the MLE that seeds var_log_nu 0,
		// yielding a degenerate prior that freezes log_nu, mean_log_nu and var_log_nu (their
		// proposals would never be accepted).
		const double log_nu_init = coretools::instances::randomGenerator().getNormalRandom(
		    ProgramOptions::LOG_NU_C, ProgramOptions::LOG_NU_C_INIT_SD);
		_log_nu_c->set(c, log_nu_init);
		_alpha_c->set(c, coretools::Probability(ProgramOptions::ALPHA));
		_nu_c[c] = std::exp(_log_nu_c->value(c));
		_cliques[c].set_lambda(_alpha_c->value(c), _nu_c[c]);
	}

	_set_initial_branch_lengths(false);
}

double TTree::getSumLogPriorDensity(const Storage &) const {
	throw coretools::TDevError("Should never be called");
}
double TTree::getDensity(const Storage &, size_t) const {
	throw coretools::TDevError("Should never be called");
}
double TTree::getLogDensityRatio(const UpdatedStorage &, size_t) const {
	throw coretools::TDevError("Should never be called");
}

void TTree::_simulateUnderPrior(Storage *) {
	using namespace coretools::instances;
	_set_initial_branch_lengths(true);
	for (size_t c = 0; c < _cliques.size(); ++c) {
		_nu_c[c] = std::exp(_log_nu_c->value(c));
		_cliques[c].set_lambda(_alpha_c->value(c), _nu_c[c]);
	}
}

void TTree::_initialize_Z(IndexArray num_leaves_per_tree) {
	num_leaves_per_tree[_dimension] = this->get_number_of_internal_nodes();

	_Z.initialize_dimensions(num_leaves_per_tree);

	std::string set_Z_cli_command = "set_" + get_tree_name() + "_Z";
	if (coretools::instances::parameters().exists(set_Z_cli_command)) {
		std::string filename = coretools::instances::parameters().get(set_Z_cli_command);
		coretools::TInputFile file(filename, coretools::FileType::Header);

		if (file.numCols() != 4) {
			throw coretools::TUserError("The file for setting Z must have 4 columns ! ");
		}
		// read each line of the file
		for (; !file.empty(); file.popFront()) {
			auto linear_index_in_Z_space = file.get<uint32_t>(2);
			bool state                   = file.get<bool>(3);
			if (state) { _Z.insert_one(linear_index_in_Z_space); }
		}
		return;
	}
}

const TStorageZMatrix &TTree::get_Z() const { return _Z; };
TStorageZMatrix &TTree::get_Z() { return _Z; };

void TTree::simulate_Z(size_t tree_index) {
	for (size_t c = 0; c < _cliques.size(); ++c) {
		auto &clique = _cliques[c];
		_simulation_prepare_cliques(c, clique);
		TCurrentState current_state(*this, clique.get_increment(), get_number_of_leaves(),
		                            get_number_of_internal_nodes());

		// we sample the roots
		if (ProgramOptions::SIMULATION_NO_Z_INITIALIZATION) { continue; }
		double proba_root = TClique::get_stationary_probability(true, _alpha_c->value(c));
		coretools::Probability p(proba_root);

		// we can also prepare the queue for the DFS
		std::queue<size_t> node_queue;
		for (const auto root_index_in_tree : this->get_root_nodes()) {
			bool root_state = coretools::instances::randomGenerator().pickOneOfTwo(p);
			if (root_state) {
				_simulate_one(clique, current_state, tree_index, root_index_in_tree);
			}
			for (const auto child : this->get_node(root_index_in_tree).children_indices_in_tree()) {
				if (!this->isLeaf(child)) { node_queue.push(child); }
			} // those are the first children of the tree (children of the roots).
		} // roots done, we go to the internal nodes

		// sampling the internal nodes
		while (!node_queue.empty()) {
			size_t node_index = node_queue.front();
			node_queue.pop();
			const TNode &node = this->get_node(node_index);

			// we want to sample the state of the node given its parent (and independently of its
			// children since we haven't sampled them yet).
			std::array<coretools::TSumLogProbability, 2> sum_log;
			clique.calculate_log_prob_parent_to_node(
			    node_index,
			    (TypeBinnedBranchLengths)_binned_branch_lengths->value(
			        _leaves_and_internal_nodes_without_roots_indices[node_index]),
			    this, 0, current_state, sum_log);
			bool internal_node_state = sample(sum_log);
			if (internal_node_state) {
				_simulate_one(clique, current_state, tree_index, node_index);
			}

			for (size_t child_index : node.children_indices_in_tree()) {
				if (!this->isLeaf(child_index)) {
					node_queue.push(child_index);
				} // as long as your are not a leaf we can continue sampling Z
			}
		} // internal nodes done, we go to the leaves
	}
}

void TTree::_simulate_one(const TClique &clique, TCurrentState &current_state, size_t tree_index,
                          size_t node_index_in_tree) {
	auto index_in_leaves_space        = clique.get_start_index_in_leaf_space();
	index_in_leaves_space[tree_index] = this->get_index_within_internal_nodes(node_index_in_tree);
	_Z.insert_one(index_in_leaves_space);
	current_state.set(node_index_in_tree, true);
}
