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
#include "tree/io/read_Z.h"
#include "tree/node_state_shape.h"

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

/// Read the topology, then bin the branch lengths it carries. The logging stays here rather than
/// moving into read_phylogeny: a phylogeny is a value, and a value that writes to a logfile cannot
/// be built in a test.
void TTree::_load_from_file(const std::string &filename, const std::string &tree_name) {
	coretools::instances::logfile().listFlush("Reading tree from file '", filename, "' ...");
	_tree_name = tree_name;
	_phylogeny.emplace(read_phylogeny(filename));
	_bin_branch_lengths_from_tree();
	coretools::instances::logfile().done();
	coretools::instances::logfile().conclude(
	    "Read ", _topology().n_nodes(), " nodes of which ", _topology().n_roots(),
	    " are roots and ", _topology().n_leaves(), " are leaves and ",
	    _topology().internal_nodes_without_roots().size(), " are internal nodes.");
}

void TTree::initialize_cliques_and_Z(const std::vector<std::unique_ptr<TTree>> &all_trees) {

	// we initialize the number of leaves we have in each tree
	IndexArray num_leaves_per_tree;
	for (size_t i = 0; i < all_trees.size(); ++i) {
		num_leaves_per_tree[i] = all_trees[i]->get_number_of_leaves();
	}

	_initialize_Z(num_leaves_per_tree, all_trees);
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
	branch_names.reserve(_topology().branches().size());
	for (size_t node_idx : _topology().branches()) {
		branch_names.push_back(get_node_id(node_idx));
	}
	_binned_branch_lengths->initStorage(this, {get_number_of_nodes() - get_number_of_roots()},
	                                    {std::make_shared<coretools::TNamesStrings>(branch_names)});
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
		_cliques[c].set_transition_grid(TTransitionGrid(_alpha_c->value(c), _nu_c[c], _grid()));
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
		_cliques[c].set_transition_grid(TTransitionGrid(_alpha_c->value(c), _nu_c[c], _grid()));
	}
}

void TTree::_initialize_Z(IndexArray num_leaves_per_tree,
                          const std::vector<std::unique_ptr<TTree>> &all_trees) {
	// The node state spans every node of this tree, leaves included. The rule lives in
	// node_state_dimensions so that the storage tests can assert over the shapes production
	// actually builds rather than over a restatement of them.
	_Z.initialize_dimensions(node_state_dimensions(num_leaves_per_tree, _dimension, _topology()));

	const std::string set_Z_cli_command = "set_" + get_tree_name() + "_Z";
	if (coretools::instances::parameters().exists(set_Z_cli_command)) {
		read_Z_from_file(coretools::instances::parameters().get(set_Z_cli_command), _Z, all_trees,
		                 _dimension);
	}
}

const TNodeStateStorage &TTree::get_Z() const { return _Z; };
TNodeStateStorage &TTree::get_Z() { return _Z; };

void TTree::simulate_Z(size_t tree_index) {
	// A stream of its own, so this draw and the chain's first update are two draws (ADR-0007).
	const TCellUniforms uniforms(run_seed(), TCellStream::node_state_at_start, 0, _dimension);

	for (size_t c = 0; c < _cliques.size(); ++c) {
		auto &clique = _cliques[c];
		_simulation_prepare_cliques(c, clique);
		// Simulation starts from nothing sampled, which is what the old current state also gave
		// here: it was sized but never filled from the storages.
		TCliqueWalkStates current_state(get_number_of_nodes());

		// we sample the roots
		if (ProgramOptions::SIMULATION_NO_Z_INITIALIZATION) { continue; }
		double proba_root = clique.transition_grid().stationary(true);
		coretools::Probability p(proba_root);

		// we can also prepare the queue for the DFS
		std::queue<size_t> node_queue;
		for (const auto root_index_in_tree : this->get_root_nodes()) {
			bool root_state =
			    sample(p, uniforms.at(clique.linear_index_of(_Z, root_index_in_tree)));
			if (root_state) {
				_simulate_one(clique, current_state, tree_index, root_index_in_tree);
			}
			for (const auto child : this->children_of(root_index_in_tree)) {
				if (!this->isLeaf(child)) { node_queue.push(child); }
			} // those are the first children of the tree (children of the roots).
		} // roots done, we go to the internal nodes

		// sampling the internal nodes
		while (!node_queue.empty()) {
			size_t node_index = node_queue.front();
			node_queue.pop();

			// we want to sample the state of the node given its parent (and independently of its
			// children since we haven't sampled them yet).
			std::array<coretools::TSumLogProbability, 2> sum_log;
			clique.calculate_log_prob_parent_to_node(
			    (TypeBinnedBranchLengths)_binned_branch_lengths->value(
			        _topology().branch_index(node_index)),
			    current_state.get(parent_of(node_index)), sum_log);
			bool internal_node_state =
			    sample(sum_log, uniforms.at(clique.linear_index_of(_Z, node_index)));
			if (internal_node_state) {
				_simulate_one(clique, current_state, tree_index, node_index);
			}

			for (size_t child_index : this->children_of(node_index)) {
				if (!this->isLeaf(child_index)) {
					node_queue.push(child_index);
				} // as long as your are not a leaf we can continue sampling Z
			}
		} // internal nodes done, we go to the leaves
	}
}

void TTree::_simulate_one(const TClique &clique, TCliqueWalkStates &current_state,
                          size_t tree_index, size_t node_index_in_tree) {
	auto index_in_leaves_space        = clique.get_start_index_in_leaf_space();
	index_in_leaves_space[tree_index] = node_index_in_tree;
	_Z.insert_one(index_in_leaves_space);
	current_state.set(node_index_in_tree, true);
}
