//
// Created by Marco Visani on 26.06.23.
//

#include "TTree.h"
#include "cli.h"
#include "constants.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Math/TSumLog.h"
#include "coretools/Types/probability.h"
#include "tree/io/node_state_columns.h"
#include "tree/io/read_Z.h"
#include "tree/node_state_draw.h"
#include "tree/node_state_shape.h"

#include <cstddef>
#include <cstdlib>
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
	_alpha_c->initStorage(this, {n_cliques()},
	                      {std::make_shared<coretools::TNamesStrings>(_clique_names)});

	// now we initialize the mu_c_1
	_log_nu_c->initStorage(this, {n_cliques()},
	                       {std::make_shared<coretools::TNamesStrings>(_clique_names)});
	_nu_c.resize(n_cliques());
	for (size_t c = 0; c < n_cliques(); ++c) { _nu_c[c] = std::exp(_log_nu_c->value(c)); }

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
	for (size_t c = 0; c < n_cliques(); ++c) {
		// Draw log_nu[c] ~ Normal(LOG_NU_C, LOG_NU_C_INIT_SD^2) instead of setting every clique to
		// the same constant. Identical initial values would make the MLE that seeds var_log_nu 0,
		// yielding a degenerate prior that freezes log_nu, mean_log_nu and var_log_nu (their
		// proposals would never be accepted).
		const double log_nu_init = coretools::instances::randomGenerator().getNormalRandom(
		    ProgramOptions::LOG_NU_C, ProgramOptions::LOG_NU_C_INIT_SD);
		_log_nu_c->set(c, log_nu_init);
		_alpha_c->set(c, coretools::Probability(ProgramOptions::ALPHA));
		_nu_c[c] = std::exp(_log_nu_c->value(c));
		_set_transition_grid(c, TTransitionGrid(_alpha_c->value(c), _nu_c[c], _grid()));
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
	for (size_t c = 0; c < n_cliques(); ++c) {
		_nu_c[c] = std::exp(_log_nu_c->value(c));
		_set_transition_grid(c, TTransitionGrid(_alpha_c->value(c), _nu_c[c], _grid()));
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
		read_Z_from_file(coretools::instances::parameters().get(set_Z_cli_command), _Z,
		                 node_state_columns(all_trees), _dimension);
	}
}

const TNodeStateStorage &TTree::get_Z() const { return _Z; };
TNodeStateStorage &TTree::get_Z() { return _Z; };

void TTree::simulate_Z() {
	// A stream of its own, so this draw and the chain's first update are two draws (ADR-0007).
	const TCellUniforms uniforms(run_seed(), TCellStream::node_state_at_start, 0, _dimension);

	// The bin each branch sits in. Never asked of a root, which has no branch.
	const auto bin_of = [this](size_t node) { return get_binned_branch_length(node); };

	// One list per clique, committed in one batch below. Not per clique: a sparse node state
	// re-sorts every row and every column of the whole matrix on commit, which is the right cost
	// to pay once over every list at once and the wrong one to pay once per clique (ADR-0006).
	std::vector<std::vector<size_t>> indices_to_insert(n_cliques());

	for (size_t c = 0; c < n_cliques(); ++c) {
		_set_transition_grid(c, TTransitionGrid(_alpha_c->value(c), _nu_c[c], _grid()));

		// This clique's cells of the node state, every node of them. Every node reads the state
		// its parent was given, which the view shows even where the node state could not take
		// the write.
		//
		// The leaves are drawn with the rest. Their block is this tree's tree field (ADR-0005).
		TNodeStateCliqueView nodes(_Z, _topology(), _clique_index(c), _dimension);
		node_state_draw::draw_clique(_topology(), transition_grid(c), bin_of, uniforms, nodes);
		indices_to_insert[c] = nodes.take_deferred_inserts();
	}

	// Cliques are disjoint runs, and a view reads back its own deferred write, so no clique
	// needed another clique's inserts to have landed first.
	_Z.insert_in_Z(indices_to_insert);
}
