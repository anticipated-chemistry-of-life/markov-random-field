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
#include "tree/io/read_Z.h"

#include <coretools/algorithms.h>
#include <cstddef>
#include <cstdlib>
#include <fmt/base.h>
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
	_alpha_c->initStorage(this, {_n_cliques},
	                      {std::make_shared<coretools::TNamesStrings>(_clique_names)});

	// now we initialize the mu_c_1
	_log_nu_c->initStorage(this, {_n_cliques},
	                       {std::make_shared<coretools::TNamesStrings>(_clique_names)});
	_nu_c.resize(_n_cliques);
	for (size_t c = 0; c < _n_cliques; ++c) { _nu_c[c] = std::exp(_log_nu_c->value(c)); }

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
	for (size_t c = 0; c < _n_cliques; ++c) {
		// Draw log_nu[c] ~ Normal(LOG_NU_C, LOG_NU_C_INIT_SD^2) instead of setting every clique to
		// the same constant. Identical initial values would make the MLE that seeds var_log_nu 0,
		// yielding a degenerate prior that freezes log_nu, mean_log_nu and var_log_nu (their
		// proposals would never be accepted).
		const double log_nu_init = coretools::instances::randomGenerator().getNormalRandom(
		    ProgramOptions::LOG_NU_C, ProgramOptions::LOG_NU_C_INIT_SD);
		_log_nu_c->set(c, log_nu_init);
		_alpha_c->set(c, coretools::Probability(ProgramOptions::ALPHA));
		_nu_c[c]                       = std::exp(_log_nu_c->value(c));
		_transition_grid_per_clique[c] = TTransitionGrid(_alpha_c->value(c), _nu_c[c], _grid());
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
	for (size_t c = 0; c < _n_cliques; ++c) {
		_nu_c[c]                       = std::exp(_log_nu_c->value(c));
		_transition_grid_per_clique[c] = TTransitionGrid(_alpha_c->value(c), _nu_c[c], _grid());
	}
}

void TTree::_initialize_Z(IndexArray num_leaves_per_tree,
                          [[maybe_unused]] const std::vector<std::unique_ptr<TTree>> &all_trees) {
	num_leaves_per_tree[_dimension] = this->get_number_of_nodes();

	_Z.initialize_dimensions(num_leaves_per_tree);

	const std::string set_Z_cli_command = "set_" + get_tree_name() + "_Z";
	if (coretools::instances::parameters().exists(set_Z_cli_command)) {
		throw std::runtime_error("set_" + get_tree_name() + "_Z is not supported yet");
		// read_Z_from_file(coretools::instances::parameters().get(set_Z_cli_command), _Z,
		// all_trees, _dimension);
	}
}

const TTreeStateStorage &TTree::get_Z() const { return _Z; };
TTreeStateStorage &TTree::get_Z() { return _Z; };

void TTree::calculate_log_prob_parent_to_node(
    size_t index_in_tree, size_t clique_index, TypeBinnedBranchLengths binned_branch_length,
    std::array<coretools::TSumLogProbability, 2> &sum_log) const {
	const auto &process       = _transition_grid_per_clique.at(clique_index);
	const size_t parent_index = parent_of(index_in_tree);
	IndexArray index          = coretools::getSubscriptsAsArray(clique_index, _dimension_cliques);
	index[_dimension]         = parent_index;
	for (size_t i = 0; i < 2; ++i) {
		const bool parent_state = _Z.is_one(index).is_one;
		sum_log[i].add(process.probability(binned_branch_length, parent_state, i));
	}
};

void TTree::simulate_Z() {
	IndexArray index{};
	for (size_t c = 0; c < _n_cliques; ++c) {
		_simulation_prepare_cliques(c);
		index = coretools::getSubscriptsAsArray(c, _dimension_cliques);

		// we sample the roots
		if (ProgramOptions::SIMULATION_NO_Z_INITIALIZATION) { continue; }
		double proba_root = _transition_grid_per_clique[c].stationary(true);
		coretools::Probability p(proba_root);

		// we can also prepare the queue for the DFS
		std::queue<size_t> node_queue;
		for (const auto root_index_in_tree : this->get_root_nodes()) {
			bool root_state = coretools::instances::randomGenerator().pickOneOfTwo(p);
			if (root_state) {
				index[_dimension] = root_index_in_tree;
				_Z.insert_one(index);
			}
			for (const auto child : this->children_of(root_index_in_tree)) {
				node_queue.push(child);
			} // those are the first children of the tree (children of the roots).
		} // roots done, we go to the internal nodes

		// sampling the internal nodes
		while (!node_queue.empty()) {
			size_t node_index = node_queue.front();
			node_queue.pop();

			// we want to sample the state of the node given its parent (and independently of its
			// children since we haven't sampled them yet).
			std::array<coretools::TSumLogProbability, 2> sum_log;
			calculate_log_prob_parent_to_node(
			    node_index, c,
			    (TypeBinnedBranchLengths)_binned_branch_lengths->value(
			        _topology().branch_index(node_index)),
			    sum_log);
			bool internal_node_state = sample(sum_log);
			if (internal_node_state) {
				index[_dimension] = node_index;
				_Z.insert_one(index);
			}

			for (size_t child_index : this->children_of(node_index)) {
				node_queue.push(child_index);
			}
		}
	}
}

size_t TTree::clique_index(const IndexArray &index_in_leaves_space) const {
	IndexArray index  = index_in_leaves_space;
	index[_dimension] = 1;
	return coretools::containerProduct(index);
}

void TTree::update_Z_clique(size_t clique_index, std::vector<double> &joint_prob_density,
                            const TFieldStorage &Y) {

	const double stationary_0 = _transition_grid_per_clique.at(clique_index).stationary(false);
	for (size_t index_in_tree = 0; index_in_tree < get_number_of_nodes(); ++index_in_tree) {
		// prepare log probabilities for the two possible states
		std::array<coretools::TSumLogProbability, 2> sum_log;

		if (this->is_root(index_in_tree)) { // calculate stationary
			_calculate_log_prob_root(stationary_0, sum_log);
		} else { // calculate P(node = 0 | parent) and P(node = 1 | parent)
			// note: for compatibility with update of Y, we need to pass
			// leaf_index_in_tree_of_last_dim, but this doesn't matter for this update (just pass 0)
			// Note: the *previous* bin, because branch lengths are proposed before the loop starts
			const auto bin_branch_len = this->get_previous_binned_branch_length(index_in_tree);
			calculate_log_prob_parent_to_node(index_in_tree, clique_index, bin_branch_len, sum_log);
		}

		// calculate P(child | node = 0) and P(child | node = 1) for all children of node
		if (!this->isLeaf(index_in_tree)) {
			_calculate_log_prob_node_to_children(index_in_tree, clique_index, sum_log);
		}
		if (this->isLeaf(index_in_tree)) {
			_calculate_log_prob_leaf_to_Y(index_in_tree, clique_index, sum_log, Y);
		}

		// sample new state and update Z accordingly
		const double log_prob_0 = sum_log[0].getSum();
		const double log_prob_1 = sum_log[1].getSum();
		bool new_state          = sample(log_prob_0, log_prob_1);

		if (new_state) {
			joint_prob_density[omp_get_thread_num()] += log_prob_1;
		} else {
			joint_prob_density[omp_get_thread_num()] += log_prob_0;
		}

		_update_state();
	}
}

void TTree::_calculate_log_prob_root(double stationary_0,
                                     std::array<coretools::TSumLogProbability, 2> &sum_log) {
	sum_log[0].add(stationary_0);
	sum_log[1].add(1.0 - stationary_0);
}

void TTree::_calculate_log_prob_node_to_children(
    size_t index_in_tree, size_t clique_index,
    std::array<coretools::TSumLogProbability, 2> &sum_log) const {
	const auto &process = _transition_grid_per_clique.at(clique_index);
	for (const auto &child_index : children_of(index_in_tree)) {
		// Note: the *previous* bin, because new values were proposed before the loop started
		auto bin_length        = get_previous_binned_branch_length(child_index);
		auto index             = coretools::getSubscriptsAsArray(clique_index, _dimension_cliques);
		index[_dimension]      = child_index;
		const bool child_state = _Z.is_one(child_index).is_one;
		for (size_t i = 0; i < 2; ++i) { // loop over possible values (0 or 1) of the node
			sum_log[i].add(process.probability(bin_length, i, child_state));
		}
	}
}

void TTree::_calculate_log_prob_leaf_to_Y(size_t index_in_tree, size_t clique_index,
                                          std::array<coretools::TSumLogProbability, 2> &sum_log,
                                          const TFieldStorage &Y) {
	auto index         = coretools::getSubscriptsAsArray(clique_index, _dimension_cliques);
	index[_dimension]  = index_in_tree;
	// since the trees are in topological order, the index of the Y state is exactly the same as the
	// leaf index of both trees.
	const bool y_state = Y.is_one(index).is_one;
	for (size_t i = 0; i < 2; ++i) { // loop over possible values (0 or 1) of the node
		sum_log[i].add(0.0);
	}
}
