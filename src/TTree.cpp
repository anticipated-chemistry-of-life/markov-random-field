//
// Created by Marco Visani on 26.06.23.
//

#include "TTree.h"
#include "TClique.h"
#include "TStorageYVector.h"
#include "TStorageZ.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Math/TSumLog.h"
#include "coretools/Storage/TDimension.h"
#include "coretools/Types/commonWeakTypes.h"
#include "coretools/Types/probability.h"
#include "coretools/algorithms.h"
#include "coretools/devtools.h"
#include "stattools/Updates/TPairIndexSampler.h"

#include <cstddef>
#include <cstdlib>
#include <list>
#include <numeric>
#include <queue>
#include <string>
#include <vector>

TTree::TTree(size_t dimension, const std::string &filename, const std::string &tree_name, TypeParamAlpha *Alpha,
             TypeParamLogNu *LogNu, TypeParamBinBranches *Binned_Branch_Lenghts)
    : _dimension(dimension), _binned_branch_lengths(Binned_Branch_Lenghts), _log_nu_c(LogNu), _alpha_c(Alpha) {

	// tell stattools that these parameters belong to a prior distribution
	this->addPriorParameter({_binned_branch_lengths, _alpha_c, _log_nu_c});

	_load_from_file(filename, tree_name);
}

TTree::~TTree() = default;

void TTree::_initialize_grid_branch_lengths(size_t number_of_branches) {
	// read a, b and K from command-line
	double default_a = 1.0 / 2; // TODO : change to min 1/1000 as per discussion with Dan
	_a               = coretools::instances::parameters().get("a", default_a);
	double default_b = 1.5;
	_b               = coretools::instances::parameters().get("b", default_b);
	_number_of_bins  = coretools::instances::parameters().get("K", 100);

	const size_t max_type = std::numeric_limits<coretools::underlyingType<TypeBinnedBranchLengths>::type>::max();
	if (_number_of_bins >= max_type) {
		UERROR("More bins (", _number_of_bins, ") required than type allows (", max_type,
		       ")! Please decrease K or change type of bins.");
	}
	TypeBinnedBranchLengths::setMax(_number_of_bins - 1);

	// calculate Delta
	_delta = ((double)_b - (double)_a) / (double)_number_of_bins;
}

std::vector<size_t> TTree::_bin_branch_lengths(const std::vector<double> &branch_lengths, bool exclude_root) const {
	std::vector<size_t> binned_branch_lengths;
	binned_branch_lengths.reserve(get_number_of_nodes() - get_number_of_roots());

	double total_branch_length = 0.0;
	for (size_t i = 0; i < branch_lengths.size(); ++i) { // loop over all nodes
		if (exclude_root && _nodes[i].isRoot()) { continue; }
		// find bin
		auto it = std::lower_bound(_grid_branch_lengths.begin(), _grid_branch_lengths.end(), branch_lengths[i]);

		if (it == _grid_branch_lengths.end()) {
			// last bin
			binned_branch_lengths.push_back(_grid_branch_lengths.size() - 1);
			total_branch_length += _grid_branch_lengths.back();
		} else {
			// take the distance between the lower bin and the value and the higher bin and the value
			// and then we take the one that is closer to the value

			auto it_next = it + 1;
			if (it_next == _grid_branch_lengths.end()) {
				// last bin
				binned_branch_lengths.push_back(std::distance(_grid_branch_lengths.begin(), it));
				total_branch_length += _grid_branch_lengths.back();
			} else if (std::abs(branch_lengths[i] - *it) < std::abs(branch_lengths[i] - *it_next)) {
				// take the lower bin
				binned_branch_lengths.push_back(std::distance(_grid_branch_lengths.begin(), it));
				total_branch_length += *it;
			} else {
				// take the higher bin
				binned_branch_lengths.push_back(std::distance(_grid_branch_lengths.begin(), it_next));
				total_branch_length += *it_next;
			}
		}
	}

	// if total branch length is smaller than one, we randomly sample some branch lengths and increase them of one until
	// we get a total branch length of one
	// Adjust total_branch_length to exactly 1.0
	double number_of_branches = double(get_number_of_nodes() - get_number_of_roots());
	while (std::abs(total_branch_length - number_of_branches) > _delta / 2.0) { // TODO: check why _delta / 2.0
		auto idx = coretools::instances::randomGenerator().getRand<size_t>(0, binned_branch_lengths.size() - 1);
		size_t new_bin;
		if (total_branch_length < number_of_branches) {
			// Increase bin index
			if (binned_branch_lengths[idx] >= _grid_branch_lengths.size() - 1) { continue; }
			new_bin = binned_branch_lengths[idx] + 1;

		} else {
			// Decrease bin index
			if (binned_branch_lengths[idx] == 0) { continue; }
			new_bin = binned_branch_lengths[idx] - 1;
		}
		total_branch_length += (_grid_branch_lengths[new_bin] - _grid_branch_lengths[binned_branch_lengths[idx]]);
		binned_branch_lengths[idx] = new_bin;
	}

	return binned_branch_lengths;
}

void TTree::_bin_branch_lengths_from_tree(std::vector<double> &branch_lengths) {
	_initialize_grid_branch_lengths(branch_lengths.size());
	// normalize such that they sum to one

	double sum = 0.0;
	int count  = 0;
	for (const double &branch_length : branch_lengths) {
		if (branch_length <= 0.0) { continue; }
		sum += branch_length;
		++count;
	}
	double average = sum / (double)count;

	for (double &branch_length : branch_lengths) {
		if (branch_length <= 0.0) { continue; }
		branch_length = branch_length / average;
	}

	_grid_branch_lengths.resize(_number_of_bins);
	for (size_t k = 0; k < _number_of_bins; ++k) { _grid_branch_lengths[k] = (_a + _delta * ((double)k + 1.0)); }

	_binned_branch_lengths_from_tree = _bin_branch_lengths(branch_lengths, true);
};

void TTree::_load_from_file(const std::string &filename, const std::string &tree_name) {
	coretools::instances::logfile().listFlush("Reading tree from file '", filename, "' ...");
	coretools::TInputFile file(filename, coretools::FileType::Header);
	this->_tree_name = tree_name;
	if (file.numCols() != 3) {
		UERROR("File '", filename, "' is expected to have 3 columns, but has ", file.numCols(), " !");
	}

	std::vector<double> branch_lengths;

	// read each line of the file
	for (; !file.empty(); file.popFront()) {
		// we read line by line the edge list from the file
		// if the child node (column 0) is not in the tree, we add the node,
		// if the parent is not in the tree, we add the node and set the parent as root
		// as long as it has no parents.
		std::string child  = std::string(file.get(0));
		std::string parent = std::string(file.get(1));
		auto branch_length = file.get<double>(2);
		if (branch_length <= 0.0) { UERROR("You can't have a negative branch length or equal to 0.0 !"); }

		if (!in_tree(child) && !in_tree(parent)) {
			// we add the parent
			_nodes.emplace_back(parent, -1, true);
			branch_lengths.push_back(0.0);
			_node_map[parent] = _nodes.size() - 1;

			// since we just added the parent, we know that it is at the last element of the vector
			size_t parent_index = _nodes.size() - 1;
			_nodes.emplace_back(child, parent_index, false);
			branch_lengths.push_back(branch_length);
			_node_map[child] = _nodes.size() - 1;
			_nodes[parent_index].addChild_index_in_tree(_nodes.size() - 1);
		} else if (!in_tree(child) && in_tree(parent)) {
			size_t parent_index = get_node_index(parent);
			// we add the child node to the tree
			_nodes.emplace_back(child, parent_index, false);
			branch_lengths.push_back(branch_length);
			_node_map[child] = _nodes.size() - 1;
			_nodes[parent_index].addChild_index_in_tree(_nodes.size() - 1);
		} else if (in_tree(child) && !in_tree(parent)) {
			// if the child node was in the tree but not the parent
			// that means that the child was a root and is now becoming
			// a child of the new parent which will be added to the tree
			size_t child_index = get_node_index(child);

			if (!_nodes[child_index].isRoot()) {
				// if the child was not a root and the parent was not in the tree
				// we throw an error because the child already has a parent
				UERROR("Node: '", child, "' has already a parent in the tree. Adding an other parent is not allowed !");
			}
			_nodes[child_index].set_is_root(false);
			branch_lengths[child_index] = branch_length;
			_nodes[child_index].set_parent_index_in_tree(_nodes.size());
			_nodes.emplace_back(parent, -1, true);
			branch_lengths.push_back(0.0);
			_node_map[parent] = _nodes.size() - 1;
			_nodes[_nodes.size() - 1].addChild_index_in_tree(child_index);
		} else {
			size_t child_index  = get_node_index(child);
			size_t parent_index = get_node_index(parent);
			auto node           = _nodes[child_index];
			auto node_parent    = _nodes[parent_index];
			if (!node.isRoot()) {
				// if the child was not a root and the parent was already in the tree
				// we throw an error because the child already has a parent
				UERROR("Node: '", child, "' has already a parent in the tree. Adding an other parent is not allowed !");

			} else {
				// if the child was a root and the parent was already in the tree
				// we set the parent as the parent of the child
				_nodes[child_index].set_is_root(false);
				branch_lengths[child_index] = branch_length;
				_nodes[child_index].set_parent_index_in_tree(parent_index);
				_nodes[parent_index].addChild_index_in_tree(child_index);
			}
		}
	}

	// identify roots and leaves
	this->_leafIndices.resize(_nodes.size(), -1);
	this->_rootIndices.resize(_nodes.size(), -1);
	this->_internalIndices.resize(_nodes.size(), -1);
	this->_internalIndicesWithoutRoots.resize(_nodes.size(), -1);
	this->_leaves_and_internal_nodes_without_roots_indices.resize(_nodes.size(), -1);
	for (auto it = _nodes.begin(); it != _nodes.end(); ++it) {
		if (it->isLeaf()) {
			this->_leafIndices[it - _nodes.begin()] = _leaves.size();
			this->_leaves.push_back(it - _nodes.begin());

			this->_leaves_and_internal_nodes_without_roots_indices[it - _nodes.begin()] =
			    _leaves_and_internal_nodes_without_roots.size();
			this->_leaves_and_internal_nodes_without_roots.push_back(it - _nodes.begin());
		} else {
			this->_internalIndices[it - _nodes.begin()] = _internal_nodes.size();
			this->_internal_nodes.push_back(it - _nodes.begin());
			if (it->isRoot()) {
				_rootIndices[it - _nodes.begin()] = _roots.size();
				_roots.push_back(it - _nodes.begin());
			} else if (it->isInternalNode()) {
				_internalIndicesWithoutRoots[it - _nodes.begin()] = _internal_nodes_without_roots.size();
				_internal_nodes_without_roots.push_back(it - _nodes.begin());

				this->_leaves_and_internal_nodes_without_roots_indices[it - _nodes.begin()] =
				    _leaves_and_internal_nodes_without_roots.size();
				this->_leaves_and_internal_nodes_without_roots.push_back(it - _nodes.begin());
			}
		}
	}

	// binned branches
	_bin_branch_lengths_from_tree(branch_lengths);

	coretools::instances::logfile().done();
	coretools::instances::logfile().conclude("Read ", _nodes.size(), " nodes of which ", _roots.size(),
	                                         " are roots and ", _leaves.size(), " are leaves and ",
	                                         _internal_nodes_without_roots.size(), " are internal nodes.");
}

const TNode &TTree::get_node(const std::string &Id) const {
	auto it = _node_map.find(Id);
	if (it == _node_map.end()) { UERROR("Node '", Id, "' does not exist!"); }
	return _nodes[it->second]; // Retrieve node from vector using the index
}

const TNode &TTree::get_node(size_t index) const { return _nodes[index]; }
bool TTree::isLeaf(size_t index) const { return _nodes[index].isLeaf(); }

size_t TTree::get_node_index(const std::string &Id) const {
	auto it = _node_map.find(Id);
	if (it == _node_map.end()) { UERROR("Node '", Id, "' does not exist in the tree !"); }
	return it->second; // Return the index from the map
}

void TTree::initialize_cliques_and_Z(const std::vector<std::unique_ptr<TTree>> &all_trees) {

	// we initialize the number of leaves we have in each tree
	std::vector<size_t> num_leaves_per_tree(all_trees.size());
	for (size_t i = 0; i < all_trees.size(); ++i) { num_leaves_per_tree[i] = all_trees[i]->get_number_of_leaves(); }

	_initialize_Z(num_leaves_per_tree);
	_initialize_cliques(num_leaves_per_tree, all_trees);
}

[[nodiscard]] std::string TTree::name() const { return "tree"; }

void TTree::initialize() {
	// stattools initialization function
	_alpha_c->initStorage(this, {_cliques.size()}, {std::make_shared<coretools::TNamesIndices>(_cliques.size())});

	// now we initialize the mu_c_1
	_log_nu_c->initStorage(this, {_cliques.size()}, {std::make_shared<coretools::TNamesIndices>(_cliques.size())});
	_nu_c.resize(_cliques.size());
	for (size_t c = 0; c < _cliques.size(); ++c) { _nu_c[c] = std::exp(_log_nu_c->value(c)); }

	// number of branches = number of leaves + number of internal nodes without roots
	_binned_branch_lengths->initStorage(
	    this, {get_number_of_nodes() - get_number_of_roots()},
	    {std::make_shared<coretools::TNamesIndices>(get_number_of_nodes() - get_number_of_roots())});
}

void TTree::guessInitialValues() {
	// TODO: initialize mu's
	for (size_t c = 0; c < _cliques.size(); ++c) {
		_log_nu_c->set(c, -0.1);
		_alpha_c->set(c, coretools::Probability(0.5));
		_nu_c[c] = std::exp(_log_nu_c->value(c));
		_cliques[c].set_lambda(_alpha_c->value(c), _nu_c[c]);
	}

	_set_initial_branch_lengths(false);
}

double TTree::getSumLogPriorDensity(const Storage &) const { DEVERROR("Should never be called"); }
double TTree::getDensity(const Storage &, size_t) const { DEVERROR("Should never be called"); }
double TTree::getLogDensityRatio(const UpdatedStorage &, size_t) const { DEVERROR("Should never be called"); }

void TTree::_set_initial_branch_lengths(bool is_simulation) {
	// overwrite simulated branch length: use branch lengths from tree
	if (_binned_branch_lengths->hasFixedInitialValue() || is_simulation) { // use from simulation
		// translate bin into actual branch lengths
		std::vector<double> vals(_binned_branch_lengths->size());
		for (size_t i = 0; i < _binned_branch_lengths->size(); ++i) {
			vals[i] = (_a + _delta * (_binned_branch_lengths->value(i) + 1.0));
		}

		// normalize such that the average branch length is 1
		double average = std::reduce(vals.begin(), vals.end(), 0.0) / (double)vals.size();
		for (double &val : vals) { val = val / average; }

		// translate back to bin
		auto binned_branch_lengths = _bin_branch_lengths(vals, false);

		// set these values (hack stattools to pretend initial values are not fixed)
		_binned_branch_lengths->fixInitialization(false);
		for (size_t i = 0; i < _binned_branch_lengths->size(); ++i) {
			_binned_branch_lengths->set(i, binned_branch_lengths[i]);

			// we have to do it a second else the oldValue is still the one not normalized.
			// Indeed we first propose the update of the branch length so then in the "update_mu"
			// function, we will have to take the oldValue of the binned branch length.
			_binned_branch_lengths->set(i, binned_branch_lengths[i]);
		}
		_binned_branch_lengths->fixInitialization(true);
	} else { // use from tree
		for (size_t i = 0; i < _binned_branch_lengths->size(); ++i) {
			_binned_branch_lengths->set(i, _binned_branch_lengths_from_tree[i]);
		}
	}
}

void TTree::_simulateUnderPrior(Storage *) {
	using namespace coretools::instances;
	_set_initial_branch_lengths(true);
	for (size_t c = 0; c < _cliques.size(); ++c) {
		_nu_c[c] = std::exp(_log_nu_c->value(c));
		_cliques[c].set_lambda(_alpha_c->value(c), _nu_c[c]);
	}
}

void TTree::_initialize_Z(std::vector<size_t> num_leaves_per_tree) {
	num_leaves_per_tree[_dimension] = this->get_number_of_internal_nodes();

	_Z.initialize_dimensions(num_leaves_per_tree);

	std::string set_Z_cli_command = "set_" + get_tree_name() + "_Z";
	if (coretools::instances::parameters().exists(set_Z_cli_command)) {
		std::string filename = coretools::instances::parameters().get(set_Z_cli_command);
		coretools::TInputFile file(filename, coretools::FileType::Header);

		if (file.numCols() != 4) { UERROR("The file for setting Z must have 4 columns ! "); }
		// read each line of the file
		for (; !file.empty(); file.popFront()) {
			auto linear_index_in_Z_space = file.get<uint32_t>(2);
			bool state                   = file.get<bool>(3);
			if (state) { _Z.insert_one(linear_index_in_Z_space); }
		}
		return;
	}
}

void TTree::_initialize_cliques(const std::vector<size_t> &num_leaves_per_tree,
                                const std::vector<std::unique_ptr<TTree>> &all_trees) {
	// clique of a tree: runs along that dimension
	// the cliques of a tree are can only contain leaves in all trees except the one we are working on.
	_dimension_cliques             = num_leaves_per_tree;
	_dimension_cliques[_dimension] = 1;

	// we then caclulate how many cliques we will have in total for that tree. Which is the product of the number of
	// leaves in each tree except the one we are working on (that is why we set it to 1 before).
	const size_t n_cliques = coretools::containerProduct(_dimension_cliques);

	// calculate increment: product of the number of leaves of all subsequent dimensions
	size_t increment = 1;
	for (size_t i = _dimension + 1; i < all_trees.size(); ++i) { increment *= all_trees[i]->get_number_of_leaves(); }

	// initialize cliques
	for (size_t i = 0; i < n_cliques; ++i) {
		// get start index of each clique in leaves space
		std::vector<size_t> start_index_in_leaves_space = coretools::getSubscripts(i, _dimension_cliques);
		_cliques.emplace_back(start_index_in_leaves_space, _dimension, _nodes.size(), increment);
		_cliques.back().initialize(_a, _delta, _number_of_bins);
	}
}

stattools::TPairIndexSampler TTree::_build_pairs_branch_lengths() const {
	const size_t num_branches = _nodes.size() - get_number_of_roots();
	stattools::TPairIndexSampler sampler(num_branches);
	sampler.sampleIndices();

	return sampler;
}

void TTree::_propose_new_branch_lengths(size_t p1, size_t p2, int val) {
	_binned_branch_lengths->set(p1, (int)_binned_branch_lengths->value(p1) + val);
	_binned_branch_lengths->set(p2, (int)_binned_branch_lengths->value(p2) - val);
}

void TTree::_propose_new_branch_lengths(const stattools::TPairIndexSampler &pairs) {
	for (size_t p = 0; p < pairs.length(); ++p) {
		auto [p1, p2] = pairs.getIndexPair(p);

		const bool both_at_min = _binned_branch_lengths->value(p1) == TypeBinnedBranchLengths::min() &&
		                         _binned_branch_lengths->value(p2) == TypeBinnedBranchLengths::min();
		const bool both_at_max = _binned_branch_lengths->value(p1) == TypeBinnedBranchLengths::max() &&
		                         _binned_branch_lengths->value(p2) == TypeBinnedBranchLengths::max();
		if (both_at_min || both_at_max) {
			// can not do update, both are at boundary
			_propose_new_branch_lengths(p1, p2, 0);
		} else if (_binned_branch_lengths->value(p1) == TypeBinnedBranchLengths::min()) {
			// p1 is at minimum (at the left-most position) -> can only go to the right
			_propose_new_branch_lengths(p1, p2, 1);
		} else if (_binned_branch_lengths->value(p2) == TypeBinnedBranchLengths::min()) {
			// p2 is at minimum (at the left-most position) -> can only go to the right
			_propose_new_branch_lengths(p1, p2, -1);
		} else if (_binned_branch_lengths->value(p1) == TypeBinnedBranchLengths::max()) {
			// p1 is at maximum (at the right-most position) -> can only go to the left
			_propose_new_branch_lengths(p1, p2, -1);
		} else if (_binned_branch_lengths->value(p2) == TypeBinnedBranchLengths::max()) {
			// p2 is at maximum (at the right-most position) -> can only go to the left
			_propose_new_branch_lengths(p1, p2, 1);
		} else {
			// we can choose to go left or right
			bool left = coretools::instances::randomGenerator().pickOneOfTwo(coretools::P(0.5));
			if (left) {
				_propose_new_branch_lengths(p1, p2, -1);
			} else {
				_propose_new_branch_lengths(p1, p2, 1);
			}
		}
	}
};

double TTree::_calculate_likelihood_ratio_branch_length(size_t index_in_binned_branch_length, const TClique &clique,
                                                        const TCurrentState &current_state) const {
	// translate index in binned branch length vector (of size leaves + internal nodes without roots) to index in
	// nodes
	const size_t index_in_tree = _leaves_and_internal_nodes_without_roots[index_in_binned_branch_length];

	// calculate probability of parent to node for old branch length
	double prob_old = clique.calculate_prob_to_parent<false>(
	    index_in_tree, this, _binned_branch_lengths->oldValue(index_in_binned_branch_length), current_state);

	// calculate probability of parent to node for new branch length
	double prob_new = clique.calculate_prob_to_parent<false>(
	    index_in_tree, this, _binned_branch_lengths->value(index_in_binned_branch_length), current_state);

	return prob_new / prob_old;
}

void TTree::_add_to_LL_branch_lengths(size_t c, const TCurrentState &current_state,
                                      std::vector<coretools::TSumLogProbability> &log_sum,
                                      const stattools::TPairIndexSampler &pairs) const {
	const auto &clique = _cliques[c];

	for (size_t p = 0; p < pairs.length(); ++p) { // loop over all possible pairs
		// get index of branches to calculate LL: p1 and p2
		// that index corresponds to the index in fake, concatenated vector of leaves and internal nodes without
		// roots
		auto [p1, p2]   = pairs.getIndexPair(p);
		double ratio_p1 = _calculate_likelihood_ratio_branch_length(p1, clique, current_state);
		double ratio_p2 = _calculate_likelihood_ratio_branch_length(p2, clique, current_state);

		// TODO: make sure we don't have issues with parallelization, probably need pragma once
		log_sum[p].add(ratio_p1);
		log_sum[p].add(ratio_p2);
	}
}

void TTree::_evalute_update_branch_length(std::vector<coretools::TSumLogProbability> &log_sum,
                                          const stattools::TPairIndexSampler &pairs) {
	for (size_t p = 0; p < pairs.length(); ++p) {
		const double LL = log_sum[p].getSum();
		auto [p1, p2]   = pairs.getIndexPair(p);
		const double log_prior_ratio =
		    _binned_branch_lengths->getLogDensityRatio(p1) + _binned_branch_lengths->getLogDensityRatio(p2);
		const double log_H = LL + log_prior_ratio;
		_binned_branch_lengths->acceptOrReject(log_H, coretools::TRange(p1, p2), coretools::TRange());
	}
}

const TStorageZVector &TTree::get_Z() const { return _Z; };
TStorageZVector &TTree::get_Z() { return _Z; };
std::vector<TClique> &TTree::get_cliques() { return _cliques; }
const TClique &TTree::get_clique(const std::vector<size_t> &index_in_leaves_space) const {
	std::vector<size_t> local_index = index_in_leaves_space;
	local_index[_dimension]         = 0; // set to start index
	const size_t ix_clique          = coretools::getLinearIndex(local_index, _dimension_cliques);
	return _cliques[ix_clique];
}

TClique &TTree::get_clique(const std::vector<size_t> &index_in_leaves_space) {
	std::vector<size_t> local_index = index_in_leaves_space;
	local_index[_dimension]         = 0; // set to start index
	const size_t ix_clique        = coretools::getLinearIndex(local_index, _dimension_cliques);
	return _cliques[ix_clique];
}

void TTree::simulate_Z(size_t tree_index) {
	for (size_t c = 0; c < _cliques.size(); ++c) {
		auto &clique = _cliques[c];
		_simulation_prepare_cliques(c, clique);
		TCurrentState current_state(*this, clique.get_increment(), get_number_of_leaves(),
		                            get_number_of_internal_nodes());

		// we sample the roots
		if (SIMULATION_NO_Z_INITIALIZATION) { continue; }
		double proba_root = clique.get_stationary_probability(true, _alpha_c->value(c));
		coretools::Probability p(proba_root);

		// we can also prepare the queue for the DFS
		std::queue<size_t> node_queue;
		for (const auto root_index_in_tree : this->get_root_nodes()) {
			bool root_state = coretools::instances::randomGenerator().pickOneOfTwo(p);
			if (root_state) { _simulate_one(clique, current_state, tree_index, root_index_in_tree); }
			for (const auto child : this->get_node(root_index_in_tree).children_indices_in_tree()) {
				if (!this->isLeaf(child)) { node_queue.push(child); }
			} // those are the first children of the tree (children of the roots).
		} // roots done, we go to the internal nodes

		// sampling the internal nodes
		while (!node_queue.empty()) {
			size_t node_index = node_queue.front();
			node_queue.pop();
			const TNode &node = this->get_node(node_index);

			// we want to sample the state of the node given its parent (and independently of its children since we
			// haven't sampled them yet).
			std::array<coretools::TSumLogProbability, 2> sum_log;
			clique.calculate_log_prob_parent_to_node(node_index,
			                                         (TypeBinnedBranchLengths)_binned_branch_lengths->value(
			                                             _leaves_and_internal_nodes_without_roots_indices[node_index]),
			                                         this, 0, current_state, sum_log);
			bool internal_node_state = sample(sum_log);
			if (internal_node_state) { _simulate_one(clique, current_state, tree_index, node_index); }

			for (size_t child_index : node.children_indices_in_tree()) {
				if (!this->isLeaf(child_index)) {
					node_queue.push(child_index);
				} // as long as your are not a leaf we can continue sampling Z
			}
		} // internal nodes done, we go to the leaves
	}
}

void TTree::_simulation_prepare_cliques(size_t c, TClique &clique) const {
	clique.initialize(this->get_a(), this->get_delta(), this->get_number_of_bins());
	clique.set_lambda(_alpha_c->value(c), _nu_c[c]);
};

void TTree::_simulate_one(const TClique &clique, TCurrentState &current_state, size_t tree_index,
                          size_t node_index_in_tree) {
	auto index_in_leaves_space        = clique.get_start_index_in_leaf_space();
	index_in_leaves_space[tree_index] = this->get_index_within_internal_nodes(node_index_in_tree);
	_Z.insert_one(index_in_leaves_space);
	current_state.set(node_index_in_tree, true);
}
