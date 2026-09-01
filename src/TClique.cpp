//
// Created by madleina on 22.10.24.
//

#include "TClique.h"
#include "constants.h"
#include "coretools/Math/TSumLog.h"
#include "storages/storage_backend.h"
#include "tree/TPhylogeny.h"
#include "tree/TTree.h"
#include <cstddef>
#include <vector>

TClique::TClique(const IndexArray &start_index_in_leaves_space, size_t variable_dimension,
                 size_t n_nodes, size_t increment) {
	_start_index_in_leaves_space = start_index_in_leaves_space;
	_variable_dimension          = variable_dimension;
	_n_nodes                     = n_nodes;
	_increment                   = increment;
}

bool TClique::state_of(const TFieldStorage &Y, const TInternalStateStorage &Z,
                       const TPhylogeny &topology, size_t node,
                       size_t leaf_index_in_tree_of_last_dim) const {
	const auto cell = cell_of(node, leaf_index_in_tree_of_last_dim);
	if (topology.is_leaf(node)) { return Y.is_one(Y.get_linear_index_in_container_space(cell)); }
	return Z.is_one(Z.get_linear_index_in_container_space(cell));
}

TCliqueWalkStates TClique::read_states(const TFieldStorage &Y, const TInternalStateStorage &Z,
                                       const TPhylogeny &topology) const {
	TCliqueWalkStates walk(topology.n_nodes());
	const size_t column = _start_index_in_leaves_space.back();
	for (size_t node = 0; node < topology.n_nodes(); ++node) {
		walk.set(node, state_of(Y, Z, topology, node, column));
	}
	return walk;
}

void TClique::_write_new_state(TInternalStateStorage &Z, TCliqueWalkStates &walk, size_t node,
                               bool new_state,
                               std::vector<size_t> &linear_indices_in_Z_space_to_insert) const {
	// Mirrors how the field is written in TMarkovField::_set_new_Y: flip a cell the storage already
	// holds, and defer one it does not, so no sparse row is reallocated while threads are reading
	// it. The walk records the new state either way, which is what makes a deferred write visible
	// to the parent that reads it later in this same walk.
	if (walk.get(node) != new_state) {
		const auto cell                      = cell_of(node, _start_index_in_leaves_space.back());
		const size_t linear_index_in_Z_space = Z.get_linear_index_in_container_space(cell);
		if (new_state && !Z.is_stored(linear_index_in_Z_space)) {
			linear_indices_in_Z_space_to_insert.emplace_back(linear_index_in_Z_space);
		} else {
			Z.set_state(linear_index_in_Z_space, new_state);
		}
	}
	walk.set(node, new_state);
}

std::vector<size_t> TClique::update_Z(std::vector<double> &joint_prob_density,
                                      TCliqueWalkStates &walk, TInternalStateStorage &Z,
                                      const TTree *tree) const {
	std::vector<size_t> linear_indices_in_Z_space_to_insert;

	const double stationary_0 = transition_grid().stationary(false);

	for (const auto index_in_tree : tree->get_internal_nodes()) {
		// prepare log probabilities for the two possible states
		std::array<coretools::TSumLogProbability, 2> sum_log;
		if (tree->is_root(index_in_tree)) { // calculate stationary
			_calculate_log_prob_root(stationary_0, sum_log);
		} else { // calculate P(node = 0 | parent) and P(node = 1 | parent)
			// Note: the *previous* bin, because branch lengths are proposed before the loop starts.
			// The parent comes after this node in post-order, so the walk has not touched it and
			// its state is still the one this update started from.
			const auto bin_branch_len = tree->get_previous_binned_branch_length(index_in_tree);
			calculate_log_prob_parent_to_node(bin_branch_len,
			                                  walk.get(tree->parent_of(index_in_tree)), sum_log);
		}

		// calculate P(child | node = 0) and P(child | node = 1) for all children of node
		_calculate_log_prob_node_to_children(index_in_tree, tree, walk, sum_log);

		// sample new state and update Z accordingly
		const double log_prob_0 = sum_log[0].getSum();
		const double log_prob_1 = sum_log[1].getSum();
		bool new_state          = sample(log_prob_0, log_prob_1);

		if (new_state) {
			joint_prob_density[omp_get_thread_num()] += log_prob_1;
		} else {
			joint_prob_density[omp_get_thread_num()] += log_prob_0;
		}

		_write_new_state(Z, walk, index_in_tree, new_state, linear_indices_in_Z_space_to_insert);
	}

	return linear_indices_in_Z_space_to_insert;
}

std::vector<size_t> TClique::initialize_Z_from_children(TCliqueWalkStates &walk,
                                                        TInternalStateStorage &Z,
                                                        const TTree *tree) const {

	// initialise vector that will insert the Z not in parallel
	std::vector<size_t> linear_indices_in_Z_space_to_insert;

	// Bottom-up update of Z, as one forward walk. The internal nodes are stored as the non-root
	// block in post-order followed by the roots (ADR-0004), so every node's children are already
	// done by the time it comes up -- leaves before all of them, and each parent after its own
	// children. This used to be a fixed-point loop revisiting a set of not-yet-ready nodes until
	// none remained, which is quadratic in the worst case because nothing about the storage order
	// guaranteed anything.
	for (const size_t node_index : tree->get_internal_nodes()) {
		_set_Z_to_MLE(node_index, walk, Z, tree, linear_indices_in_Z_space_to_insert);
	}

	return linear_indices_in_Z_space_to_insert;
}

void TClique::_set_Z_to_MLE(size_t node_index, TCliqueWalkStates &walk, TInternalStateStorage &Z,
                            const TTree *tree,
                            std::vector<size_t> &linear_indices_in_Z_space_to_insert) const {
	std::array<coretools::TSumLogProbability, 2> sum_log;

	_calculate_log_prob_node_to_children(node_index, tree, walk, sum_log);

	// sample new state and update Z accordingly
	const double log_prob_0 = sum_log[0].getSum();
	const double log_prob_1 = sum_log[1].getSum();

	bool new_state = log_prob_1 > log_prob_0;
	_write_new_state(Z, walk, node_index, new_state, linear_indices_in_Z_space_to_insert);
}

void TClique::_calculate_log_prob_root(double stationary_0,
                                       std::array<coretools::TSumLogProbability, 2> &sum_log) {
	sum_log[0].add(stationary_0);
	sum_log[1].add(1.0 - stationary_0);
}

void TClique::_calculate_log_prob_node_to_children(
    size_t index_in_tree, const TTree *tree, const TCliqueWalkStates &walk,
    std::array<coretools::TSumLogProbability, 2> &sum_log) const {
	const auto &process = transition_grid();
	for (const auto &child_index : tree->children_of(index_in_tree)) {
		// Note: the *previous* bin, because new values were proposed before the loop started.
		// Children come before their parent in post-order, so this reads a state the walk has
		// already assigned -- which is the whole reason the walk carries them.
		auto bin_length        = tree->get_previous_binned_branch_length(child_index);
		const bool child_state = walk.get(child_index);
		for (size_t i = 0; i < 2; ++i) { // loop over possible values (0 or 1) of the node
			sum_log[i].add(process.probability(bin_length, i, child_state));
		}
	}
}

bool sample(std::array<coretools::TSumLogProbability, 2> &sum_log) {
	const double log_Q = sum_log[1].getSum() - sum_log[0].getSum();
	return coretools::TAcceptOddsRatio::accept(log_Q);
}

bool sample(double log_prob_0, double log_prob_1) {
	const double log_Q = log_prob_1 - log_prob_0;
	return coretools::TAcceptOddsRatio::accept(log_Q);
}

size_t TClique::_get_parent_index(size_t index_in_tree, const TTree *tree) {
	return tree->parent_of(index_in_tree);
}
