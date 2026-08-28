//
// Created by madleina on 22.10.24.
//

#include "TClique.h"
#include "TCurrentState.h"
#include "constants.h"
#include "coretools/Math/TSumLog.h"
#include "storages/storage_backend.h"
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

TCurrentState TClique::create_current_state(const TFieldStorage &Y, const TInternalStateStorage &Z,
                                            const TPhylogeny &topology) {
	TCurrentState current_state(topology, this->_increment, topology.n_leaves(),
	                            topology.n_internal_nodes());
	current_state.fill(_start_index_in_leaves_space, Y, Z);

	return current_state;
}

std::vector<size_t> TClique::update_Z(std::vector<double> &joint_prob_density,
                                      TCurrentState &current_state, TInternalStateStorage &Z,
                                      const TTree *tree) const {
	std::vector<size_t> linear_indices_in_Z_space_to_insert;

	const double stationary_0 = transition_grid().stationary(false);

	for (const auto index_in_tree : tree->get_internal_nodes()) {
		// prepare log probabilities for the two possible states
		std::array<coretools::TSumLogProbability, 2> sum_log;
		if (tree->is_root(index_in_tree)) { // calculate stationary
			_calculate_log_prob_root(stationary_0, sum_log);
		} else { // calculate P(node = 0 | parent) and P(node = 1 | parent)
			// note: for compatibility with update of Y, we need to pass
			// leaf_index_in_tree_of_last_dim, but this doesn't matter for this update (just pass 0)
			// Note: the *previous* bin, because branch lengths are proposed before the loop starts
			const auto bin_branch_len = tree->get_previous_binned_branch_length(index_in_tree);
			calculate_log_prob_parent_to_node(index_in_tree, bin_branch_len, tree, 0, current_state,
			                                  sum_log);
		}

		// calculate P(child | node = 0) and P(child | node = 1) for all children of node
		_calculate_log_prob_node_to_children(index_in_tree, tree, current_state, sum_log);

		// sample new state and update Z accordingly
		const double log_prob_0 = sum_log[0].getSum();
		const double log_prob_1 = sum_log[1].getSum();
		bool new_state          = sample(log_prob_0, log_prob_1);

		if (new_state) {
			joint_prob_density[omp_get_thread_num()] += log_prob_1;
		} else {
			joint_prob_density[omp_get_thread_num()] += log_prob_0;
		}

		_update_current_state(Z, current_state, index_in_tree, new_state,
		                      linear_indices_in_Z_space_to_insert, tree);
	}

	return linear_indices_in_Z_space_to_insert;
}

void TClique::_update_current_state(TInternalStateStorage &Z, TCurrentState &current_state,
                                    size_t index_in_tree, bool new_state,
                                    std::vector<size_t> &linear_indices_in_Z_space_to_insert,
                                    const TTree * /*tree*/) const {
	// The Z slot of current_state holds the linear index in Z space of this node (filled by
	// TInternalStateStorage::fill_current_state). This mirrors how Y is updated in
	// TMarkovField::_set_new_Y: in-place state flips for cells that already exist, and deferred
	// bulk insertion for new cells (so the shared sparse matrix is not reallocated in parallel).
	const size_t linear_index_in_Z_space = current_state.get_linear_index_in_storage(index_in_tree);
	const bool cur_state                 = current_state.get(index_in_tree);
	const bool exists                    = current_state.exists_in_storage(index_in_tree);

	if (cur_state && !new_state) { // 1 -> 0: cell exists -> flip state in place
		Z.set_state(linear_index_in_Z_space, false);
	} else if (!cur_state && new_state) { // 0 -> 1
		if (exists) {                     // already stored -> flip state in place
			Z.set_state(linear_index_in_Z_space, true);
		} else { // not stored yet -> defer the insert until after the parallel region
			linear_indices_in_Z_space_to_insert.emplace_back(linear_index_in_Z_space);
		}
	}
	current_state.set(index_in_tree, new_state);
}

std::vector<size_t> TClique::initialize_Z_from_children(TCurrentState &current_state,
                                                        TInternalStateStorage &Z,
                                                        const TTree *tree) const {

	// initialise vector that will insert the Z not in parallel
	std::vector<size_t> linear_indices_in_Z_space_to_insert;

	// Bottom-up update of Z, as one forward walk. The internal nodes are stored as the non-root
	// block in post-order followed by the roots (ADR-0004), so every node's children are already
	// done by the time it comes up -- leaves before all of them, and each parent after its own
	// children. This used to be a fixed-point loop re-sweeping a set of not-yet-ready nodes until
	// none remained, which is quadratic in the worst case because nothing about the storage order
	// guaranteed anything.
	for (const size_t node_index : tree->get_internal_nodes()) {
		_set_Z_to_MLE(node_index, current_state, Z, tree, linear_indices_in_Z_space_to_insert);
	}

	return linear_indices_in_Z_space_to_insert;
}

void TClique::_set_Z_to_MLE(size_t node_index, TCurrentState &current_state,
                            TInternalStateStorage &Z, const TTree *tree,
                            std::vector<size_t> &linear_indices_in_Z_space_to_insert) const {
	std::array<coretools::TSumLogProbability, 2> sum_log;

	_calculate_log_prob_node_to_children(node_index, tree, current_state, sum_log);

	// sample new state and update Z accordingly
	const double log_prob_0 = sum_log[0].getSum();
	const double log_prob_1 = sum_log[1].getSum();

	bool new_state = log_prob_1 > log_prob_0;
	_update_current_state(Z, current_state, node_index, new_state,
	                      linear_indices_in_Z_space_to_insert, tree);
}

void TClique::_calculate_log_prob_root(double stationary_0,
                                       std::array<coretools::TSumLogProbability, 2> &sum_log) {
	sum_log[0].add(stationary_0);
	sum_log[1].add(1.0 - stationary_0);
}

void TClique::_calculate_log_prob_node_to_children(
    size_t index_in_tree, const TTree *tree, const TCurrentState &current_state,
    std::array<coretools::TSumLogProbability, 2> &sum_log) const {
	const auto &process = transition_grid();
	for (const auto &child_index : tree->children_of(index_in_tree)) {
		// Note: the *previous* bin, because new values were proposed before the loop started
		auto bin_length        = tree->get_previous_binned_branch_length(child_index);
		const bool child_state = current_state.get(child_index);
		for (size_t i = 0; i < 2; ++i) { // loop over possible values (0 or 1) of the node
			sum_log[i].add(process.probability(bin_length, i, child_state));
		}
	}
}

void TClique::update_counter_leaves_state_1(bool new_state, bool old_state) {
	const int diff = (int)new_state - (int)old_state;
	update_counter_leaves_state_1(diff);
}
void TClique::update_counter_leaves_state_1(int difference) {
	assert((int)_counter_leaves_state_1 + difference >= 0);
	_counter_leaves_state_1 += difference;
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
