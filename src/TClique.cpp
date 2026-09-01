//
// Created by madleina on 22.10.24.
//

#include "TClique.h"
#include "constants.h"
#include "coretools/Math/TSumLog.h"
#include "storages/storage_backend.h"
#include "tree/TPhylogeny.h"
#include "tree/TTree.h"
#include <cmath>
#include <cstddef>
#include <vector>

TClique::TClique(const IndexArray &start_index_in_leaves_space, size_t variable_dimension,
                 size_t n_nodes, size_t increment) {
	_start_index_in_leaves_space = start_index_in_leaves_space;
	_variable_dimension          = variable_dimension;
	_n_nodes                     = n_nodes;
	_increment                   = increment;
}

// Both windows open at the same cell and step by the same stride. The field and the leaf block of
// a node state address a leaf pair at the same (row, column) (ADR-0005). Only the clique's own
// dimension differs in extent between the two containers, so the step along the other one is the
// same in both. Each storage turns that subscript into a linear index of its own, and the two
// differ.
TCliqueStates::TCliqueStates(TFieldStorage &Y, TNodeStateStorage &Z, const TClique &clique,
                             const TPhylogeny &topology)
    : _leaves(Y.open_window(clique.first_cell(), topology.n_leaves(), clique.get_increment())),
      _nodes(
          Z.open_window(clique.first_cell(), clique.get_number_of_nodes(), clique.get_increment())),
      _topology(&topology) {}

TCliqueStates TClique::open_states(TFieldStorage &Y, TNodeStateStorage &Z,
                                   const TPhylogeny &topology) const {
	return {Y, Z, *this, topology};
}

TNodeStateStorage::TWindow TClique::open_node_state_window(TNodeStateStorage &Z) const {
	return Z.open_window(first_cell(), _n_nodes, _increment);
}

void TClique::update_Z(std::vector<double> &joint_prob_density, TCliqueStates &states,
                       const TTree *tree, const TCellUniforms &uniforms) const {
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
			                                  states.get(tree->parent_of(index_in_tree)), sum_log);
		}

		// calculate P(child | node = 0) and P(child | node = 1) for all children of node
		_calculate_log_prob_node_to_children(index_in_tree, tree, states, sum_log);

		// sample new state and update Z accordingly. The cell decides which uniform it draws, so
		// the state this node gets does not depend on which thread walked this clique.
		const double log_prob_0 = sum_log[0].getSum();
		const double log_prob_1 = sum_log[1].getSum();
		bool new_state =
		    sample(log_prob_0, log_prob_1, uniforms.at(states.linear_index_in_Z(index_in_tree)));

		if (new_state) {
			joint_prob_density[omp_get_thread_num()] += log_prob_1;
		} else {
			joint_prob_density[omp_get_thread_num()] += log_prob_0;
		}

		// The window writes the cell it already holds in place, and buffers the one it does not,
		// because inserting reallocates a sparse row. A later read on this window sees the new
		// state either way, which is what the parent of this node needs.
		states.set(index_in_tree, new_state);
	}
}

void TClique::initialize_Z_from_children(TCliqueStates &states, const TTree *tree) const {
	// Bottom-up update of Z, as one forward walk. The internal nodes are stored as the non-root
	// block in post-order followed by the roots (ADR-0004), so every node's children are already
	// done by the time it comes up -- leaves before all of them, and each parent after its own
	// children. This used to be a fixed-point loop revisiting a set of not-yet-ready nodes until
	// none remained, which is quadratic in the worst case because nothing about the storage order
	// guaranteed anything.
	for (const size_t node_index : tree->get_internal_nodes()) {
		_set_Z_to_MLE(node_index, states, tree);
	}
}

void TClique::_set_Z_to_MLE(size_t node_index, TCliqueStates &states, const TTree *tree) const {
	std::array<coretools::TSumLogProbability, 2> sum_log;

	_calculate_log_prob_node_to_children(node_index, tree, states, sum_log);

	// sample new state and update Z accordingly
	const double log_prob_0 = sum_log[0].getSum();
	const double log_prob_1 = sum_log[1].getSum();

	bool new_state = log_prob_1 > log_prob_0;
	states.set(node_index, new_state);
}

void TClique::_calculate_log_prob_root(double stationary_0,
                                       std::array<coretools::TSumLogProbability, 2> &sum_log) {
	sum_log[0].add(stationary_0);
	sum_log[1].add(1.0 - stationary_0);
}

void TClique::_calculate_log_prob_node_to_children(
    size_t index_in_tree, const TTree *tree, const TCliqueStates &states,
    std::array<coretools::TSumLogProbability, 2> &sum_log) const {
	const auto &process = transition_grid();
	for (const auto &child_index : tree->children_of(index_in_tree)) {
		// Note: the *previous* bin, because new values were proposed before the loop started.
		// Children come before their parent in post-order, so this reads a state the walk has
		// already assigned -- which the window shows even where it could not write it in place.
		auto bin_length        = tree->get_previous_binned_branch_length(child_index);
		const bool child_state = states.get(child_index);
		for (size_t i = 0; i < 2; ++i) { // loop over possible values (0 or 1) of the node
			sum_log[i].add(process.probability(bin_length, i, child_state));
		}
	}
}

bool sample(coretools::Probability probability_of_one, double uniform) {
	return uniform < probability_of_one;
}

namespace {
/// `1 / (1 + exp(-log_odds))` is the probability of state 1. Log odds of NaN make the comparison
/// false, so a pair of probabilities that says nothing reads as state 0.
///
/// The tails need no branch of their own: the exponential saturates at either end and gives the
/// answer that tail asks for. ADR-0007 says why this is written here rather than taken from the
/// odds-ratio helper it replaces.
///
/// The comparison is written out rather than handed to the probability overload, because a weak
/// probability type rejects the NaN this has to let through.
bool state_one(double log_odds, double uniform) {
	return uniform < 1.0 / (1.0 + std::exp(-log_odds));
}
} // namespace

bool sample(std::array<coretools::TSumLogProbability, 2> &sum_log, double uniform) {
	return state_one(sum_log[1].getSum() - sum_log[0].getSum(), uniform);
}

bool sample(double log_prob_0, double log_prob_1, double uniform) {
	return state_one(log_prob_1 - log_prob_0, uniform);
}

size_t TClique::_get_parent_index(size_t index_in_tree, const TTree *tree) {
	return tree->parent_of(index_in_tree);
}
