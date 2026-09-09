//
// What the link says about one leaf of one tree's node state.
//
// A tree draws every node of its node state, leaves included. Every node above the leaf block is
// conditionally independent of everything outside its own tree, and a leaf is not: it is one of the
// two tree field cells the link reconciles into the field cell at that leaf pair (ADR-0005). So a
// leaf takes one term the rest of the walk does not, and that term is what this file answers.
//
// The walk asks it through a seam (tree/node_state_walk.h), and the tree only forwards what its
// caller bound. That is what keeps a tree's update free of the field, of the other tree and of the
// error probability -- and what lets the walk run in a test with none of the three.
//
// It reads two storages and writes neither. Every thread of a tree's update asks one of these, and
// the field and the other tree's node state stand still for the whole of that update: each is drawn
// by an update of its own, and the three run one after another.
//

#pragma once

#include "constants.h"
#include "field/TFieldMath.h"
#include "storages/storage_concepts.h"
#include "tree/clique/TCliqueView.h"

#include <array>
#include <cstddef>

namespace tree_field_link {

// The link reads one tree field cell per tree, and the two are told apart by dimension. A third
// tree would be a third state in the bucket, so this says so here rather than looping over a
// dimension count that cannot change.
static_assert(NUMBER_OF_TREES == 2,
              "The tree field link is written for one species tree and one molecule tree.");

/// The link at the leaves of one clique of one tree.
///
/// A clique carries a leaf in every dimension but its own tree's, so its index plus a leaf names
/// one leaf pair -- and a leaf's node index is its index in leaf space (ADR-0004), so the same
/// subscript addresses the field, this tree's node state and the other tree's. The linear indices
/// differ, because the container shapes do.
template<field_math::LinkPolicy Policy, BinaryStorage Field, BinaryStorage NodeState>
class TCliqueLink {
private:
	const Field *_Y                    = nullptr;
	const NodeState *_other_tree_field = nullptr;
	field_math::TErrorProbability _omega;
	IndexArray _clique_index{};
	size_t _dimension = 0;

public:
	TCliqueLink(const Field &Y, const NodeState &other_tree_field,
	            const field_math::TErrorProbability &omega, const IndexArray &clique_index,
	            size_t dimension)
	    : _Y(&Y), _other_tree_field(&other_tree_field), _omega(omega), _clique_index(clique_index),
	      _dimension(dimension) {}

	/// `{ P(Y | this leaf = 0), P(Y | this leaf = 1) }` at the leaf pair this leaf occupies.
	///
	/// The field cell and the other tree's cell are read and never written. The walk multiplies
	/// this into the term its own process gives the leaf, which is the leaf's full conditional.
	[[nodiscard]] std::array<double, 2> prob_of_leaf_states(size_t leaf) const {
		const IndexArray cell = clique_cell(_clique_index, _dimension, leaf);
		const bool y          = _Y->is_one(_Y->get_linear_index_in_container_space(cell));
		const bool other =
		    _other_tree_field->is_one(_other_tree_field->get_linear_index_in_container_space(cell));

		std::array<double, 2> probability{};
		for (size_t state = 0; state < 2; ++state) {
			const bool mine = state != 0;
			// The link takes the species tree field first. Which of the two this tree is, is what
			// its dimension says. It matters even while one error probability serves both, because
			// the link is already written for one per tree.
			const double link = _dimension == 0 ? Policy::prob_y_is_one(mine, other, _omega)
			                                    : Policy::prob_y_is_one(other, mine, _omega);
			probability[state] = y ? link : 1.0 - link;
		}
		return probability;
	}
};

/// The field, the other tree's node state and the error probability, bound for the length of one
/// tree's update. It hands out one link per clique.
///
/// One of these is built per tree per iteration, outside the parallel region, and handed to that
/// tree. The tree names it by concept alone.
template<field_math::LinkPolicy Policy, BinaryStorage Field, BinaryStorage NodeState>
class TLeafLinks {
private:
	const Field *_Y                    = nullptr;
	const NodeState *_other_tree_field = nullptr;
	field_math::TErrorProbability _omega;

public:
	TLeafLinks(const Field &Y, const NodeState &other_tree_field,
	           const field_math::TErrorProbability &omega)
	    : _Y(&Y), _other_tree_field(&other_tree_field), _omega(omega) {}

	/// The link at the leaves of the clique `clique_index` of the tree that owns `dimension`.
	///
	/// The tree passes both, so the dimension the link reads the clique with is the one the tree
	/// walks it with, rather than a second copy of the same number.
	[[nodiscard]] TCliqueLink<Policy, Field, NodeState> for_clique(const IndexArray &clique_index,
	                                                              size_t dimension) const {
		return {*_Y, *_other_tree_field, _omega, clique_index, dimension};
	}
};

} // namespace tree_field_link
