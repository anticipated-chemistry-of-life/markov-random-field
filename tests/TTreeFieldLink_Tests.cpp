//
// The seam a tree's leaves are drawn through.
//
// The link's arithmetic is pure and tested against brute force (TFieldMath_Tests.cpp), and the
// walk that consumes this is tested against a link a test writes by hand (TNodeStateWalk_Tests.cpp).
// What is left belongs to this file alone: that the two cells it reads are the ones at the leaf
// pair the leaf occupies, that the tree it stands for is the one its dimension names, and that it
// hands the walk `{ P(Y | leaf = 0), P(Y | leaf = 1) }` and not something proportional to it.
//
// Those questions are all about indexing, so the suite runs over every pairing of the two storages
// and over tree pairs whose shapes differ -- a chain has a single leaf, which makes a container one
// cell wide in that dimension, and that is where an index property stops being a coincidence.
//

#include "backend_pairings.h"
#include "constants.h"
#include "field/TFieldMath.h"
#include "field/link_backend.h"
#include "field/tree_field_link.h"
#include "tree/TPhylogeny.h"
#include "gtest/gtest.h"

#include <array>
#include <cstddef>
#include <string>

namespace {

using backends::AllBackends;
using backends::field_shape;
using backends::make_storage;
using backends::molecule_shape;
using backends::seed_ones;
using backends::species_shape;
using backends::tree_pairs;

/// Small enough that the link never runs degenerate, and recognisable in a failure message.
constexpr double OMEGA = 0.125;

/// The clique of one tree named by a leaf of the other: a leaf in every dimension but this tree's
/// own, which carries a 0 (tree/clique/TCliqueView.h).
IndexArray clique_of(size_t dimension, size_t other_leaf) {
	IndexArray clique{};
	clique[dimension]     = 0;
	clique[1 - dimension] = other_leaf;
	return clique;
}

template<typename Backends> class TreeFieldLink : public ::testing::Test {
public:
	using Field     = typename Backends::field;
	using NodeState = typename Backends::node_state;
};

TYPED_TEST_SUITE(TreeFieldLink, AllBackends);

/// At every leaf of every clique of both trees, the link is the link table read at that leaf pair:
/// the field cell there, against this tree's two states and the other tree's state.
///
/// The expectation is written out from `prob_y_is_one` and the two cells the test looked up itself,
/// so a link that read the wrong cell -- the transposed leaf pair, the wrong clique, or its own
/// tree's node state instead of the other's -- gives a different number.
TYPED_TEST(TreeFieldLink, is_the_link_table_at_its_own_leaf_pair) {
	using Field     = typename TestFixture::Field;
	using NodeState = typename TestFixture::NodeState;

	const field_math::TErrorProbability omega(OMEGA);
	for (const auto &pair : tree_pairs()) {
		SCOPED_TRACE(pair.name);
		auto Y          = make_storage<Field>(field_shape(pair));
		auto Z_species  = make_storage<NodeState>(species_shape(pair));
		auto Z_molecule = make_storage<NodeState>(molecule_shape(pair));
		seed_ones(Y, 1);
		seed_ones(Z_species, 2);
		seed_ones(Z_molecule, 3);

		// One source of links per tree, each carrying the *other* tree's node state.
		const tree_field_link::TLeafLinks<TLinkPolicy, Field, NodeState> for_species(Y, Z_molecule,
		                                                                            omega);
		const tree_field_link::TLeafLinks<TLinkPolicy, Field, NodeState> for_molecule(Y, Z_species,
		                                                                             omega);

		for (size_t s = 0; s < pair.species.n_leaves(); ++s) {
			for (size_t m = 0; m < pair.molecule.n_leaves(); ++m) {
				SCOPED_TRACE("leaf pair " + std::to_string(s) + "," + std::to_string(m));
				const bool y = Y.is_one(Y.get_linear_index_in_container_space({s, m}));
				const bool z_s =
				    Z_species.is_one(Z_species.get_linear_index_in_container_space({s, m}));
				const bool z_m =
				    Z_molecule.is_one(Z_molecule.get_linear_index_in_container_space({s, m}));

				// The species tree's clique is named by the molecule leaf, and vice versa.
				const auto species_link  = for_species.for_clique(clique_of(0, m), 0);
				const auto molecule_link = for_molecule.for_clique(clique_of(1, s), 1);

				const std::array<double, 2> got_species = species_link.prob_of_leaf_states(s);
				const std::array<double, 2> got_molecule = molecule_link.prob_of_leaf_states(m);

				for (size_t state = 0; state < 2; ++state) {
					const bool mine = state != 0;
					// The species tree's leaf varies the first argument, the molecule tree's the
					// second, and each reads the other tree's cell as it stands.
					const double link_species =
					    TLinkPolicy::prob_y_is_one(mine, z_m, omega);
					const double link_molecule =
					    TLinkPolicy::prob_y_is_one(z_s, mine, omega);
					EXPECT_DOUBLE_EQ(got_species[state], y ? link_species : 1.0 - link_species)
					    << "species leaf state " << state;
					EXPECT_DOUBLE_EQ(got_molecule[state], y ? link_molecule : 1.0 - link_molecule)
					    << "molecule leaf state " << state;
				}
			}
		}
	}
}

/// What a link reports is a probability of the field cell, and not a weight with a constant
/// dropped. Two things say so, and both are checked here.
///
/// Every value is strictly inside (0, 1), which is what the walk needs: it sums logs, so a zero
/// would rule a leaf state out rather than make it unlikely.
///
/// And the value at a cell whose field state is one, plus the value at the same cell with the field
/// state flipped, is exactly one. The two are `P(Y = 1 | leaf)` and `P(Y = 0 | leaf)`, so they must
/// be. A pair of unnormalised weights would not be.
TYPED_TEST(TreeFieldLink, reports_a_probability_of_the_field_cell_under_each_leaf_state) {
	using Field     = typename TestFixture::Field;
	using NodeState = typename TestFixture::NodeState;

	const field_math::TErrorProbability omega(OMEGA);
	const auto &pair = tree_pairs().front();
	auto Z_species   = make_storage<NodeState>(species_shape(pair));
	auto Z_molecule  = make_storage<NodeState>(molecule_shape(pair));
	seed_ones(Z_species, 2);
	seed_ones(Z_molecule, 3);

	// The same configuration twice, once with an empty field and once with a full one.
	auto empty_field = make_storage<Field>(field_shape(pair));
	auto full_field  = make_storage<Field>(field_shape(pair));
	for (size_t i = 0; i < full_field.total_size_of_container_space(); ++i) {
		full_field.insert_one(i);
	}

	const tree_field_link::TLeafLinks<TLinkPolicy, Field, NodeState> at_zero(empty_field,
	                                                                        Z_molecule, omega);
	const tree_field_link::TLeafLinks<TLinkPolicy, Field, NodeState> at_one(full_field, Z_molecule,
	                                                                       omega);

	for (size_t m = 0; m < pair.molecule.n_leaves(); ++m) {
		const auto zero = at_zero.for_clique(clique_of(0, m), 0);
		const auto one  = at_one.for_clique(clique_of(0, m), 0);
		for (size_t s = 0; s < pair.species.n_leaves(); ++s) {
			SCOPED_TRACE("leaf pair " + std::to_string(s) + "," + std::to_string(m));
			const std::array<double, 2> below_zero = zero.prob_of_leaf_states(s);
			const std::array<double, 2> below_one  = one.prob_of_leaf_states(s);
			for (size_t state = 0; state < 2; ++state) {
				EXPECT_GT(below_zero[state], 0.0);
				EXPECT_LT(below_zero[state], 1.0);
				EXPECT_DOUBLE_EQ(below_zero[state] + below_one[state], 1.0);
			}
		}
	}
}

/// A leaf's state moves the link in the direction the AND says it should: a leaf at one makes a
/// field cell at one more likely than a leaf at zero does, whatever the other tree holds.
TYPED_TEST(TreeFieldLink, a_leaf_at_one_makes_the_field_cell_at_one_more_likely) {
	using Field     = typename TestFixture::Field;
	using NodeState = typename TestFixture::NodeState;

	const field_math::TErrorProbability omega(OMEGA);
	const auto &pair = tree_pairs().front();
	auto Y           = make_storage<Field>(field_shape(pair));
	auto Z_species   = make_storage<NodeState>(species_shape(pair));
	auto Z_molecule  = make_storage<NodeState>(molecule_shape(pair));
	seed_ones(Y, 1);
	seed_ones(Z_species, 2);
	seed_ones(Z_molecule, 3);

	const tree_field_link::TLeafLinks<TLinkPolicy, Field, NodeState> for_species(Y, Z_molecule,
	                                                                            omega);
	for (size_t m = 0; m < pair.molecule.n_leaves(); ++m) {
		const auto link = for_species.for_clique(clique_of(0, m), 0);
		for (size_t s = 0; s < pair.species.n_leaves(); ++s) {
			SCOPED_TRACE("leaf pair " + std::to_string(s) + "," + std::to_string(m));
			const std::array<double, 2> below = link.prob_of_leaf_states(s);
			if (Y.is_one(Y.get_linear_index_in_container_space({s, m}))) {
				EXPECT_GT(below[1], below[0]) << "a field cell at one";
			} else {
				EXPECT_LT(below[1], below[0]) << "a field cell at zero";
			}
		}
	}
}

} // namespace
