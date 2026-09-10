//
// Properties of the clique <-> cell correspondence.
//
// Before TCliqueSpace, the two directions were two methods on TTree written in opposite
// conventions: one read a clique number as a row-major subscript, the other walked a column-major
// stride. They agreed, because clique space has one dimension above 1 when there are two trees.
// Neither could be asked from a test, because both were private to a class that stands up three
// live stattools parameters.
//
// The round trip over three dimensions below is the case that tells the two conventions apart. It
// is a property of this module and not a claim about the model: two trees is what the block update
// and the eight-state table are written for. See ADR-0011.
//

#include "tree/clique/TCliqueSpace.h"
#include "gtest/gtest.h"

#include <array>
#include <cstddef>
#include <stdexcept>

namespace {

constexpr size_t SPECIES   = 0;
constexpr size_t MOLECULES = 1;

// The two trees, at the shape the fixtures run: a clique of the species tree is named by a molecule
// leaf, and a clique of the molecule tree by a species leaf.
constexpr size_t N_SPECIES_LEAVES   = 7;
constexpr size_t N_MOLECULE_LEAVES  = 5;

TCliqueSpace<> species_tree_space() {
	return TCliqueSpace<>({N_SPECIES_LEAVES, N_MOLECULE_LEAVES}, SPECIES);
}

TCliqueSpace<> molecule_tree_space() {
	return TCliqueSpace<>({N_SPECIES_LEAVES, N_MOLECULE_LEAVES}, MOLECULES);
}

//-----------------------------------
// How many cliques a tree has
//-----------------------------------

TEST(CliqueSpace, counts_one_clique_per_leaf_of_every_other_tree) {
	// A clique runs along the owning tree, so that tree's leaf count does not enter the product.
	EXPECT_EQ(species_tree_space().n_cliques(), N_MOLECULE_LEAVES);
	EXPECT_EQ(molecule_tree_space().n_cliques(), N_SPECIES_LEAVES);
}

TEST(CliqueSpace, counts_the_product_over_three_dimensions) {
	const TCliqueSpace<3> space({7, 5, 3}, /*dimension=*/1);
	EXPECT_EQ(space.n_cliques(), 7 * 3);
}

//-----------------------------------
// The owning dimension
//-----------------------------------

TEST(CliqueSpace, puts_a_zero_in_the_owning_dimension_of_every_index) {
	const auto species = species_tree_space();
	for (size_t clique = 0; clique < species.n_cliques(); ++clique) {
		EXPECT_EQ(species.index_of(clique)[SPECIES], 0u) << "clique " << clique;
	}

	const auto molecules = molecule_tree_space();
	for (size_t clique = 0; clique < molecules.n_cliques(); ++clique) {
		EXPECT_EQ(molecules.index_of(clique)[MOLECULES], 0u) << "clique " << clique;
	}
}

TEST(CliqueSpace, reads_the_owning_dimension_of_a_cell_as_nothing) {
	// Every node of a clique sits in the same clique. So whatever the owning dimension carries --
	// a leaf, an internal node, a root -- the answer does not move. This is what lets a caller hand
	// over the whole cell instead of naming the slot to blank.
	const auto species = species_tree_space();
	for (size_t molecule_leaf = 0; molecule_leaf < N_MOLECULE_LEAVES; ++molecule_leaf) {
		const size_t expected = species.clique_of({0, molecule_leaf});
		for (size_t node = 0; node < 13; ++node) {
			EXPECT_EQ(species.clique_of({node, molecule_leaf}), expected)
			    << "node " << node << " of the clique at molecule leaf " << molecule_leaf;
		}
	}
}

//-----------------------------------
// The round trip
//-----------------------------------

TEST(CliqueSpace, walks_from_a_clique_to_its_cell_and_back) {
	for (const auto &space : {species_tree_space(), molecule_tree_space()}) {
		for (size_t clique = 0; clique < space.n_cliques(); ++clique) {
			EXPECT_EQ(space.clique_of(space.index_of(clique)), clique)
			    << "dimension " << space.dimension() << ", clique " << clique;
		}
	}
}

TEST(CliqueSpace, walks_from_a_clique_to_its_cell_and_back_over_three_dimensions) {
	// The case the two old conventions disagreed on. With one dimension above 1 a row-major
	// subscript and a column-major stride land on the same number; with two, they do not.
	for (size_t dimension = 0; dimension < 3; ++dimension) {
		const TCliqueSpace<3> space({7, 5, 3}, dimension);
		for (size_t clique = 0; clique < space.n_cliques(); ++clique) {
			EXPECT_EQ(space.clique_of(space.index_of(clique)), clique)
			    << "dimension " << dimension << ", clique " << clique;
		}
	}
}

TEST(CliqueSpace, gives_every_clique_its_own_number) {
	// The round trip alone does not say the map is one to one over cells. This does: every leaf of
	// every other tree names a different clique.
	const TCliqueSpace<3> space({7, 5, 3}, /*dimension=*/1);
	std::array<bool, 7 * 3> seen{};
	for (size_t species_leaf = 0; species_leaf < 7; ++species_leaf) {
		for (size_t tissue_leaf = 0; tissue_leaf < 3; ++tissue_leaf) {
			const size_t clique = space.clique_of({species_leaf, 0, tissue_leaf});
			ASSERT_LT(clique, seen.size());
			EXPECT_FALSE(seen[clique]) << "two cells reached clique " << clique;
			seen[clique] = true;
		}
	}
	for (size_t clique = 0; clique < seen.size(); ++clique) {
		EXPECT_TRUE(seen[clique]) << "no cell reached clique " << clique;
	}
}

//-----------------------------------
// What it refuses
//-----------------------------------

TEST(CliqueSpace, refuses_a_dimension_the_clique_space_does_not_have) {
	EXPECT_THROW((TCliqueSpace<>({N_SPECIES_LEAVES, N_MOLECULE_LEAVES}, 2)), std::invalid_argument);
	EXPECT_THROW((TCliqueSpace<3>({7, 5, 3}, 3)), std::invalid_argument);
}

} // namespace
