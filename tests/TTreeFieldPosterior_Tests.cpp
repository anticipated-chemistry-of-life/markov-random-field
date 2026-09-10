//
// The posterior of a tree field, and the file it is written to.
//
// A tree field is the leaf block of a node state, and a node state carries no counter of its own,
// so the counter here has to reach into the node state through the *field's* index. That is the
// property the suite is mostly about: it counts the leaf block and nothing above it, whichever
// storage backs either container.
//
// The bodies are written against the storage concepts and instantiated for every backend pairing,
// the way the node-state file suite is.
//

#include "backend_pairings.h"
#include "constants.h"
#include "field/tree_field_posterior.h"
#include "storages/storage_concepts.h"
#include "temp_file.h"
#include "tree/TPhylogeny.h"
#include "tree/io/node_state_columns.h"
#include "tree/io/write_tree_field.h"
#include "gtest/gtest.h"

#include <cstddef>
#include <sstream>
#include <string>
#include <vector>

namespace {

using backends::AllBackends;
using backends::make_storage;
using phylo::edge;

constexpr size_t SPECIES   = 0;
constexpr size_t MOLECULES = 1;

/// Two leaves per tree, and an internal node above each pair, so the species node state has rows
/// the tree field must not be counted from.
struct TTwoTrees {
	TPhylogeny species  = build_phylogeny({edge("frog", "s_root"), edge("newt", "s_root")});
	TPhylogeny molecule = build_phylogeny(
	    {edge("acid", "m_inner"), edge("base", "m_inner"), edge("m_inner", "m_root")});

	[[nodiscard]] std::vector<TNodeStateColumn> columns() const {
		return {TNodeStateColumn{&species, "species"}, TNodeStateColumn{&molecule, "molecules"}};
	}
	[[nodiscard]] IndexArray field_shape() const {
		return IndexArray{species.n_leaves(), molecule.n_leaves()};
	}
	[[nodiscard]] IndexArray species_shape() const {
		return node_state_dimensions(field_shape(), SPECIES, species);
	}
	[[nodiscard]] IndexArray cell(const std::string &species_node,
	                              const std::string &molecule_node) const {
		return IndexArray{species.index_of(species_node), molecule.index_of(molecule_node)};
	}
};

template<typename Backends> class TreeFieldPosterior : public ::testing::Test {};
TYPED_TEST_SUITE(TreeFieldPosterior, AllBackends);

// -------------------------------------------------------------------------
// Counting
// -------------------------------------------------------------------------

TYPED_TEST(TreeFieldPosterior, counts_a_cell_the_tree_field_holds_a_one_at) {
	const TTwoTrees trees;
	auto Y = make_storage<typename TypeParam::field>(trees.field_shape());
	auto Z = make_storage<typename TypeParam::node_state>(trees.species_shape());
	Z.insert_one(Z.get_linear_index_in_container_space(trees.cell("frog", "acid")));

	TTreeFieldPosterior posterior;
	posterior.initialize(Y.total_size_of_container_space(), 1);
	for (size_t iteration = 0; iteration < 4; ++iteration) {
		posterior.add_to_counter(iteration, Y, Z);
	}

	const size_t counted = Y.get_linear_index_in_container_space(trees.cell("frog", "acid"));
	const size_t other   = Y.get_linear_index_in_container_space(trees.cell("newt", "base"));
	EXPECT_EQ(posterior.get_total_counts(), 4u);
	EXPECT_EQ(posterior.get_counter(counted), 4u);
	EXPECT_EQ(posterior.get_counter(other), 0u);
	EXPECT_DOUBLE_EQ(posterior.get_fraction_of_ones(counted), 1.0);
	EXPECT_DOUBLE_EQ(posterior.get_fraction_of_ones(other), 0.0);
}

/// The node state reaches above the leaves, and the tree field is the leaf block alone. A one at
/// an internal node is not a cell of the tree field, and no cell of the posterior counts it.
TYPED_TEST(TreeFieldPosterior, counts_no_cell_above_the_leaves) {
	const TTwoTrees trees;
	auto Y = make_storage<typename TypeParam::field>(trees.field_shape());
	auto Z = make_storage<typename TypeParam::node_state>(trees.species_shape());
	for (const auto &molecule_leaf : {"acid", "base"}) {
		Z.insert_one(Z.get_linear_index_in_container_space(trees.cell("s_root", molecule_leaf)));
	}

	TTreeFieldPosterior posterior;
	posterior.initialize(Y.total_size_of_container_space(), 1);
	posterior.add_to_counter(0, Y, Z);

	for (size_t cell = 0; cell < posterior.size(); ++cell) {
		EXPECT_EQ(posterior.get_counter(cell), 0u) << "cell " << cell;
	}
}

/// The counter is thinned the way the field's is, so both posteriors read the same iterations.
TYPED_TEST(TreeFieldPosterior, counts_one_iteration_in_the_thinning_factor) {
	const TTwoTrees trees;
	auto Y = make_storage<typename TypeParam::field>(trees.field_shape());
	auto Z = make_storage<typename TypeParam::node_state>(trees.species_shape());
	Z.insert_one(Z.get_linear_index_in_container_space(trees.cell("frog", "acid")));

	TTreeFieldPosterior posterior;
	posterior.initialize(Y.total_size_of_container_space(), 3);
	for (size_t iteration = 0; iteration < 9; ++iteration) {
		posterior.add_to_counter(iteration, Y, Z);
	}

	EXPECT_EQ(posterior.get_total_counts(), 3u);
	EXPECT_EQ(
	    posterior.get_counter(Y.get_linear_index_in_container_space(trees.cell("frog", "acid"))),
	    3u);
}

TYPED_TEST(TreeFieldPosterior, reports_no_posterior_before_it_has_counted_anything) {
	const TTwoTrees trees;
	auto Y = make_storage<typename TypeParam::field>(trees.field_shape());

	TTreeFieldPosterior posterior;
	posterior.initialize(Y.total_size_of_container_space(), 1);

	EXPECT_EQ(posterior.get_total_counts(), 0u);
	EXPECT_DOUBLE_EQ(posterior.get_fraction_of_ones(0), 0.0);
}

/// Burn-in clears the counters, numerator and denominator together, so what follows is a
/// posterior of the chain and not of the burn-in.
TYPED_TEST(TreeFieldPosterior, forgets_everything_it_counted_when_the_counts_are_reset) {
	const TTwoTrees trees;
	auto Y = make_storage<typename TypeParam::field>(trees.field_shape());
	auto Z = make_storage<typename TypeParam::node_state>(trees.species_shape());
	Z.insert_one(Z.get_linear_index_in_container_space(trees.cell("frog", "acid")));

	TTreeFieldPosterior posterior;
	posterior.initialize(Y.total_size_of_container_space(), 1);
	posterior.add_to_counter(0, Y, Z);
	posterior.reset_counts();

	EXPECT_EQ(posterior.get_total_counts(), 0u);
	EXPECT_EQ(
	    posterior.get_counter(Y.get_linear_index_in_container_space(trees.cell("frog", "acid"))),
	    0u);
}

/// Every backend pairing, over the shapes the storage suites use: a fraction is a probability.
TYPED_TEST(TreeFieldPosterior, every_fraction_is_a_probability) {
	for (const auto &pair : backends::tree_pairs()) {
		SCOPED_TRACE(pair.name);
		auto Y = make_storage<typename TypeParam::field>(backends::field_shape(pair));
		auto Z = make_storage<typename TypeParam::node_state>(backends::species_shape(pair));
		backends::seed_ones(Z, 7);

		TTreeFieldPosterior posterior;
		posterior.initialize(Y.total_size_of_container_space(), 2);
		for (size_t iteration = 0; iteration < 10; ++iteration) {
			posterior.add_to_counter(iteration, Y, Z);
		}

		for (size_t cell = 0; cell < posterior.size(); ++cell) {
			const double fraction = posterior.get_fraction_of_ones(cell);
			EXPECT_GE(fraction, 0.0) << "cell " << cell;
			EXPECT_LE(fraction, 1.0) << "cell " << cell;
		}
	}
}

// -------------------------------------------------------------------------
// The file
// -------------------------------------------------------------------------

TYPED_TEST(TreeFieldPosterior, writes_a_row_per_counted_leaf_pair_naming_both_leaves) {
	const TTwoTrees trees;
	auto Y = make_storage<typename TypeParam::field>(trees.field_shape());
	auto Z = make_storage<typename TypeParam::node_state>(trees.species_shape());
	Z.insert_one(Z.get_linear_index_in_container_space(trees.cell("newt", "base")));

	TTreeFieldPosterior posterior;
	posterior.initialize(Y.total_size_of_container_space(), 1);
	posterior.add_to_counter(0, Y, Z);
	posterior.add_to_counter(1, Y, Z);

	const TTempFile file("tree_field_posterior.txt");
	write_tree_field_posterior(file.path(), Y, Z, posterior, trees.columns());

	// Only the one counted cell has a row; the header names the two trees between the state and
	// the fraction.
	EXPECT_EQ(file.content(), "position\tZ_state\tspecies\tmolecules\tfraction_of_one\n"
	                          "3\t1\tnewt\tbase\t1\n");
}

/// A cell the chain never held a one at carries no posterior, so it is left out -- the rule the
/// field's own posterior file follows.
TYPED_TEST(TreeFieldPosterior, leaves_out_a_leaf_pair_that_was_never_a_one) {
	const TTwoTrees trees;
	auto Y = make_storage<typename TypeParam::field>(trees.field_shape());
	auto Z = make_storage<typename TypeParam::node_state>(trees.species_shape());

	TTreeFieldPosterior posterior;
	posterior.initialize(Y.total_size_of_container_space(), 1);
	posterior.add_to_counter(0, Y, Z);

	const TTempFile file("empty_tree_field_posterior.txt");
	write_tree_field_posterior(file.path(), Y, Z, posterior, trees.columns());

	EXPECT_EQ(file.content(), "position\tZ_state\tspecies\tmolecules\tfraction_of_one\n");
}

} // namespace
