//
// What a node-state file may say, and what comes back when it says it.
//
// The reader and the writer used to take trees, so neither had a test: building a tree means
// standing up three live stattools parameters and a command line. They take a topology and a name
// per column instead (node_state_columns.h), and a phylogeny is a value a test can build.
//
// The file this suite cares most about is one written *before* ADR-0005 extended the node state
// down to the leaves. Such a file names only internal nodes in its own column and holds no leaf
// row at all. It has to keep loading, because a long run must not be thrown away by a model
// change (#40). Where its leaf rows come from is the chain start, so the last test runs the
// reader and the start together.
//
// The bodies are written against the storage concepts and instantiated for every backend pairing,
// the way the leaf-layer suite is: a node state a file is read into is a storage like any other.
//

#include "backend_pairings.h"
#include "constants.h"
#include "coretools/Main/TError.h"
#include "field/leaf_layer_start.h"
#include "field/link_backend.h"
#include "phylogeny_generators.h"
#include "storages/storage_concepts.h"
#include "temp_file.h"
#include "tree/TPhylogeny.h"
#include "tree/io/node_state_columns.h"
#include "tree/io/read_Z.h"
#include "tree/io/write_Z.h"
#include "gtest/gtest.h"

#include <cstddef>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

using backends::AllBackends;
using backends::make_storage;
using phylo::edge;

/// Which column belongs to which tree. The species tree is dimension 0 and the molecule tree
/// dimension 1, as they are everywhere else.
constexpr size_t SPECIES   = 0;
constexpr size_t MOLECULES = 1;

// -------------------------------------------------------------------------
// The two trees every fixture below is written against
// -------------------------------------------------------------------------

/// A named pair of trees, small enough that a file over them can be written out by hand.
///
/// The species tree is a star: two leaves under one root, and no node between them. The molecule
/// tree has a node of every kind: two leaves, an internal non-root above them, and a root above
/// that. So the molecule column has both an internal node that is not a root and one that is,
/// which is the distinction the old writer made and the reader no longer does.
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
	[[nodiscard]] IndexArray molecule_shape() const {
		return node_state_dimensions(field_shape(), MOLECULES, molecule);
	}

	/// The cell two node names address, in either node state: every coordinate is a node index.
	[[nodiscard]] IndexArray cell(const std::string &species_node,
	                              const std::string &molecule_node) const {
		return IndexArray{species.index_of(species_node), molecule.index_of(molecule_node)};
	}
};

/// A species node state as the old model wrote one: every row names an internal node of the
/// species tree, and no row names a leaf, because the old node state had no leaf rows to write.
///
/// The position column holds what the old writer put there -- an index into the old internal-node
/// space -- and s_root's cell does not sit at either of those indices now.
const std::string OLD_SPECIES_FILE = "species\tmolecules\tposition\tZ_state\n"
                                     "s_root\tacid\t0\t1\n"
                                     "s_root\tbase\t1\t0\n";

/// The same, for the molecule node state: the molecule column names the internal node and the
/// root, and the species column names leaves.
const std::string OLD_MOLECULE_FILE = "species\tmolecules\tposition\tZ_state\n"
                                      "frog\tm_inner\t0\t1\n"
                                      "newt\tm_root\t3\t1\n";

// -------------------------------------------------------------------------
// Helpers
// -------------------------------------------------------------------------

/// Whether a node state holds this cell, read the way its own index reads it.
template<typename NodeState> bool holds(const NodeState &Z, const IndexArray &cell) {
	return Z.is_one(Z.get_linear_index_in_container_space(cell));
}

template<typename NodeState> size_t number_of_ones(const NodeState &Z) {
	size_t ones = 0;
	for (size_t i = 0; i < Z.total_size_of_container_space(); ++i) {
		if (Z.is_one(i)) { ++ones; }
	}
	return ones;
}

/// coretools' user error is a CTAD class template, so it cannot be named in an exception
/// declaration and EXPECT_THROW cannot be used on it. Catching the err::TError base and asserting
/// it is *not* a dev error is what "the user gets a clear error" actually means here.
template<typename F> void expect_user_error(F &&f) {
	try {
		f();
	} catch (coretools::err::TError &e) {
		EXPECT_FALSE(e.isDevError()) << "expected a user error, got a dev error: " << e.what();
		return;
	}
	ADD_FAILURE() << "expected a user error, but nothing was thrown";
}

/// A written file read back as what it says about each cell: the two node names to the state.
/// The position column is deliberately dropped -- nothing reads it, and this suite says so.
std::map<std::pair<std::string, std::string>, bool> states_by_name(const std::string &text) {
	std::map<std::pair<std::string, std::string>, bool> states;
	std::istringstream lines(text);
	std::string line;
	std::getline(lines, line); // the header
	while (std::getline(lines, line)) {
		if (line.empty()) { continue; }
		std::istringstream fields(line);
		std::string species_node;
		std::string molecule_node;
		std::string position;
		std::string state;
		std::getline(fields, species_node, '\t');
		std::getline(fields, molecule_node, '\t');
		std::getline(fields, position, '\t');
		std::getline(fields, state, '\t');
		states[{species_node, molecule_node}] = (state == "1");
	}
	return states;
}

// -------------------------------------------------------------------------
// The suite, over all four storage pairings
// -------------------------------------------------------------------------

template<typename Backends> class NodeStateFile : public ::testing::Test {
public:
	using Field     = typename Backends::field;
	using NodeState = typename Backends::node_state;
};

TYPED_TEST_SUITE(NodeStateFile, AllBackends);

/// A row is resolved through its node names. The index the writer also emits is not read, so a
/// file survives any change to how nodes are numbered.
TYPED_TEST(NodeStateFile, resolves_a_row_by_node_name_and_not_by_the_written_index) {
	using NodeState = typename TestFixture::NodeState;

	const TTwoTrees trees;
	// The position column says 99, which is outside the container space altogether. A reader that
	// trusted it could not have arrived at the right cell.
	const TTempFile file("stale_index_Z.txt", "species\tmolecules\tposition\tZ_state\n"
	                                          "s_root\tbase\t99\t1\n");

	auto Z = make_storage<NodeState>(trees.species_shape());
	read_Z_from_file(file.path(), Z, trees.columns(), SPECIES);

	EXPECT_TRUE(holds(Z, trees.cell("s_root", "base")));
	EXPECT_EQ(number_of_ones(Z), 1u);
}

/// The check that a column naming an internal node really is one is relaxed to "this name exists
/// in the tree this column belongs to", because the node state now contains leaves (ADR-0005).
TYPED_TEST(NodeStateFile, accepts_a_leaf_in_the_column_of_its_own_tree) {
	using NodeState = typename TestFixture::NodeState;

	const TTwoTrees trees;
	const TTempFile file("leaf_row_Z.txt", "species\tmolecules\tposition\tZ_state\n"
	                                       "frog\tacid\t0\t1\n");

	auto Z = make_storage<NodeState>(trees.species_shape());
	read_Z_from_file(file.path(), Z, trees.columns(), SPECIES);

	EXPECT_TRUE(holds(Z, trees.cell("frog", "acid")));
	EXPECT_EQ(number_of_ones(Z), 1u);
}

/// A file written before the model change loads, and says what it always said: the internal rows
/// it holds, and nothing about the leaves.
TYPED_TEST(NodeStateFile, loads_a_file_written_before_the_model_change) {
	using NodeState = typename TestFixture::NodeState;

	const TTwoTrees trees;
	const TTempFile species_file("old_species_Z.txt", OLD_SPECIES_FILE);
	const TTempFile molecule_file("old_molecule_Z.txt", OLD_MOLECULE_FILE);

	auto Z_species = make_storage<NodeState>(trees.species_shape());
	read_Z_from_file(species_file.path(), Z_species, trees.columns(), SPECIES);
	EXPECT_TRUE(holds(Z_species, trees.cell("s_root", "acid")));
	EXPECT_FALSE(holds(Z_species, trees.cell("s_root", "base")));
	EXPECT_EQ(number_of_ones(Z_species), 1u); // no leaf row, so no other cell moved

	auto Z_molecule = make_storage<NodeState>(trees.molecule_shape());
	read_Z_from_file(molecule_file.path(), Z_molecule, trees.columns(), MOLECULES);
	EXPECT_TRUE(holds(Z_molecule, trees.cell("frog", "m_inner")));
	EXPECT_TRUE(holds(Z_molecule, trees.cell("newt", "m_root")));
	EXPECT_EQ(number_of_ones(Z_molecule), 2u);
}

/// The relaxation goes exactly this far: a name still has to be a node of the tree its column
/// belongs to.
TYPED_TEST(NodeStateFile, rejects_a_name_that_is_not_in_that_columns_tree) {
	using NodeState = typename TestFixture::NodeState;

	const TTwoTrees trees;
	// "acid" is a molecule, so it names nothing in the species column.
	const TTempFile file("foreign_name_Z.txt", "species\tmolecules\tposition\tZ_state\n"
	                                           "acid\tacid\t0\t1\n");

	auto Z = make_storage<NodeState>(trees.species_shape());
	expect_user_error([&] { read_Z_from_file(file.path(), Z, trees.columns(), SPECIES); });
}

/// Only the node state's own column reaches past the leaves. Every other column indexes a leaf,
/// so an internal node there is still an error.
TYPED_TEST(NodeStateFile, rejects_an_internal_node_in_a_foreign_column) {
	using NodeState = typename TestFixture::NodeState;

	const TTwoTrees trees;
	const TTempFile file("foreign_internal_Z.txt", "species\tmolecules\tposition\tZ_state\n"
	                                               "frog\tm_inner\t0\t1\n");

	// Read as the species node state, the molecule column is foreign, and m_inner is not a leaf.
	auto Z_species = make_storage<NodeState>(trees.species_shape());
	expect_user_error([&] { read_Z_from_file(file.path(), Z_species, trees.columns(), SPECIES); });

	// The same row read as the molecule node state is fine: m_inner is then in its own column.
	auto Z_molecule = make_storage<NodeState>(trees.molecule_shape());
	read_Z_from_file(file.path(), Z_molecule, trees.columns(), MOLECULES);
	EXPECT_TRUE(holds(Z_molecule, trees.cell("frog", "m_inner")));
}

TYPED_TEST(NodeStateFile, rejects_a_file_with_the_wrong_number_of_columns) {
	using NodeState = typename TestFixture::NodeState;

	const TTwoTrees trees;
	const TTempFile file("three_columns_Z.txt", "species\tmolecules\tZ_state\n"
	                                            "frog\tacid\t1\n");

	auto Z = make_storage<NodeState>(trees.species_shape());
	expect_user_error([&] { read_Z_from_file(file.path(), Z, trees.columns(), SPECIES); });
}

/// A file written now carries the leaf rows, because the node state spans them, and it reads back
/// as the node state it was written from.
TYPED_TEST(NodeStateFile, a_written_file_carries_the_leaf_rows_and_reads_back_as_it_was_written) {
	using NodeState = typename TestFixture::NodeState;

	const TTwoTrees trees;
	auto Z = make_storage<NodeState>(trees.species_shape());
	Z.insert_one(Z.get_linear_index_in_container_space(trees.cell("frog", "acid")));   // a leaf
	Z.insert_one(Z.get_linear_index_in_container_space(trees.cell("s_root", "base"))); // a root

	const TTempFile file("written_species_Z.txt");
	write_Z_to_file(file.path(), Z, trees.columns(), /*write_full_Z =*/true);

	const auto states = states_by_name(file.content());
	EXPECT_EQ(states.size(), Z.total_size_of_container_space());
	EXPECT_TRUE(states.at({"frog", "acid"}));
	EXPECT_TRUE(states.at({"s_root", "base"}));
	EXPECT_FALSE(states.at({"newt", "acid"}));
	EXPECT_FALSE(states.at({"frog", "base"}));

	auto read_back = make_storage<NodeState>(trees.species_shape());
	read_Z_from_file(file.path(), read_back, trees.columns(), SPECIES);
	for (size_t i = 0; i < Z.total_size_of_container_space(); ++i) {
		EXPECT_EQ(read_back.is_one(i), Z.is_one(i)) << "at cell " << i;
	}
}

/// The leaf rows an old file does not have come from the chain start, which runs after the file is
/// read: the field takes the records, and both tree fields take the field. The internal rows the
/// file gave stand, which is the whole of what such a file ever held.
TYPED_TEST(NodeStateFile, the_chain_start_fills_the_leaf_rows_an_old_file_does_not_have) {
	using Field     = typename TestFixture::Field;
	using NodeState = typename TestFixture::NodeState;

	const TTwoTrees trees;
	const TTempFile species_file("old_species_Z.txt", OLD_SPECIES_FILE);
	const TTempFile molecule_file("old_molecule_Z.txt", OLD_MOLECULE_FILE);

	auto Z_species = make_storage<NodeState>(trees.species_shape());
	read_Z_from_file(species_file.path(), Z_species, trees.columns(), SPECIES);
	auto Z_molecule = make_storage<NodeState>(trees.molecule_shape());
	read_Z_from_file(molecule_file.path(), Z_molecule, trees.columns(), MOLECULES);

	// One record, so the start has both a leaf pair to write a one at and one to leave at zero.
	auto records = make_storage<Field>(trees.field_shape());
	records.insert_one(records.get_linear_index_in_container_space(trees.cell("newt", "acid")));

	auto Y = make_storage<Field>(trees.field_shape());
	leaf_layer_start::start_the_field_at(records, Y);
	const auto counters =
	    leaf_layer_start::hold_tree_fields_at_the_field<TLinkPolicy>(Y, Z_species, Z_molecule);
	EXPECT_EQ(counters.total(), Y.total_size_of_container_space());

	for (const auto &species_leaf : {"frog", "newt"}) {
		for (const auto &molecule_leaf : {"acid", "base"}) {
			const IndexArray cell = trees.cell(species_leaf, molecule_leaf);
			SCOPED_TRACE(std::string(species_leaf) + "," + molecule_leaf);
			const bool reported = holds(records, cell);
			EXPECT_EQ(holds(Y, cell), reported);
			EXPECT_EQ(holds(Z_species, cell), reported);
			EXPECT_EQ(holds(Z_molecule, cell), reported);
		}
	}

	// The rows above the leaves are the file's, untouched.
	EXPECT_TRUE(holds(Z_species, trees.cell("s_root", "acid")));
	EXPECT_FALSE(holds(Z_species, trees.cell("s_root", "base")));
	EXPECT_TRUE(holds(Z_molecule, trees.cell("frog", "m_inner")));
	EXPECT_TRUE(holds(Z_molecule, trees.cell("newt", "m_root")));
}

} // namespace
