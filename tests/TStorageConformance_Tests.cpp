//
// One suite, both backends.
//
// Every body below is written once, against the concepts in storage_concepts.h, and instantiated
// for the dense and the sparse implementation of each storage. Two parallel suites would drift:
// a property somebody adds to one is a property the other silently never gets asked. Here a
// property is asked of every implementation by construction, and a body that passes for one and
// fails for the other names which one in the test's own name.
//
// The shapes are not a fixture either. They come from the phylogeny generator the tree tests
// already use -- multi-root forests, a deep chain, a wide star -- because the shape is what
// decides whether an index property is a property or a coincidence. A balanced fixture satisfies
// "linear index round-trips" whatever the layout; a container one cell wide, which is what a tree
// with a single leaf asks for, does not.
//
// What is deliberately *not* asserted equal between the backends is `exists`: dense stores the
// whole container space and reports every cell of a clique as existing, where sparse reports only
// the cells it holds. The sweep reads that answer to choose between flipping a cell in place and
// deferring an insert (TClique::_update_current_state), and both routes end at the same state --
// which is what the equivalence tests below check instead. `empty()` is the same story and has a
// test of its own.
//

#include "constants.h"
#include "phylogeny_generators.h"
#include "storages/storage_concepts.h"
#include "storages/y_storage/TStorageYDense.h"
#include "storages/y_storage/TStorageYMatrix.h"
#include "storages/z_storage/TStorageZDense.h"
#include "storages/z_storage/TStorageZMatrix.h"
#include "tree/TPhylogeny.h"
#include "tree/node_state_shape.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <random>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

// -------------------------------------------------------------------------
// The shapes a pair of trees asks the storages to hold
// -------------------------------------------------------------------------

/// The container spaces two trees imply: the field is indexed in leaf space, one dimension per
/// tree, and each tree's node state is that same space with its own dimension spanning every node
/// of that tree instead (TMarkovField::initialize and TTree::_initialize_Z).
///
/// The node state reaches its leaves, so its own dimension is n_nodes and not n_internal_nodes.
/// That is what puts a leaf pair at the same (row, column) in the field and in either node state,
/// and it is asserted as a property below.
///
/// A storage knows nothing of a tree beyond the size of each dimension, so a shape is all a
/// storage test needs from one -- and taking it from a real phylogeny rather than from a literal
/// is what keeps the degenerate cases in: a chain has one leaf, a star has one internal node.
std::vector<IndexArray> shapes_of(const TPhylogeny &first, const TPhylogeny &second) {
	const IndexArray field{first.n_leaves(), second.n_leaves()};
	// through node_state_dimensions, the same function TTree::_initialize_Z sizes the real node
	// state with -- so changing the rule changes what these tests are asserted over
	return {field, node_state_dimensions(field, 0, first), node_state_dimensions(field, 1, second)};
}

/// Every shape the properties below are asserted over: the pairings of four multi-root forests,
/// a deep chain and a wide star, each pairing taken in both orders because the two dimensions of
/// a container are not interchangeable -- one is walked as a row and the other as a column.
/// The forests every shape is taken from. Built once: a phylogeny is immutable, and the tests
/// below pair each with each, so rebuilding them per pairing would be the same six trees six
/// times over.
const std::vector<TPhylogeny> &generated_phylogenies() {
	static const std::vector<TPhylogeny> trees = [] {
		std::mt19937_64 rng(20260828);
		std::vector<TPhylogeny> forests;
		for (size_t n_roots = 1; n_roots <= 4; ++n_roots) {
			forests.push_back(
			    build_phylogeny(phylo::random_forest(rng, 2 * n_roots + 14, n_roots)));
		}
		forests.push_back(build_phylogeny(phylo::chain(60))); // one leaf, 59 internal nodes
		forests.push_back(build_phylogeny(phylo::star(199))); // 199 leaves, one internal node
		return forests;
	}();
	return trees;
}

std::vector<IndexArray> generated_shapes() {
	const auto &trees = generated_phylogenies();

	std::vector<IndexArray> shapes;
	for (const auto &first : trees) {
		for (const auto &second : trees) {
			const auto three = shapes_of(first, second);
			shapes.insert(shapes.end(), three.begin(), three.end());
		}
	}
	// The pairings repeat shapes -- every tree paired with itself gives the same field shape twice
	// over -- and a shape asserted twice says nothing the first time did not.
	std::sort(shapes.begin(), shapes.end());
	shapes.erase(std::unique(shapes.begin(), shapes.end()), shapes.end());
	return shapes;
}

/// How long a chain a field under test is sized for. Short enough that both implementations thin
/// by one -- the sparse counter has 15 bits and the dense one 16, so they thin differently for a
/// chain that is long enough to need it -- which is what lets the counters be compared exactly.
constexpr size_t N_ITERATIONS = 300;

/// Both storages take their dimensions at construction; a field takes the chain length too,
/// because that is what its counter is sized from. Which of the two a storage is, is a question
/// the concepts answer -- the same distinction the sampler makes.
template<typename Storage> Storage make_storage(const IndexArray &dimensions) {
	if constexpr (FieldStorage<Storage>) {
		return Storage(N_ITERATIONS, dimensions);
	} else {
		return Storage(dimensions);
	}
}

/// How often a cell was counted a one. The two fields spell this differently -- the sparse one
/// hands out the packed entry the counter shares a word with, the dense one holds the counter in
/// an array of its own -- and neither spelling is in the concept, so the equivalence between them
/// runs through this pair of overloads rather than through a member.
uint16_t counter_of(const TStorageYMatrix &field, size_t linear_index) {
	return field[linear_index].get_counter();
}
uint16_t counter_of(const TStorageYDense &field, size_t linear_index) {
	return field.get_counter(linear_index);
}

// -------------------------------------------------------------------------
// A sequence of writes, and what it is supposed to leave behind
// -------------------------------------------------------------------------

/// One write, in the vocabulary of the shared surface: the two ways a cell's state is written --
/// in place, or by an insert that also starts the cell's counter over -- and the compaction the
/// sampler runs between iterations.
struct TWrite {
	enum class Kind : uint8_t { set_state, insert, remove_zeros };
	Kind kind           = Kind::set_state;
	size_t linear_index = 0;
	bool state          = false;
};

/// A script over `total_size` cells. Ones are drawn more often than zeros, so that a small
/// container fills several times over -- a script that leaves every cell zero agrees between the
/// backends for no reason worth having -- while the largest shapes are left sparsely populated,
/// which is the regime the sparse backend exists for.
std::vector<TWrite> random_writes(std::mt19937_64 &rng, size_t total_size, size_t n_writes) {
	std::vector<TWrite> writes;
	writes.reserve(n_writes);

	std::uniform_int_distribution<size_t> cell(0, total_size - 1);
	std::uniform_real_distribution<double> unit(0.0, 1.0);
	for (size_t i = 0; i < n_writes; ++i) {
		if (unit(rng) < 0.05) {
			writes.push_back({TWrite::Kind::remove_zeros, 0, false});
			continue;
		}
		const auto kind = unit(rng) < 0.25 ? TWrite::Kind::insert : TWrite::Kind::set_state;
		writes.push_back({kind, cell(rng), unit(rng) < 0.6});
	}
	return writes;
}

/// What both implementations are supposed to hold: a state per cell and, for a field, how often
/// that cell was counted a one. Not a third implementation -- it is what says *which* backend is
/// wrong when the two disagree, which a pairwise comparison alone cannot.
struct TExpectedCells {
	std::vector<uint8_t> states;
	std::vector<size_t> counts;
	explicit TExpectedCells(size_t total_size) : states(total_size, 0), counts(total_size, 0) {}
};

/// Replays a script into a storage and into the expectation beside it.
///
/// The two places the counter enters are the two the implementations had to be made to agree on:
/// an insert writes the cell anew, counter included, and `remove_zeros` erases a cell that is no
/// longer a one -- sparse by dropping it outright, dense by zeroing the counter it keeps.
template<typename Storage>
void replay(Storage &storage, const std::vector<TWrite> &writes, TExpectedCells &expected) {
	for (const auto &write : writes) {
		switch (write.kind) {
		case TWrite::Kind::set_state:
			storage.set_state(write.linear_index, write.state);
			expected.states[write.linear_index] = static_cast<uint8_t>(write.state);
			break;
		case TWrite::Kind::insert:
			if (write.state) {
				storage.insert_one(write.linear_index);
			} else {
				storage.insert_zero(write.linear_index);
			}
			expected.states[write.linear_index] = static_cast<uint8_t>(write.state);
			expected.counts[write.linear_index] = 0;
			break;
		case TWrite::Kind::remove_zeros:
			storage.remove_zeros();
			for (size_t i = 0; i < expected.states.size(); ++i) {
				if (expected.states[i] == 0) { expected.counts[i] = 0; }
			}
			break;
		}
	}
}

/// One script per iteration, so that both backends can be driven through the *same* chain rather
/// than through two chains drawn from the same generator.
std::vector<std::vector<TWrite>> random_script(std::mt19937_64 &rng, size_t total_size,
                                               size_t n_iterations) {
	std::vector<std::vector<TWrite>> script;
	script.reserve(n_iterations);
	for (size_t iteration = 0; iteration < n_iterations; ++iteration) {
		script.push_back(random_writes(rng, total_size, 3));
	}
	return script;
}

/// A chain, as the sampler runs one: the iteration's writes, then the counter.
///
/// The expectation counts every iteration, which is a rule of its own and not the field's:
/// `add_to_counter` counts one iteration in `get_thinning_factor()`, and repeating that expression
/// here would be checking each implementation against itself. The two rules coincide only while
/// the field does not thin, which is what N_ITERATIONS is chosen for and what the assertion below
/// pins. Thinning itself is asserted where the two factors differ, at the end of the file.
template<typename Field>
void run_chain(Field &field, const std::vector<std::vector<TWrite>> &script,
               TExpectedCells &expected) {
	ASSERT_EQ(field.get_thinning_factor(), 1u);
	for (size_t iteration = 0; iteration < script.size(); ++iteration) {
		replay(field, script[iteration], expected);
		field.add_to_counter(iteration);
		for (size_t i = 0; i < expected.states.size(); ++i) {
			if (expected.states[i] != 0) { ++expected.counts[i]; }
		}
	}
}

/// How many writes a shape is worth: enough to fill a small container several times over, capped
/// so that the largest generated shape does not turn the suite into a benchmark.
size_t n_writes_for(size_t total_size) { return std::min<size_t>(4 * total_size, 1000); }

// -------------------------------------------------------------------------
// One body per storage, instantiated for both backends
// -------------------------------------------------------------------------

class StorageNames {
public:
	template<typename Storage> static std::string GetName(int /*index*/) {
		if constexpr (std::is_same_v<Storage, TStorageZMatrix>) {
			return "sparse_internal_state";
		} else if constexpr (std::is_same_v<Storage, TStorageZDense>) {
			return "dense_internal_state";
		} else if constexpr (std::is_same_v<Storage, TStorageYMatrix>) {
			return "sparse_field";
		} else {
			return "dense_field";
		}
	}
};

/// One tree pairing, checked through one backend: the coordinates of a leaf pair address the same
/// cell in the field and in either node state, with nothing converted in between.
template<typename FieldT, typename NodeStateT>
void check_the_leaf_block_needs_no_conversion(const TPhylogeny &first, const TPhylogeny &second,
                                              size_t &n_pairings_with_rows_the_field_lacks) {
	const auto shapes                 = shapes_of(first, second);
	const IndexArray field_shape      = shapes[0];
	const IndexArray node_state_shape = shapes[1]; // the first tree's

	// Checked against the phylogenies rather than against each other: the owning dimension spans
	// every node of its tree, and the other one stays in leaf space. Comparing the two shapes
	// would only compare two outputs of one function.
	ASSERT_EQ(node_state_shape[0], first.n_nodes());
	ASSERT_EQ(node_state_shape[1], second.n_leaves());
	ASSERT_EQ(field_shape[0], first.n_leaves());
	if (node_state_shape[0] > field_shape[0]) { ++n_pairings_with_rows_the_field_lacks; }

	const FieldT field          = make_storage<FieldT>(field_shape);
	const NodeStateT node_state = make_storage<NodeStateT>(node_state_shape);

	for (size_t leaf = 0; leaf < first.n_leaves(); ++leaf) {
		// ADR-0004: a leaf's index in leaf space *is* its node index. That is the whole reason
		// the row needs no arithmetic -- without it this would be a subtraction that changed.
		ASSERT_TRUE(first.is_leaf(leaf)) << "node " << leaf << " should be a leaf";
		ASSERT_EQ(first.leaf_index(leaf), leaf) << "leaf " << leaf;

		for (size_t column = 0; column < second.n_leaves(); ++column) {
			const IndexArray cell{leaf, column};
			// the same coordinates are addressable in both, and mean the same leaf pair in both
			const size_t in_field      = field.get_linear_index_in_container_space(cell);
			const size_t in_node_state = node_state.get_linear_index_in_container_space(cell);
			ASSERT_EQ(field.get_multi_dimensional_index(in_field), cell);
			ASSERT_EQ(node_state.get_multi_dimensional_index(in_node_state), cell);
			// Stronger, and specific to the tree that owns the row dimension: only its *row* count
			// grew, so the field and its node state have the same column count and a leaf pair
			// lands on the same linear index in both. Nothing has to be translated, not even the
			// flat offset.
			ASSERT_EQ(in_node_state, in_field) << "leaf " << leaf << ", column " << column;
		}
	}
}

template<typename Storage> class StorageConformance : public ::testing::Test {};
using Storages = ::testing::Types<TStorageZMatrix, TStorageZDense, TStorageYMatrix, TStorageYDense>;
TYPED_TEST_SUITE(StorageConformance, Storages, StorageNames);

// -------------------------------------------------------------------------
// The leaf block of a node state is the field, index for index
// -------------------------------------------------------------------------

TEST(NodeStateLeafBlock,
     a_leaf_pair_sits_at_the_same_row_and_column_in_the_field_and_the_node_state) {
	// ADR-0005 extends each tree's node state down to its leaves and leaves the other dimension in
	// leaf space, so the field and either node state address one leaf pair identically and the
	// subtract-the-leaf-count conversion disappears rather than getting more complex.
	//
	// Over generated forests rather than a fixture, because the shape is what makes this a property
	// instead of a coincidence: a star has one internal node, so its node state is barely taller
	// than its field, and a chain has one leaf, so its field is one cell wide.
	size_t n_pairings_with_rows_the_field_lacks = 0;

	for (const auto &first : generated_phylogenies()) {
		for (const auto &second : generated_phylogenies()) {
			ASSERT_NO_FATAL_FAILURE(
			    (check_the_leaf_block_needs_no_conversion<TStorageYMatrix, TStorageZMatrix>(
			        first, second, n_pairings_with_rows_the_field_lacks)));
			ASSERT_NO_FATAL_FAILURE(
			    (check_the_leaf_block_needs_no_conversion<TStorageYDense, TStorageZDense>(
			        first, second, n_pairings_with_rows_the_field_lacks)));
		}
	}

	// Without this the body above would pass vacuously on a node state that never had rows the
	// field lacks. Since the shapes come from node_state_dimensions, this guards production's rule
	// and not a copy of it: revert that rule and this test goes red.
	EXPECT_GT(n_pairings_with_rows_the_field_lacks, 0u);
}

TYPED_TEST(StorageConformance, index_conversion_round_trips_over_every_generated_shape) {
	for (const auto &shape : generated_shapes()) {
		const auto storage = make_storage<TypeParam>(shape);
		ASSERT_EQ(storage.total_size_of_container_space(), shape[0] * shape[1]);

		for (size_t linear = 0; linear < storage.total_size_of_container_space(); ++linear) {
			const auto multidim = storage.get_multi_dimensional_index(linear);
			ASSERT_LT(multidim[0], shape[0]) << "shape " << shape[0] << "x" << shape[1];
			ASSERT_LT(multidim[1], shape[1]) << "shape " << shape[0] << "x" << shape[1];
			ASSERT_EQ(storage.get_linear_index_in_container_space(multidim), linear)
			    << "shape " << shape[0] << "x" << shape[1];
			// Row-major: the linear index is the cell's position in a row-by-row walk, which is
			// what lets a caller hold one number where the storage holds two.
			ASSERT_EQ(multidim[0] * shape[1] + multidim[1], linear)
			    << "shape " << shape[0] << "x" << shape[1];
		}
	}
}

TYPED_TEST(StorageConformance, every_cell_holds_the_state_it_was_last_written) {
	std::mt19937_64 rng(20260828);
	for (const auto &shape : generated_shapes()) {
		auto storage         = make_storage<TypeParam>(shape);
		const size_t n_cells = storage.total_size_of_container_space();
		TExpectedCells expected(n_cells);

		replay(storage, random_writes(rng, n_cells, n_writes_for(n_cells)), expected);
		for (size_t i = 0; i < n_cells; ++i) {
			ASSERT_EQ(storage.is_one(i), expected.states[i] != 0)
			    << "cell " << i << " of shape " << shape[0] << "x" << shape[1];
		}
	}
}

// The two ways a clique runs through a container: along the last dimension, where the cells are
// consecutive, and along the first, where they are a whole row apart. Everything the sweep reads
// comes through this call, so it has to answer what the point lookups answer -- for the sparse
// implementations that means the line walk has to find every stored cell of the clique.
TYPED_TEST(StorageConformance, fill_current_state_agrees_with_the_point_lookups) {
	std::mt19937_64 rng(20260828);
	std::vector<uint8_t> state;
	std::vector<uint8_t> exists;
	std::vector<size_t> linear;

	for (const auto &shape : generated_shapes()) {
		auto storage         = make_storage<TypeParam>(shape);
		const size_t n_cells = storage.total_size_of_container_space();
		TExpectedCells expected(n_cells);
		replay(storage, random_writes(rng, n_cells, n_writes_for(n_cells)), expected);

		auto check_clique = [&](const IndexArray &start, size_t K, size_t increment) {
			storage.fill_current_state(start, K, increment, state, exists, linear);
			const size_t start_linear = start[0] * shape[1] + start[1];
			ASSERT_EQ(state.size(), K);
			ASSERT_EQ(linear.size(), K);
			for (size_t k = 0; k < K; ++k) {
				ASSERT_EQ(linear[k], start_linear + k * increment)
				    << "clique cell " << k << " of shape " << shape[0] << "x" << shape[1];
				ASSERT_EQ(state[k] != 0, storage.is_one(linear[k]))
				    << "clique cell " << k << " of shape " << shape[0] << "x" << shape[1];
			}
		};

		// Along the last dimension: one clique per row, and one starting halfway along it, because
		// a clique need not start at the beginning of its line.
		for (size_t row = 0; row < shape[0]; ++row) {
			check_clique(IndexArray{row, 0}, shape[1], 1);
			const size_t half = shape[1] / 2;
			check_clique(IndexArray{row, half}, shape[1] - half, 1);
			if (::testing::Test::HasFailure()) { FAIL() << "row " << row; }
		}
		// Along the first dimension: one clique per column, the increment being the width of a row.
		for (size_t col = 0; col < shape[1]; ++col) {
			check_clique(IndexArray{0, col}, shape[0], shape[1]);
			const size_t half = shape[0] / 2;
			check_clique(IndexArray{half, col}, shape[0] - half, shape[1]);
			if (::testing::Test::HasFailure()) { FAIL() << "column " << col; }
		}
	}
}

// A clique along the first dimension steps by the width of a row, so in a container one cell wide
// it steps by one -- the same increment a clique along the *last* dimension has. The two are told
// apart by the shape and not by the increment: where there is a single column, a clique of more
// than one cell can only be a column. A container gets a single column whenever the other tree
// has a single leaf, which is exactly what the generator's deep chain is.
TYPED_TEST(StorageConformance, a_clique_runs_down_a_container_that_is_one_cell_wide) {
	auto storage = make_storage<TypeParam>(IndexArray{4, 1});
	storage.insert_one(1);
	storage.insert_one(3);

	std::vector<uint8_t> state;
	std::vector<uint8_t> exists;
	std::vector<size_t> linear;
	storage.fill_current_state(IndexArray{0, 0}, /*K=*/4, /*increment=*/1, state, exists, linear);

	EXPECT_EQ(linear, (std::vector<size_t>{0, 1, 2, 3}));
	EXPECT_EQ(state, (std::vector<uint8_t>{0, 1, 0, 1}));
}

TYPED_TEST(StorageConformance, an_insert_outside_the_container_space_throws) {
	for (const auto &shape : generated_shapes()) {
		auto storage         = make_storage<TypeParam>(shape);
		const size_t n_cells = storage.total_size_of_container_space();
		// Every tree has at least one leaf and at least one internal node, so no generated shape
		// is empty and the last cell below is a cell.
		ASSERT_GT(n_cells, 0u);

		EXPECT_ANY_THROW(storage.insert_one(n_cells));
		EXPECT_ANY_THROW(storage.insert_zero(n_cells));
		EXPECT_NO_THROW(storage.insert_one(n_cells - 1));
		if (::testing::Test::HasFailure()) { FAIL() << "shape " << shape[0] << "x" << shape[1]; }
	}
}

// Compaction is a memory decision, not a model one: sparse reclaims the cells that are no longer
// ones and dense keeps them, and no cell may change state either way.
TYPED_TEST(StorageConformance, remove_zeros_changes_no_cell_state) {
	std::mt19937_64 rng(20260828);
	for (const auto &shape : generated_shapes()) {
		auto storage         = make_storage<TypeParam>(shape);
		const size_t n_cells = storage.total_size_of_container_space();
		TExpectedCells expected(n_cells);
		replay(storage, random_writes(rng, n_cells, n_writes_for(n_cells)), expected);

		std::vector<uint8_t> before(n_cells, 0);
		for (size_t i = 0; i < n_cells; ++i) {
			before[i] = static_cast<uint8_t>(storage.is_one(i));
		}
		storage.remove_zeros();
		for (size_t i = 0; i < n_cells; ++i) {
			ASSERT_EQ(storage.is_one(i), before[i] != 0)
			    << "cell " << i << " of shape " << shape[0] << "x" << shape[1];
		}
	}
}

// -------------------------------------------------------------------------
// The same, for what the field adds: a posterior counter
// -------------------------------------------------------------------------

class FieldNames {
public:
	template<typename Field> static std::string GetName(int /*index*/) {
		if constexpr (std::is_same_v<Field, TStorageYMatrix>) {
			return "sparse";
		} else {
			return "dense";
		}
	}
};

template<typename Field> class FieldConformance : public ::testing::Test {};
using Fields = ::testing::Types<TStorageYMatrix, TStorageYDense>;
TYPED_TEST_SUITE(FieldConformance, Fields, FieldNames);

TYPED_TEST(FieldConformance, the_counter_counts_the_iterations_a_cell_was_a_one) {
	std::mt19937_64 rng(20260828);
	for (const auto &shape : generated_shapes()) {
		auto field           = make_storage<TypeParam>(shape);
		const size_t n_cells = field.total_size_of_container_space();
		TExpectedCells expected(n_cells);

		run_chain(field, random_script(rng, n_cells, N_ITERATIONS), expected);
		for (size_t i = 0; i < n_cells; ++i) {
			ASSERT_EQ(counter_of(field, i), expected.counts[i])
			    << "cell " << i << " of shape " << shape[0] << "x" << shape[1];
			// The fraction is that count over the number of iterations the chain counted, which is
			// what turns a thinned counter back into a posterior probability.
			ASSERT_DOUBLE_EQ(field.get_fraction_of_ones(i),
			                 static_cast<double>(expected.counts[i]) /
			                     static_cast<double>(field.get_total_counts()))
			    << "cell " << i << " of shape " << shape[0] << "x" << shape[1];
		}
	}
}

// One shape is enough here where the tests above sweep them all: resetting is a pass over the
// whole container, and there is no index arithmetic in it for a shape to catch out.
TYPED_TEST(FieldConformance, reset_counts_clears_every_counter_and_keeps_every_state) {
	std::mt19937_64 rng(20260828);
	auto field           = make_storage<TypeParam>(IndexArray{4, 5});
	const size_t n_cells = field.total_size_of_container_space();
	TExpectedCells expected(n_cells);
	run_chain(field, random_script(rng, n_cells, N_ITERATIONS), expected);

	field.reset_counts();
	for (size_t i = 0; i < n_cells; ++i) {
		EXPECT_EQ(counter_of(field, i), 0) << "cell " << i;
		EXPECT_DOUBLE_EQ(field.get_fraction_of_ones(i), 0.0) << "cell " << i;
		EXPECT_EQ(field.is_one(i), expected.states[i] != 0) << "cell " << i;
	}
}

// -------------------------------------------------------------------------
// The two backends together
// -------------------------------------------------------------------------

/// Cell for cell, and clique for clique: what one implementation answers, the other answers.
template<typename First, typename Second>
void expect_same_cells(const First &first, const Second &second, const IndexArray &shape) {
	ASSERT_EQ(first.total_size_of_container_space(), second.total_size_of_container_space());
	const size_t n_cells = first.total_size_of_container_space();

	for (size_t i = 0; i < n_cells; ++i) {
		ASSERT_EQ(first.is_one(i), second.is_one(i))
		    << "cell " << i << " of shape " << shape[0] << "x" << shape[1];
		ASSERT_EQ(first.get_multi_dimensional_index(i), second.get_multi_dimensional_index(i))
		    << "cell " << i << " of shape " << shape[0] << "x" << shape[1];
	}

	std::vector<uint8_t> first_state;
	std::vector<uint8_t> second_state;
	std::vector<uint8_t> ignored_exists;
	std::vector<size_t> first_linear;
	std::vector<size_t> second_linear;
	auto same_clique = [&](const IndexArray &start, size_t K, size_t increment) {
		first.fill_current_state(start, K, increment, first_state, ignored_exists, first_linear);
		second.fill_current_state(start, K, increment, second_state, ignored_exists, second_linear);
		ASSERT_EQ(first_linear, second_linear) << "clique at " << start[0] << "," << start[1]
		                                       << " of shape " << shape[0] << "x" << shape[1];
		ASSERT_EQ(first_state, second_state) << "clique at " << start[0] << "," << start[1]
		                                     << " of shape " << shape[0] << "x" << shape[1];
	};

	for (size_t row = 0; row < shape[0]; ++row) { same_clique(IndexArray{row, 0}, shape[1], 1); }
	for (size_t col = 0; col < shape[1]; ++col) {
		same_clique(IndexArray{0, col}, shape[0], shape[1]);
	}
}

TEST(StorageEquivalence, the_backends_agree_cell_for_cell_after_the_same_writes) {
	std::mt19937_64 rng(20260828);
	for (const auto &shape : generated_shapes()) {
		const size_t n_cells = shape[0] * shape[1];
		const auto writes    = random_writes(rng, n_cells, n_writes_for(n_cells));

		TStorageZMatrix sparse_Z(shape);
		TStorageZDense dense_Z(shape);
		TExpectedCells sparse_expected(n_cells);
		TExpectedCells dense_expected(n_cells);
		replay(sparse_Z, writes, sparse_expected);
		replay(dense_Z, writes, dense_expected);
		expect_same_cells(sparse_Z, dense_Z, shape);

		TStorageYMatrix sparse_Y(N_ITERATIONS, shape);
		TStorageYDense dense_Y(N_ITERATIONS, shape);
		TExpectedCells sparse_field_expected(n_cells);
		TExpectedCells dense_field_expected(n_cells);
		replay(sparse_Y, writes, sparse_field_expected);
		replay(dense_Y, writes, dense_field_expected);
		expect_same_cells(sparse_Y, dense_Y, shape);

		// All four are also what the script says they should be, so a disagreement names the
		// backend that is wrong rather than only the fact that they differ.
		for (size_t i = 0; i < n_cells; ++i) {
			ASSERT_EQ(sparse_Z.is_one(i), sparse_expected.states[i] != 0) << "cell " << i;
			ASSERT_EQ(dense_Z.is_one(i), dense_expected.states[i] != 0) << "cell " << i;
			ASSERT_EQ(sparse_Y.is_one(i), sparse_field_expected.states[i] != 0) << "cell " << i;
			ASSERT_EQ(dense_Y.is_one(i), dense_field_expected.states[i] != 0) << "cell " << i;
		}
		if (::testing::Test::HasFailure()) { FAIL() << "shape " << shape[0] << "x" << shape[1]; }
	}
}

TEST(StorageEquivalence, the_backends_agree_on_the_counter_and_the_fraction_of_ones) {
	std::mt19937_64 rng(20260828);
	for (const auto &shape : generated_shapes()) {
		const size_t n_cells = shape[0] * shape[1];
		const auto script    = random_script(rng, n_cells, N_ITERATIONS);

		TStorageYMatrix sparse(N_ITERATIONS, shape);
		TStorageYDense dense(N_ITERATIONS, shape);
		// The comparison is exact only while neither field thins, which is what N_ITERATIONS is
		// chosen for; the test below is the one that says what happens when they do.
		ASSERT_EQ(sparse.get_thinning_factor(), 1u);
		ASSERT_EQ(dense.get_thinning_factor(), 1u);
		ASSERT_EQ(sparse.get_total_counts(), dense.get_total_counts());

		TExpectedCells sparse_expected(n_cells);
		TExpectedCells dense_expected(n_cells);
		run_chain(sparse, script, sparse_expected);
		run_chain(dense, script, dense_expected);

		for (size_t i = 0; i < n_cells; ++i) {
			ASSERT_EQ(counter_of(sparse, i), counter_of(dense, i))
			    << "cell " << i << " of shape " << shape[0] << "x" << shape[1];
			ASSERT_DOUBLE_EQ(sparse.get_fraction_of_ones(i), dense.get_fraction_of_ones(i))
			    << "cell " << i << " of shape " << shape[0] << "x" << shape[1];
		}
	}
}

// The counters themselves cannot agree over a long chain: the sparse counter shares its word with
// the state bit and holds 15 bits, the dense one has all 16, so the same chain is thinned twice as
// hard on one side as on the other. What the two still have to agree on is the fraction, which is
// the counter over the number of iterations that were counted -- and that is what the chain is
// read for.
TEST(StorageEquivalence, the_fraction_of_ones_survives_a_chain_the_two_thin_differently) {
	// 65534 == 2 * 32767: the sparse counter needs one iteration in two, the dense one every
	// iteration, and both divide the chain exactly.
	constexpr size_t n_iterations = 65534;
	TStorageYMatrix sparse(n_iterations, IndexArray{1, 2});
	TStorageYDense dense(n_iterations, IndexArray{1, 2});
	ASSERT_EQ(sparse.get_thinning_factor(), 2u);
	ASSERT_EQ(dense.get_thinning_factor(), 1u);

	sparse.insert_one(0);
	dense.insert_one(0);
	for (size_t iteration = 0; iteration < n_iterations; ++iteration) {
		sparse.add_to_counter(iteration);
		dense.add_to_counter(iteration);
	}

	// A cell that was a one for the whole chain has posterior probability one, whichever counter
	// was used to say so; a cell that never was has zero.
	EXPECT_DOUBLE_EQ(sparse.get_fraction_of_ones(0), 1.0);
	EXPECT_DOUBLE_EQ(dense.get_fraction_of_ones(0), 1.0);
	EXPECT_DOUBLE_EQ(sparse.get_fraction_of_ones(1), 0.0);
	EXPECT_DOUBLE_EQ(dense.get_fraction_of_ones(1), 0.0);
	// The counters are what differ, and by exactly the ratio of the two thinning factors.
	EXPECT_EQ(counter_of(dense, 0), 2 * counter_of(sparse, 0));
}

// The bulk paths, which the storage concept deliberately leaves out (storage_concepts.h) and
// which therefore have nothing but this to hold them together: the deferred insert the field sweep
// commits after its parallel region, and the whole-space dump the per-iteration traces are written
// from. Both are named after the storage they belong to rather than after what they do, so there
// is one block per storage here rather than one templated body.
TEST(StorageEquivalence, the_backends_agree_on_the_bulk_insert_and_the_whole_space_dump) {
	std::mt19937_64 rng(20260828);
	for (const auto &shape : generated_shapes()) {
		const size_t n_cells = shape[0] * shape[1];

		// Batches, not one list: the sweep accumulates one vector per thread and hands them over
		// together, and a cell may be named by more than one batch.
		std::uniform_int_distribution<size_t> cell(0, n_cells - 1);
		std::vector<std::vector<size_t>> batches(3);
		std::vector<uint8_t> expected(n_cells, 0);
		for (auto &batch : batches) {
			for (size_t i = 0; i < 1 + n_cells / 8; ++i) {
				const size_t linear_index = cell(rng);
				batch.push_back(linear_index);
				expected[linear_index] = 1;
			}
		}

		TStorageZMatrix sparse_Z(shape);
		TStorageZDense dense_Z(shape);
		sparse_Z.insert_in_Z(batches);
		dense_Z.insert_in_Z(batches);
		expect_same_cells(sparse_Z, dense_Z, shape);
		const std::vector<size_t> expected_Z(expected.begin(), expected.end());
		ASSERT_EQ(sparse_Z.get_full_Z_binary_vector(), expected_Z);
		ASSERT_EQ(dense_Z.get_full_Z_binary_vector(), expected_Z);

		TStorageYMatrix sparse_Y(N_ITERATIONS, shape);
		TStorageYDense dense_Y(N_ITERATIONS, shape);
		sparse_Y.insert_in_Y(batches);
		dense_Y.insert_in_Y(batches);
		expect_same_cells(sparse_Y, dense_Y, shape);
		ASSERT_EQ(sparse_Y.get_full_Y_binary_vector(), expected);
		ASSERT_EQ(dense_Y.get_full_Y_binary_vector(), expected);
		ASSERT_EQ(sparse_Y.number_of_ones(), dense_Y.number_of_ones());
		ASSERT_EQ(sparse_Y.dimensions(), dense_Y.dimensions());

		if (::testing::Test::HasFailure()) { FAIL() << "shape " << shape[0] << "x" << shape[1]; }
	}
}

// Which cells a field reports as *stored* is the one place where the two differ and it still
// reaches a file: the sparse matrix holds the cells it was given, the dense array holds the whole
// container space, and the posterior field is written by walking that list.
//
// What makes the two files agree anyway is that a cell only earns a line when it is a one now or
// was counted a one at least once (TMarkovField::_write_only_values_in_Y_vector) -- a cell that is
// neither has every column at its default, and is exactly the kind of cell the two backends
// disagree about holding. That filter is what this asserts, over a chain rather than over a few
// writes, because the counter is half of the rule.
TEST(StorageEquivalence, the_stored_cells_that_carry_a_posterior_are_the_same_ones) {
	std::mt19937_64 rng(20260828);
	size_t n_shapes_where_the_backends_hold_different_cells = 0;
	for (const auto &shape : generated_shapes()) {
		const size_t n_cells = shape[0] * shape[1];
		const auto script    = random_script(rng, n_cells, N_ITERATIONS);

		TStorageYMatrix sparse(N_ITERATIONS, shape);
		TStorageYDense dense(N_ITERATIONS, shape);
		TExpectedCells sparse_expected(n_cells);
		TExpectedCells dense_expected(n_cells);
		run_chain(sparse, script, sparse_expected);
		run_chain(dense, script, dense_expected);

		// The two entry types are different -- the sparse field hands out its packed entry, the
		// dense one a state and a full 16-bit count -- but both answer the two questions the
		// filter asks, which is the whole of what the writer needs from them.
		const auto reported = [](const auto &field) {
			std::vector<std::pair<size_t, bool>> lines;
			for (const auto &[linear_index, cell] : field.get_stored_entries()) {
				if (!cell.is_one() && cell.get_counter() == 0) { continue; }
				lines.emplace_back(linear_index, cell.is_one());
			}
			return lines;
		};
		ASSERT_EQ(reported(sparse), reported(dense)) << "shape " << shape[0] << "x" << shape[1];

		// The dense field really does report the whole space, and the sparse one a subset of it:
		// the difference above is a filter doing its job, not two lists that happen to coincide.
		ASSERT_EQ(dense.get_stored_entries().size(), n_cells);
		ASSERT_LE(sparse.get_stored_entries().size(), n_cells);
		if (sparse.get_stored_entries().size() < dense.get_stored_entries().size()) {
			++n_shapes_where_the_backends_hold_different_cells;
		}

		if (::testing::Test::HasFailure()) { FAIL() << "shape " << shape[0] << "x" << shape[1]; }
	}

	// Without this the test would pass just as happily on a script that leaves the two lists
	// identical, which is the one case in which it proves nothing.
	EXPECT_GT(n_shapes_where_the_backends_hold_different_cells, 0u)
	    << "no generated shape left the two backends holding different cells, so the filter was "
	       "never asked to reconcile anything";
}

// The cursor the likelihoods merge-join through has to yield the *same cells in the same order*
// under both backends, and not merely cells that add up to the same answer.
//
// TLotus and the simple error model both split their sum in two: the cells the cursors yield, term
// by term through an accumulator, and every other cell folded into one closed-form product. Where
// that split falls decides the rounding, so a cursor that yielded what a backend happens to store
// -- the cells it was given, for the sparse field; all of them, for the dense one -- would make the
// same chain reach answers a bit apart under the two. A Metropolis ratio turns that into two
// different chains a few iterations later, which is exactly how this was found.
TEST(StorageEquivalence, the_cursor_the_likelihoods_walk_yields_the_same_cells_under_both) {
	std::mt19937_64 rng(20260828);
	size_t n_shapes_with_a_stored_zero = 0;
	for (const auto &shape : generated_shapes()) {
		const size_t n_cells = shape[0] * shape[1];
		const auto writes    = random_writes(rng, n_cells, n_writes_for(n_cells));

		TStorageYMatrix sparse(N_ITERATIONS, shape);
		TStorageYDense dense(N_ITERATIONS, shape);
		TExpectedCells sparse_expected(n_cells);
		TExpectedCells dense_expected(n_cells);
		replay(sparse, writes, sparse_expected);
		replay(dense, writes, dense_expected);

		// What the cursor is supposed to yield, from the expectation rather than from either
		// field: the cells that are one, in ascending linear-index order.
		std::vector<size_t> expected_ones;
		for (size_t i = 0; i < n_cells; ++i) {
			if (sparse_expected.states[i] != 0) { expected_ones.push_back(i); }
		}

		const auto walked = [](const auto &field) {
			std::vector<size_t> indices;
			for (auto cursor = field.ones_cursor(); cursor.valid(); cursor.advance()) {
				indices.push_back(cursor.linear_index());
			}
			return indices;
		};
		ASSERT_EQ(walked(sparse), expected_ones) << "shape " << shape[0] << "x" << shape[1];
		ASSERT_EQ(walked(dense), expected_ones) << "shape " << shape[0] << "x" << shape[1];

		// The sparse field holds cells that are not ones -- otherwise the two would agree here for
		// no reason worth having, since the interesting case is exactly the cell one backend holds
		// and the other does not.
		if (sparse.get_stored_entries().size() > expected_ones.size()) {
			++n_shapes_with_a_stored_zero;
		}
	}
	EXPECT_GT(n_shapes_with_a_stored_zero, 0u)
	    << "no generated shape left the sparse field holding a zero, so the cursors were never "
	       "asked to agree about one";
}

// The one question the two are allowed to answer differently, and the reason they may: a cell
// inserted as a zero is something the sparse matrix stores and the dense array cannot distinguish
// from any other zero. The caller asking -- TMarkovField, checking that a field it was told to
// hold fixed was in fact read in -- means "is any cell a one", which is what dense answers, and
// the compaction the sampler runs anyway brings the two back into step.
TEST(StorageEquivalence, empty_is_the_one_answer_the_backends_may_differ_on) {
	TStorageZMatrix sparse(IndexArray{2, 3});
	TStorageZDense dense(IndexArray{2, 3});
	EXPECT_EQ(sparse.empty(), dense.empty());

	sparse.insert_one(4);
	dense.insert_one(4);
	EXPECT_FALSE(sparse.empty());
	EXPECT_FALSE(dense.empty());

	sparse.set_state(4, false);
	dense.set_state(4, false);
	EXPECT_FALSE(sparse.empty()); // the cell is still stored, holding a zero
	EXPECT_TRUE(dense.empty());   // no cell is a one

	sparse.remove_zeros();
	dense.remove_zeros();
	EXPECT_TRUE(sparse.empty());
	EXPECT_TRUE(dense.empty());
}

} // namespace
