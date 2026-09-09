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
// One list, because every storage answers the same surface: the run of cells an update writes, the
// handle it writes each cell through, and the cursor a merge join walks. Adding a storage means
// adding a type to it. The suites below split that surface by what they ask about, not by which
// storage answers.
//
// What is deliberately *not* asserted equal between the backends is which cells a storage holds.
// Dense holds the whole container space, sparse holds what it was given. It reaches an update
// through one answer only: whether a write lands in place or waits for a bulk insert. What is
// left to check is that both backends end at the same state, which is what the equivalence tests
// below do. `empty()` is the same story and has a test of its own.
//

#include "constants.h"
#include "phylogeny_generators.h"
#include "storages/TDenseStateArray.h"
#include "storages/TSparseBinaryArray.h"
#include "storages/cell_handle.h"
#include "storages/cell_write.h"
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
/// by one, which is what lets the counters be compared exactly. Both hold their counter in the
/// same 15 bits, so a longer chain would thin them equally rather than differently -- which is a
/// property of its own, and asserted at the end of this file.
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

/// How often a cell was counted a one. Both fields hand out the packed cell the counter shares a
/// word with, and the spelling is not in the concept, so the tests read it through here rather
/// than through a member.
template<typename Field> uint16_t counter_of(const Field &field, size_t linear_index) {
	return field[linear_index].get_counter();
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

/// The same script, written the way an update inside a parallel region writes.
///
/// A `set_state` becomes one call to `write_or_defer`: a handle where the storage gives one, and
/// the storage's own write where it does not. Everything else is the direct path, because an
/// insert and a compaction are not what that write is for.
///
/// The commit sits inside the loop. The sampler commits after a pass that visits each cell once,
/// where a script can write one cell twice, and a deferred cell has to be in the storage before
/// the next write finds it.
template<typename Storage>
void replay_through_write_or_defer(Storage &storage, const std::vector<TWrite> &writes) {
	std::vector<size_t> deferred_inserts;
	for (const auto &write : writes) {
		switch (write.kind) {
		case TWrite::Kind::set_state:
			write_or_defer(storage, write.linear_index, write.state, deferred_inserts);
			break;
		case TWrite::Kind::insert:
			if (write.state) {
				storage.insert_one(write.linear_index);
			} else {
				storage.insert_zero(write.linear_index);
			}
			break;
		case TWrite::Kind::remove_zeros: storage.remove_zeros(); break;
		}
		for (const size_t linear_index : deferred_inserts) { storage.insert_one(linear_index); }
		deferred_inserts.clear();
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
// The runs of cells a shape offers
// -------------------------------------------------------------------------
//
// A run is what an update walks: a start, a count and a stride. A clique is a run down one
// dimension of a node state, and a species leaf's share of the field update is a run along the
// other. Nothing addresses a run as an object any more -- the caller does the arithmetic and asks
// the storage cell by cell -- so a run here is a test fixture and no longer a type under test.

/// One run: where it starts, which dimension it varies, how many cells it holds, and how far
/// apart they are.
///
/// The dimension is carried rather than read off the stride, because the stride cannot say. A run
/// along the first dimension steps by the width of a row, which is one cell in a container one
/// column wide -- the stride a run along the last dimension has. A caller always knows which
/// dimension it is walking; the shape is what turns that into a stride.
struct TRunShape {
	IndexArray start{0, 0};
	size_t varying_dimension = 1;
	size_t n_cells           = 0;
	size_t stride            = 1;
};

/// Every line of a container, as a run: one per row and one per column, each taken whole and
/// again from halfway along, because a run need not start at the beginning of its line.
std::vector<TRunShape> every_run_over(const IndexArray &shape) {
	std::vector<TRunShape> runs;
	for (size_t row = 0; row < shape[0]; ++row) {
		const size_t half = shape[1] / 2;
		runs.push_back({IndexArray{row, 0}, 1, shape[1], 1});
		runs.push_back({IndexArray{row, half}, 1, shape[1] - half, 1});
	}
	for (size_t col = 0; col < shape[1]; ++col) {
		const size_t half = shape[0] / 2;
		runs.push_back({IndexArray{0, col}, 0, shape[0], shape[1]});
		runs.push_back({IndexArray{half, col}, 0, shape[0] - half, shape[1]});
	}
	return runs;
}

/// Four of the runs above: a whole row, a row from halfway along, a whole column, and a column
/// from halfway down. The tests that *write* take these rather than every line, because writing
/// every line of the largest generated shape turns the suite into a benchmark.
std::vector<TRunShape> a_few_runs_over(const IndexArray &shape) {
	const size_t middle_row = shape[0] / 2;
	const size_t middle_col = shape[1] / 2;
	return {
	    {IndexArray{0, 0}, 1, shape[1], 1},
	    {IndexArray{middle_row, middle_col}, 1, shape[1] - middle_col, 1},
	    {IndexArray{0, 0}, 0, shape[0], shape[1]},
	    {IndexArray{middle_row, middle_col}, 0, shape[0] - middle_row, shape[1]},
	};
}

/// Writes one run of cells the way an update inside a parallel region writes one: one call to
/// `write_or_defer` per cell, and the cells the storage could not take handed back so that the
/// caller can commit them once the region ends.
template<typename Storage>
std::vector<size_t> write_run(Storage &storage, const TRunShape &request,
                              const std::vector<uint8_t> &states) {
	const size_t start = storage.get_linear_index_in_container_space(request.start);
	std::vector<size_t> deferred_inserts;
	for (size_t k = 0; k < request.n_cells; ++k) {
		write_or_defer(storage, start + k * request.stride, states[k] != 0, deferred_inserts);
	}
	return deferred_inserts;
}

/// Whether a storage holds every cell of its container space from the moment it is sized. That is
/// what "dense" means here, and it is the one thing a write to a run behaves differently for: a
/// storage that holds every cell never defers.
template<typename Storage>
constexpr bool holds_every_cell =
    std::is_same_v<Storage, TStorageYDense> || std::is_same_v<Storage, TStorageZDense> ||
    std::is_same_v<Storage, TDenseStateArray>;

/// A state per cell of a run, drawn so that both states are asked for.
std::vector<uint8_t> random_states(std::mt19937_64 &rng, size_t n_cells) {
	std::uniform_real_distribution<double> unit(0.0, 1.0);
	std::vector<uint8_t> states(n_cells, 0);
	for (auto &state : states) { state = static_cast<uint8_t>(unit(rng) < 0.5); }
	return states;
}

// -------------------------------------------------------------------------
// One body per storage, instantiated for both backends
// -------------------------------------------------------------------------

class StorageNames {
public:
	template<typename Storage> static std::string GetName(int /*index*/) {
		if constexpr (std::is_same_v<Storage, TStorageZMatrix>) {
			return "sparse_node_state";
		} else if constexpr (std::is_same_v<Storage, TStorageZDense>) {
			return "dense_node_state";
		} else if constexpr (std::is_same_v<Storage, TStorageYMatrix>) {
			return "sparse_field";
		} else if constexpr (std::is_same_v<Storage, TSparseBinaryArray>) {
			return "sparse_binary_array";
		} else if constexpr (std::is_same_v<Storage, TDenseStateArray>) {
			return "dense_state_array";
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
/// Every storage the sampler holds: the two node states, the two fields, and the two arrays the
/// observed data is held in.
using Storages = ::testing::Types<TStorageZMatrix, TStorageZDense, TStorageYMatrix, TStorageYDense,
                                  TSparseBinaryArray, TDenseStateArray>;
TYPED_TEST_SUITE(StorageConformance, Storages, StorageNames);

template<typename Storage> class HandleConformance : public ::testing::Test {};
/// Every storage points an updater at one of its cells, so this is the list above under another
/// name. It is a suite of its own because what it asks about is the handle rather than the state
/// the handle writes.
TYPED_TEST_SUITE(HandleConformance, Storages, StorageNames);

template<typename Storage> class OnesCursorConformance : public ::testing::Test {};
/// The cursor comes with the array a storage is built on, so every storage offers one. Only the
/// field and the observed data are merge-joined through it; a node state answers the same contract
/// because it is the same code, and asking it costs nothing.
TYPED_TEST_SUITE(OnesCursorConformance, Storages, StorageNames);

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

// -------------------------------------------------------------------------
// The run of cells an update walks
// -------------------------------------------------------------------------
//
// These were the window's conformance tests. A window owned the traversal, so it owned the
// answers; now the caller does the arithmetic and the storage answers one cell at a time, and the
// same properties are asked of the storage directly. Which of the two spellings of a write a
// storage brings -- a handle, or a write it makes itself -- is `write_or_defer`'s business and no
// test names it here.

TYPED_TEST(StorageConformance, a_run_addresses_the_cells_the_multidimensional_index_addresses) {
	// The arithmetic a clique and a field row are both written as: start plus k strides. It is the
	// storage's own index conversion that says which cell that is, and a run has to agree with it
	// -- over every line of every generated shape, both ways round. Along the last dimension the
	// cells are consecutive; along the first they are a whole row apart, and in a container one
	// cell wide those two are the same step, which is why the shapes and not the stride decide.
	std::mt19937_64 rng(20260828);

	for (const auto &shape : generated_shapes()) {
		auto storage         = make_storage<TypeParam>(shape);
		const size_t n_cells = storage.total_size_of_container_space();
		TExpectedCells expected(n_cells);
		replay(storage, random_writes(rng, n_cells, n_writes_for(n_cells)), expected);

		for (const auto &request : every_run_over(shape)) {
			const size_t start = storage.get_linear_index_in_container_space(request.start);
			for (size_t k = 0; k < request.n_cells; ++k) {
				const size_t linear = start + k * request.stride;
				IndexArray cell     = request.start;
				cell[request.varying_dimension] += k;
				ASSERT_EQ(storage.get_multi_dimensional_index(linear), cell) << "run cell " << k;
				ASSERT_EQ(storage.is_one(linear), expected.states[linear] != 0) << "run cell " << k;
			}
			if (::testing::Test::HasFailure()) {
				FAIL() << "run at " << request.start[0] << "," << request.start[1] << " of "
				       << request.n_cells << " cells, stride " << request.stride << ", shape "
				       << shape[0] << "x" << shape[1];
			}
		}
	}
}

TYPED_TEST(StorageConformance, a_write_that_lands_is_read_back_and_a_deferred_one_is_not) {
	// What the window's readback contract became. A write the storage takes is there the moment it
	// is made; a write it could not take is not, and stays absent until the deferred list is
	// committed. A walk that reads what it has just written -- a post-order walk does, at every
	// parent -- has to carry that itself, and TCliqueView is where it does (ADR-0006).
	std::mt19937_64 rng(20260828);
	for (const auto &shape : generated_shapes()) {
		auto storage         = make_storage<TypeParam>(shape);
		const size_t n_cells = storage.total_size_of_container_space();
		TExpectedCells expected(n_cells);
		// Written first, so that a run covers cells the sparse storage holds and cells it does
		// not. The second kind is what a write can be deferred for.
		replay(storage, random_writes(rng, n_cells, n_writes_for(n_cells)), expected);

		for (const auto &request : a_few_runs_over(shape)) {
			const auto written  = random_states(rng, request.n_cells);
			const size_t start  = storage.get_linear_index_in_container_space(request.start);
			const auto deferred = write_run(storage, request, written);

			for (size_t k = 0; k < request.n_cells; ++k) {
				const size_t linear = start + k * request.stride;
				const bool waiting =
				    std::find(deferred.begin(), deferred.end(), linear) != deferred.end();
				ASSERT_EQ(storage.is_one(linear), written[k] != 0 && !waiting) << "run cell " << k;
			}
			for (const size_t linear : deferred) { storage.insert_one(linear); }
			for (size_t k = 0; k < request.n_cells; ++k) {
				ASSERT_EQ(storage.is_one(start + k * request.stride), written[k] != 0)
				    << "run cell " << k;
			}
			if (::testing::Test::HasFailure()) {
				FAIL() << "run at " << request.start[0] << "," << request.start[1] << ", shape "
				       << shape[0] << "x" << shape[1];
			}
		}
	}
}

TYPED_TEST(StorageConformance, a_run_along_a_row_writes_a_held_cell_at_once_and_defers_the_rest) {
	// The one place the backends are allowed to differ, and the reason they may: a dense storage
	// holds every cell of its container space, so its write is already in the storage. A sparse
	// storage cannot insert a cell it does not hold without restructuring the container, so that
	// write waits for the bulk insert.
	auto storage = make_storage<TypeParam>(IndexArray{3, 4});
	storage.insert_zero(4); // held by both, in state 0

	const std::vector<uint8_t> written{1, 1, 0, 0};
	const auto deferred = write_run(storage, {IndexArray{1, 0}, 1, 4, 1}, written);

	EXPECT_TRUE(storage.is_one(4)) << "a write to a held cell goes in place";
	if constexpr (holds_every_cell<TypeParam>) {
		EXPECT_TRUE(deferred.empty()) << "a dense storage writes every cell in place";
		EXPECT_TRUE(storage.is_one(5));
	} else {
		EXPECT_EQ(deferred, (std::vector<size_t>{5}))
		    << "a sparse storage defers the cell it does not hold, and only that one";
		EXPECT_FALSE(storage.is_one(5)) << "the write waits rather than restructuring the storage";
	}
	EXPECT_FALSE(storage.is_one(6)) << "an absent cell written to zero already reads as zero";

	for (const size_t linear : deferred) { storage.insert_one(linear); }
	EXPECT_TRUE(storage.is_one(4));
	EXPECT_TRUE(storage.is_one(5));
	EXPECT_FALSE(storage.is_one(6));
}

TYPED_TEST(StorageConformance, a_run_down_a_column_writes_a_held_cell_at_once_and_defers_the_rest) {
	// The same split, down a column. Every storage is keyed by the linear index, so a stride is
	// arithmetic the caller does and the storage never sees; a column run that named the wrong
	// cell would defer a write the storage could have taken, and the commit writes a whole new
	// cell where an in-place write keeps what the cell already carries. That is the
	// clique-to-cell mapping the parity gate's non-square shapes catch.
	auto storage = make_storage<TypeParam>(IndexArray{3, 4});
	storage.insert_zero(4); // (row 1, column 0), held by both, in state 0

	const std::vector<uint8_t> written{0, 1, 1};
	const auto deferred = write_run(storage, {IndexArray{0, 0}, 0, 3, 4}, written); // column 0

	EXPECT_TRUE(storage.is_one(4)) << "a write to a held cell goes in place";
	if constexpr (holds_every_cell<TypeParam>) {
		EXPECT_TRUE(deferred.empty()) << "a dense storage writes every cell in place";
		EXPECT_TRUE(storage.is_one(8));
	} else {
		EXPECT_EQ(deferred, (std::vector<size_t>{8}))
		    << "a sparse storage defers the cell it does not hold, and only that one";
		EXPECT_FALSE(storage.is_one(8)) << "the write waits rather than restructuring the storage";
	}

	for (const size_t linear : deferred) { storage.insert_one(linear); }
	EXPECT_TRUE(storage.is_one(4));
	EXPECT_TRUE(storage.is_one(8));
}

TYPED_TEST(StorageConformance, a_run_down_a_container_one_cell_wide_steps_by_one) {
	// A run along the first dimension steps by the width of a row, so in a container one cell wide
	// it steps by one -- the stride a run along the last dimension has. A chain gives its
	// container a single leaf, and so a single column.
	auto storage = make_storage<TypeParam>(IndexArray{4, 1});
	storage.insert_one(1);
	storage.insert_one(3);

	for (size_t k = 0; k < 4; ++k) {
		EXPECT_EQ(storage.is_one(k), k == 1 || k == 3) << "cell " << k;
	}

	const std::vector<uint8_t> written{1, 1, 0, 0};
	for (const size_t linear : write_run(storage, {IndexArray{0, 0}, 0, 4, 1}, written)) {
		storage.insert_one(linear);
	}
	EXPECT_TRUE(storage.is_one(0));
	EXPECT_TRUE(storage.is_one(1));
	EXPECT_FALSE(storage.is_one(2));
	EXPECT_FALSE(storage.is_one(3));
}

TYPED_TEST(StorageConformance, a_run_over_no_cells_writes_nothing_and_defers_nothing) {
	// A tree with one node gives a clique of one cell, and a run of none is one step further.
	auto storage = make_storage<TypeParam>(IndexArray{3, 4});
	EXPECT_TRUE(write_run(storage, {IndexArray{1, 2}, 1, 0, 1}, {}).empty());
	for (size_t i = 0; i < storage.total_size_of_container_space(); ++i) {
		EXPECT_FALSE(storage.is_one(i)) << "cell " << i;
	}
}

TYPED_TEST(StorageConformance, an_untouched_cell_reads_as_zero) {
	// Dense holds the whole container space; sparse holds what it was given. Which of the two a
	// cell falls under is no longer a question the interface answers, and this is why it need not:
	// an untouched cell reads as zero either way.
	auto storage = make_storage<TypeParam>(IndexArray{3, 4});
	EXPECT_FALSE(storage.is_one(5));
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
// The handle a storage hands an updater
// -------------------------------------------------------------------------

TYPED_TEST(HandleConformance, a_handle_says_of_every_cell_what_is_one_says) {
	// `locate` is `is_one` for a caller that is about to write, so the two have to agree cell for
	// cell. Over a storage that has been written, because an untouched one agrees for the one
	// reason worth nothing: every cell is zero.
	std::mt19937_64 rng(20260828);
	for (const auto &shape : generated_shapes()) {
		auto storage         = make_storage<TypeParam>(shape);
		const size_t n_cells = storage.total_size_of_container_space();
		TExpectedCells expected(n_cells);
		replay(storage, random_writes(rng, n_cells, n_writes_for(n_cells)), expected);

		for (size_t i = 0; i < n_cells; ++i) {
			const auto by_linear_index = storage.locate(i);
			ASSERT_EQ(by_linear_index.is_one, storage.is_one(i))
			    << "cell " << i << " of shape " << shape[0] << "x" << shape[1];
			ASSERT_EQ(by_linear_index.linear_index, i);
			// The one invariant that binds the two halves of the handle: a caller tests the flag
			// and dereferences the pointer, so a flag that is true has to come with a cell.
			ASSERT_EQ(by_linear_index.in_container, by_linear_index.cell != nullptr)
			    << "cell " << i;

			// The same cell, addressed the other way. A caller that already holds a
			// multidimensional index should not have to convert it first.
			const auto by_multidim = storage.locate(storage.get_multi_dimensional_index(i));
			ASSERT_EQ(by_multidim.is_one, by_linear_index.is_one) << "cell " << i;
			ASSERT_EQ(by_multidim.in_container, by_linear_index.in_container) << "cell " << i;
			ASSERT_EQ(by_multidim.linear_index, i);
			ASSERT_EQ(by_multidim.cell, by_linear_index.cell) << "cell " << i;
		}
	}
}

TYPED_TEST(HandleConformance, a_dense_storage_holds_every_cell_and_a_sparse_one_holds_what_it_got) {
	// The one place the handles differ, and the reason they may: a dense storage holds every cell
	// of its container space, so it can always point at one. A sparse storage holds the cells it
	// was given. Both read an untouched cell as zero, which is what keeps the difference off the
	// read-only paths.
	auto storage = make_storage<TypeParam>(IndexArray{3, 4});
	storage.insert_one(2);

	const auto held  = storage.locate(2);
	const auto other = storage.locate(3);
	EXPECT_TRUE(held.in_container);
	EXPECT_NE(held.cell, nullptr);
	EXPECT_TRUE(held.is_one);

	if constexpr (!holds_every_cell<TypeParam>) {
		EXPECT_FALSE(other.in_container) << "the sparse array was never given this cell";
		EXPECT_EQ(other.cell, nullptr);
	} else {
		EXPECT_TRUE(other.in_container) << "a dense storage holds every cell";
		EXPECT_NE(other.cell, nullptr);
	}
	EXPECT_FALSE(other.is_one) << "a cell the storage does not hold reads as zero";
}

TYPED_TEST(HandleConformance, the_helper_writes_a_held_cell_at_once_and_defers_the_rest) {
	// The branch the update loops take. A held cell is written where it lies. A cell the storage
	// does not hold waits, because the insert restructures the container and the loop runs in
	// parallel. The dense storages defer nothing, so one loop body serves both backends.
	auto storage = make_storage<TypeParam>(IndexArray{3, 4});
	storage.insert_zero(4); // held by both, in state 0

	std::vector<size_t> deferred_inserts;
	write_or_defer(storage.locate(4), true, deferred_inserts);  // held by both
	write_or_defer(storage.locate(5), true, deferred_inserts);  // held by the dense storages alone
	write_or_defer(storage.locate(6), false, deferred_inserts); // absent, and written to zero

	EXPECT_TRUE(storage.is_one(4)) << "a write to a held cell goes in place";
	if constexpr (!holds_every_cell<TypeParam>) {
		EXPECT_EQ(deferred_inserts, (std::vector<size_t>{5}))
		    << "the sparse array defers the cell it does not hold, and only that one";
		EXPECT_FALSE(storage.is_one(5)) << "the helper deferred the insert rather than making it";
	} else {
		EXPECT_TRUE(deferred_inserts.empty()) << "a dense storage writes every cell in place";
		EXPECT_TRUE(storage.is_one(5));
	}
	EXPECT_FALSE(storage.is_one(6)) << "an absent cell written to zero already reads as zero";

	for (const size_t linear_index : deferred_inserts) { storage.insert_one(linear_index); }
	EXPECT_TRUE(storage.is_one(4));
	EXPECT_TRUE(storage.is_one(5));
	EXPECT_FALSE(storage.is_one(6));
}

TYPED_TEST(HandleConformance, a_handle_does_not_survive_a_bulk_insert_or_a_zero_removal) {
	// The pointer-validity convention. A handle points into the storage, so it is good until the
	// storage is restructured, and these two calls are what restructure it. A window's lifetime
	// used to make this unsayable -- there was no pointer to keep -- so it is written down and
	// checked here instead.
	//
	// What is shown is that the convention is load-bearing and not decorative: across each call a
	// handle taken before it answers something that is no longer true, and a handle taken after it
	// is right again. Using the stale one is what the convention forbids, and no test can assert
	// that.
	auto storage = make_storage<TypeParam>(IndexArray{3, 4});
	storage.insert_zero(4);

	// A bulk insert gives a sparse storage a cell it did not hold.
	const auto before_insert = storage.locate(5);
	storage.insert_one(5);
	const auto after_insert = storage.locate(5);
	EXPECT_TRUE(after_insert.in_container) << "the cell was just inserted";
	EXPECT_NE(after_insert.cell, nullptr);
	EXPECT_TRUE(after_insert.is_one);
	if constexpr (!holds_every_cell<TypeParam>) {
		EXPECT_FALSE(before_insert.in_container)
		    << "a handle taken before the insert says the storage does not hold the cell";
	}

	// A zero removal takes one away again.
	const auto before_removal = storage.locate(4);
	ASSERT_TRUE(before_removal.in_container) << "every storage here holds a cell it was given";
	storage.remove_zeros();
	const auto after_removal = storage.locate(4);
	EXPECT_FALSE(after_removal.is_one) << "the cell was a zero either way";
	if constexpr (!holds_every_cell<TypeParam>) {
		EXPECT_FALSE(after_removal.in_container) << "the zero was reclaimed";
		EXPECT_EQ(after_removal.cell, nullptr);
	} else {
		EXPECT_TRUE(after_removal.in_container) << "a dense storage holds every cell";
	}
}

TYPED_TEST(StorageConformance,
           a_write_through_write_or_defer_leaves_the_storage_where_a_direct_write_does) {
	// The property the update's write exists for: `write_or_defer` plus the deferred inserts leave
	// every cell in the state `set_state` leaves it in, and leave a field's counters where
	// `set_state` leaves them -- a state write must not disturb a counter, through that write no
	// more than through `set_state`.
	//
	// Which cells the storage *holds* is not compared, and the two paths do differ over one: a
	// zero written to a cell the sparse array does not hold stores nothing, where `set_state`
	// stores a zero. Storedness is a backend decision throughout this suite, and `empty()` is the
	// one answer that reads it. See the head of this file.
	std::mt19937_64 rng(20260828);
	for (const auto &shape : generated_shapes()) {
		auto written_directly  = make_storage<TypeParam>(shape);
		auto written_by_update = make_storage<TypeParam>(shape);
		const size_t n_cells   = written_directly.total_size_of_container_space();
		TExpectedCells expected(n_cells);

		// Several short passes rather than one long one, so that a field counts iterations between
		// them and its counters differ from cell to cell by the time they are compared.
		constexpr size_t N_PASSES = 8;
		for (size_t iteration = 0; iteration < N_PASSES; ++iteration) {
			const auto writes = random_writes(rng, n_cells, 1 + n_writes_for(n_cells) / N_PASSES);
			replay(written_directly, writes, expected);
			replay_through_write_or_defer(written_by_update, writes);
			if constexpr (FieldStorage<TypeParam>) {
				written_directly.add_to_counter(iteration);
				written_by_update.add_to_counter(iteration);
			}
		}

		for (size_t i = 0; i < n_cells; ++i) {
			// Against the direct write, and against what the script says the cell holds. The
			// second is what names the guilty path when the two disagree.
			ASSERT_EQ(written_by_update.is_one(i), written_directly.is_one(i))
			    << "cell " << i << " of shape " << shape[0] << "x" << shape[1];
			ASSERT_EQ(written_directly.is_one(i), expected.states[i] != 0)
			    << "cell " << i << " of shape " << shape[0] << "x" << shape[1];
			if constexpr (FieldStorage<TypeParam>) {
				ASSERT_EQ(counter_of(written_by_update, i), counter_of(written_directly, i))
				    << "cell " << i << " of shape " << shape[0] << "x" << shape[1];
			}
		}
	}
}

// -------------------------------------------------------------------------
// The cursor the merge joins walk
// -------------------------------------------------------------------------

/// The ones of a storage, read one cell at a time. The cursor is checked against this rather than
/// against a second cursor, so a cursor that skipped a one has something to be wrong about.
template<typename Storage> std::vector<size_t> ones_by_point_lookup(const Storage &storage) {
	std::vector<size_t> ones;
	for (size_t i = 0; i < storage.total_size_of_container_space(); ++i) {
		if (storage.is_one(i)) { ones.push_back(i); }
	}
	return ones;
}

/// The ones the cursor yields, in the order it yields them.
template<typename Storage> std::vector<size_t> ones_by_cursor(const Storage &storage) {
	std::vector<size_t> ones;
	for (auto cursor = storage.ones_cursor(); cursor.valid(); cursor.advance()) {
		ones.push_back(cursor.linear_index());
	}
	return ones;
}

TYPED_TEST(OnesCursorConformance, the_cursor_yields_exactly_the_ones_in_ascending_order) {
	// The three merge joins depend on this and on nothing else about a storage. A cursor that
	// yielded the *stored* cells would hand the caller's closed-form term a different share of the
	// same sum under each backend, which is the defect the parity gate caught first.
	std::mt19937_64 rng(20260828);
	for (const auto &shape : generated_shapes()) {
		auto storage         = make_storage<TypeParam>(shape);
		const size_t n_cells = storage.total_size_of_container_space();
		TExpectedCells expected(n_cells);
		replay(storage, random_writes(rng, n_cells, n_writes_for(n_cells)), expected);

		const auto yielded = ones_by_cursor(storage);
		ASSERT_EQ(yielded, ones_by_point_lookup(storage))
		    << "shape " << shape[0] << "x" << shape[1];
		ASSERT_TRUE(std::is_sorted(yielded.begin(), yielded.end()))
		    << "shape " << shape[0] << "x" << shape[1];
	}
}

TYPED_TEST(OnesCursorConformance, the_cursor_follows_every_mutator_with_nothing_called_by_hand) {
	// The sparse array serves the cursor from a sorted vector it keeps beside the map. That cache
	// is the storage's own business: a caller asks for a cursor and gets the ones as they are now.
	// So each mutator is followed by a cursor, and nothing else is called in between.
	auto storage = make_storage<TypeParam>(IndexArray{3, 4});
	EXPECT_TRUE(ones_by_cursor(storage).empty()) << "a storage nothing was written to holds no one";

	storage.insert_one(7);
	EXPECT_EQ(ones_by_cursor(storage), std::vector<size_t>({7})) << "after insert_one";

	storage.set_state(2, true);
	EXPECT_EQ(ones_by_cursor(storage), std::vector<size_t>({2, 7})) << "after set_state(true)";

	storage.set_state(7, false);
	EXPECT_EQ(ones_by_cursor(storage), std::vector<size_t>({2})) << "after set_state(false)";

	storage.insert_zero(2);
	EXPECT_TRUE(ones_by_cursor(storage).empty()) << "after insert_zero";

	storage.insert_one(11);
	storage.remove_zeros();
	EXPECT_EQ(ones_by_cursor(storage), std::vector<size_t>({11})) << "after remove_zeros";
}

TYPED_TEST(HandleConformance, the_ones_cursor_yields_what_a_handle_wrote) {
	// A write through a handle does not pass through the storage, so a storage that caches the
	// order of its ones cannot see the write and has to stop trusting the cache. Reading the
	// cursor between the `locate` and the write is what catches a cache that clears itself at the
	// `locate` instead.
	auto storage = make_storage<TypeParam>(IndexArray{3, 4});
	storage.insert_one(2);
	storage.insert_zero(5);

	std::vector<size_t> deferred_inserts;
	const auto handle = storage.locate(5);
	ASSERT_TRUE(handle.in_container) << "every storage here holds a cell it was given";
	ASSERT_EQ(ones_by_cursor(storage), (std::vector<size_t>{2})) << "the cursor sorts, and caches";

	write_or_defer(handle, true, deferred_inserts);
	EXPECT_TRUE(deferred_inserts.empty()) << "the cell is held, so the write goes in place";
	EXPECT_EQ(ones_by_cursor(storage), (std::vector<size_t>{2, 5}));
	EXPECT_EQ(ones_by_cursor(storage), ones_by_point_lookup(storage));
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

// One shape is enough here where the tests above update them all: resetting is a pass over the
// whole container, and there is no index arithmetic in it for a shape to catch out.
// -------------------------------------------------------------------------
// The posterior fraction is a probability
// -------------------------------------------------------------------------

/// Chain lengths the thinning factor does not divide. Below 32768 neither counter needs thinning
/// at all, so a shorter chain cannot show this: numerator and denominator agree by accident there.
/// Each of these makes both fields thin, which is what the criterion asks for.
const std::vector<size_t> &chain_lengths_that_thin() {
	static const std::vector<size_t> lengths = {32769, 65537, 100003};
	return lengths;
}

/// A field with one cell held at one for the whole chain, counted over `n_iterations` iterations
/// starting at `first_iteration`.
template<typename Field> Field field_held_at_one(size_t n_iterations, size_t first_iteration) {
	Field field(n_iterations, IndexArray{2, 3});
	field.insert_one(0);
	for (size_t k = 0; k < n_iterations; ++k) { field.add_to_counter(first_iteration + k); }
	return field;
}

TYPED_TEST(FieldConformance, every_posterior_fraction_is_a_probability) {
	// The bound, over cells that are one for all, some and none of the chain -- where the test
	// below takes the one case that pins the denominator exactly. A cell switched off part way
	// through is what a real chain mostly holds, and it must land strictly inside (0, 1).
	for (const size_t n_iterations : chain_lengths_that_thin()) {
		TypeParam field(n_iterations, IndexArray{2, 3});
		field.insert_one(0);  // one throughout
		field.insert_one(1);  // one for the first half
		field.insert_zero(2); // never a one

		for (size_t iteration = 0; iteration < n_iterations; ++iteration) {
			if (iteration == n_iterations / 2) { field.set_state(1, false); }
			field.add_to_counter(iteration);
		}

		for (size_t cell = 0; cell < 3; ++cell) {
			const double fraction = field.get_fraction_of_ones(cell);
			EXPECT_GE(fraction, 0.0) << "cell " << cell << ", n_iterations = " << n_iterations;
			EXPECT_LE(fraction, 1.0) << "cell " << cell << ", n_iterations = " << n_iterations;
		}
		EXPECT_DOUBLE_EQ(field.get_fraction_of_ones(2), 0.0) << "n_iterations = " << n_iterations;
		EXPECT_GT(field.get_fraction_of_ones(1), 0.0) << "n_iterations = " << n_iterations;
		EXPECT_LT(field.get_fraction_of_ones(1), 1.0) << "n_iterations = " << n_iterations;
	}
}

TYPED_TEST(FieldConformance, a_cell_that_is_always_one_has_a_posterior_fraction_of_exactly_one) {
	for (const size_t n_iterations : chain_lengths_that_thin()) {
		const auto field = field_held_at_one<TypeParam>(n_iterations, 0);
		EXPECT_DOUBLE_EQ(field.get_fraction_of_ones(0), 1.0) << "n_iterations = " << n_iterations;
	}
}

TYPED_TEST(FieldConformance, the_posterior_fraction_does_not_depend_on_where_the_chain_started) {
	// Production never counts from zero. TDataModel hands TMarkovField an index that keeps
	// climbing through burn-in, and burninHasFinished resets the per-cell counters without
	// resetting that index, so the main chain's first counted iteration is wherever burn-in left
	// off. A denominator computed from the chain length alone is only right when that offset
	// happens to line up with the thinning factor.
	for (const size_t first_iteration : {0u, 1u, 2u, 7u, 1000u}) {
		const auto field = field_held_at_one<TypeParam>(65537, first_iteration);
		EXPECT_DOUBLE_EQ(field.get_fraction_of_ones(0), 1.0)
		    << "first counted iteration = " << first_iteration;
	}
}

TYPED_TEST(FieldConformance, reset_counts_restarts_the_denominator_with_the_counters) {
	// burninHasFinished throws the burn-in's counts away, so the fraction afterwards has to be
	// over the main chain alone -- denominator included.
	constexpr size_t n_burnin     = 500;
	constexpr size_t n_iterations = 65537;
	TypeParam field(n_iterations, IndexArray{2, 3});
	field.insert_one(0);
	for (size_t k = 0; k < n_burnin; ++k) { field.add_to_counter(k); }
	field.reset_counts();
	for (size_t k = 0; k < n_iterations; ++k) { field.add_to_counter(n_burnin + k); }

	EXPECT_DOUBLE_EQ(field.get_fraction_of_ones(0), 1.0);
}

// The thinning factor is a divisor before it is a number: add_to_counter takes the iteration
// modulo it. ceil(n / capacity) is zero for a chain with no iterations in it, so a field sized for
// one would divide by zero on its first count -- which is why the factor carries a floor of one.
TYPED_TEST(FieldConformance, the_thinning_factor_is_at_least_one_for_any_chain_length) {
	// Zero is the case the floor exists for; the lengths above it say the floor changes nothing
	// there, on both sides of each counter's capacity.
	for (const size_t n_iterations : {0u, 1u, 300u, 32769u, 100003u}) {
		const TypeParam field(n_iterations, IndexArray{2, 3});
		EXPECT_GE(field.get_thinning_factor(), 1u) << "n_iterations = " << n_iterations;
	}

	// And the field sized for no iterations survives being counted, which is what the floor is
	// really for.
	TypeParam field(0, IndexArray{2, 3});
	field.insert_one(0);
	field.add_to_counter(0);
	EXPECT_EQ(field.get_total_counts(), 1u);
	EXPECT_DOUBLE_EQ(field.get_fraction_of_ones(0), 1.0);
}

TYPED_TEST(FieldConformance, a_field_that_has_counted_nothing_reports_no_posterior) {
	TypeParam field(1000, IndexArray{2, 3});
	field.insert_one(0);
	EXPECT_DOUBLE_EQ(field.get_fraction_of_ones(0), 0.0);
}

TYPED_TEST(FieldConformance, a_write_through_write_or_defer_leaves_the_counter_alone) {
	// An update writes a state, and a state write keeps the cell's counter -- the same rule
	// set_state follows, and the opposite of the one an insert follows. The counter is what the
	// posterior is read off, so a write that reset it would throw away the chain so far.
	std::mt19937_64 rng(20260828);
	auto field           = make_storage<TypeParam>(IndexArray{4, 5});
	const size_t n_cells = field.total_size_of_container_space();
	TExpectedCells expected(n_cells);
	run_chain(field, random_script(rng, n_cells, N_ITERATIONS), expected);

	// The whole container, one row at a time and then one column at a time, each cell written to
	// the state it already holds. Both ways round, because the stride is the caller's arithmetic
	// and either way round can name the wrong cell. A column run that called a held cell absent
	// would defer the write, and the commit writes a whole new cell -- counter and all.
	std::vector<size_t> deferred;
	for (size_t row = 0; row < 4; ++row) {
		for (size_t k = 0; k < 5; ++k) {
			const size_t linear = row * 5 + k;
			write_or_defer(field, linear, field.is_one(linear), deferred);
		}
	}
	for (size_t col = 0; col < 5; ++col) {
		for (size_t k = 0; k < 4; ++k) {
			const size_t linear = k * 5 + col;
			write_or_defer(field, linear, field.is_one(linear), deferred);
		}
	}
	EXPECT_TRUE(deferred.empty()) << "no cell was written to a state it did not already hold";
	for (size_t i = 0; i < n_cells; ++i) {
		EXPECT_EQ(counter_of(field, i), expected.counts[i]) << "cell " << i;
		EXPECT_EQ(field.is_one(i), expected.states[i] != 0) << "cell " << i;
	}
}

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

/// Cell for cell: what one implementation answers, the other answers.
///
/// Every cell, and not a sample of them, because a backend that lost a cell has to be caught
/// wherever it lost it. A line of cells needs no pass of its own any more: a run is read one cell
/// at a time, so a line is cells this loop has already compared.
template<typename First, typename Second>
void expect_same_cells(First &first, Second &second, const IndexArray &shape) {
	ASSERT_EQ(first.total_size_of_container_space(), second.total_size_of_container_space());
	const size_t n_cells = first.total_size_of_container_space();

	for (size_t i = 0; i < n_cells; ++i) {
		ASSERT_EQ(first.is_one(i), second.is_one(i))
		    << "cell " << i << " of shape " << shape[0] << "x" << shape[1];
		ASSERT_EQ(first.get_multi_dimensional_index(i), second.get_multi_dimensional_index(i))
		    << "cell " << i << " of shape " << shape[0] << "x" << shape[1];
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

TEST(StorageEquivalence, the_backends_agree_when_the_handed_out_inserts_are_committed_in_bulk) {
	// The shape an update takes, end to end: every run of cells writes, hands out what it could
	// not insert, and one bulk insert commits the lot afterwards. A dense storage hands out
	// nothing, so the same loop body drives both backends.
	//
	// The runs of one pass do not overlap, exactly as one clique per column and one species leaf
	// per row do not. Overlapping runs would be a different question: a deferred insert commits
	// after a later write has set the same cell to zero, and the two backends would then be asked
	// to agree about an order the sampler never asks for.
	std::mt19937_64 rng(20260828);
	size_t n_cells_the_sparse_storages_handed_out = 0;
	for (const auto &shape : generated_shapes()) {
		const size_t n_cells = shape[0] * shape[1];
		const auto writes    = random_writes(rng, n_cells, n_writes_for(n_cells));

		TStorageZMatrix sparse_Z(shape);
		TStorageZDense dense_Z(shape);
		TStorageYMatrix sparse_Y(N_ITERATIONS, shape);
		TStorageYDense dense_Y(N_ITERATIONS, shape);
		TExpectedCells ignored(n_cells);
		replay(sparse_Z, writes, ignored);
		replay(dense_Z, writes, ignored);
		replay(sparse_Y, writes, ignored);
		replay(dense_Y, writes, ignored);

		const auto run_one_pass = [&](const std::vector<TRunShape> &requests) {
			std::vector<std::vector<uint8_t>> written;
			written.reserve(requests.size());
			for (const auto &request : requests) {
				written.push_back(random_states(rng, request.n_cells));
			}
			const auto write_and_hand_out = [&](auto &storage) {
				std::vector<std::vector<size_t>> batches;
				batches.reserve(requests.size());
				for (size_t r = 0; r < requests.size(); ++r) {
					batches.push_back(write_run(storage, requests[r], written[r]));
				}
				return batches;
			};
			const auto sparse_Z_batches = write_and_hand_out(sparse_Z);
			const auto dense_Z_batches  = write_and_hand_out(dense_Z);
			const auto sparse_Y_batches = write_and_hand_out(sparse_Y);
			const auto dense_Y_batches  = write_and_hand_out(dense_Y);

			for (const auto &batch : sparse_Z_batches) {
				n_cells_the_sparse_storages_handed_out += batch.size();
			}
			for (const auto &batch : dense_Z_batches) {
				ASSERT_TRUE(batch.empty()) << "a dense storage writes every cell in place";
			}
			for (const auto &batch : dense_Y_batches) {
				ASSERT_TRUE(batch.empty()) << "a dense storage writes every cell in place";
			}

			sparse_Z.insert_in_Z(sparse_Z_batches);
			dense_Z.insert_in_Z(dense_Z_batches);
			sparse_Y.insert_in_Y(sparse_Y_batches);
			dense_Y.insert_in_Y(dense_Y_batches);
			expect_same_cells(sparse_Z, dense_Z, shape);
			expect_same_cells(sparse_Y, dense_Y, shape);
		};

		// One pass down the rows and one along the columns, four runs each at most, because a
		// pass over every line of the largest generated shape turns the suite into a benchmark.
		std::vector<TRunShape> rows;
		for (size_t row = 0; row < std::min<size_t>(shape[0], 4); ++row) {
			rows.push_back({IndexArray{row, 0}, 1, shape[1], 1});
		}
		std::vector<TRunShape> columns;
		for (size_t col = 0; col < std::min<size_t>(shape[1], 4); ++col) {
			columns.push_back({IndexArray{0, col}, 0, shape[0], shape[1]});
		}
		run_one_pass(rows);
		run_one_pass(columns);

		if (::testing::Test::HasFailure()) { FAIL() << "shape " << shape[0] << "x" << shape[1]; }
	}

	// Without this the body above would pass just as happily on runs that handed out nothing,
	// which is the one case in which the bulk commit proves nothing.
	EXPECT_GT(n_cells_the_sparse_storages_handed_out, 0u)
	    << "no sparse storage handed an insert out, so the bulk commit was never asked to do "
	       "anything";
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

// The two fields hold the same cell, so they thin the same chain by the same factor and their
// counters agree cell for cell however long it runs. That is what makes a posterior field written
// by one comparable with a posterior field written by the other, and it is what decides which
// iterations get a trace line, so it also keeps two traces the same length.
//
// The counters used to differ by a factor of two: the sparse counter shared its 16-bit word with
// the state bit and the dense one did not. A dense posterior field written then is at twice the
// resolution of one written now.
TEST(StorageEquivalence, the_two_backends_thin_every_chain_by_the_same_factor) {
	// Around, on and well past the capacity of the counter, because the factor is a ceiling: it
	// steps at each multiple of 32767 and nowhere else.
	for (const size_t n_iterations : {0u, 1u, 300u, 32766u, 32767u, 32768u, 65534u, 100003u}) {
		const TStorageYMatrix sparse(n_iterations, IndexArray{1, 2});
		const TStorageYDense dense(n_iterations, IndexArray{1, 2});
		EXPECT_EQ(sparse.get_thinning_factor(), dense.get_thinning_factor())
		    << "n_iterations = " << n_iterations;
	}

	// And a chain long enough to thin leaves both fields holding the same counter, not merely the
	// same fraction. 65534 == 2 * 32767, so both count one iteration in two.
	constexpr size_t n_iterations = 65534;
	TStorageYMatrix sparse(n_iterations, IndexArray{1, 2});
	TStorageYDense dense(n_iterations, IndexArray{1, 2});
	ASSERT_EQ(sparse.get_thinning_factor(), 2u);
	ASSERT_EQ(dense.get_thinning_factor(), 2u);

	sparse.insert_one(0);
	dense.insert_one(0);
	for (size_t iteration = 0; iteration < n_iterations; ++iteration) {
		sparse.add_to_counter(iteration);
		dense.add_to_counter(iteration);
	}

	EXPECT_EQ(counter_of(sparse, 0), counter_of(dense, 0));
	EXPECT_EQ(sparse.get_total_counts(), dense.get_total_counts());
	// A cell that was a one for the whole chain has posterior probability one; a cell that never
	// was has zero.
	EXPECT_DOUBLE_EQ(sparse.get_fraction_of_ones(0), 1.0);
	EXPECT_DOUBLE_EQ(dense.get_fraction_of_ones(0), 1.0);
	EXPECT_DOUBLE_EQ(sparse.get_fraction_of_ones(1), 0.0);
	EXPECT_DOUBLE_EQ(dense.get_fraction_of_ones(1), 0.0);
}

// The bulk paths, which the storage concept deliberately leaves out (storage_concepts.h) and
// which therefore have nothing but this to hold them together: the deferred insert the field update
// commits after its parallel region, and the whole-space dump the per-iteration traces are written
// from. Both are named after the storage they belong to rather than after what they do, so there
// is one block per storage here rather than one templated body.
TEST(StorageEquivalence, the_backends_agree_on_the_bulk_insert_and_the_whole_space_dump) {
	std::mt19937_64 rng(20260828);
	for (const auto &shape : generated_shapes()) {
		const size_t n_cells = shape[0] * shape[1];

		// Batches, not one list: the update accumulates one vector per thread and hands them over
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

		// Both fields hand out the same packed cell, and it answers the two questions the filter
		// asks, which is the whole of what the writer needs from it.
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
