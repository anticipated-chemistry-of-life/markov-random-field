//
// The cells of one clique, as the tree that owns them addresses them.
//
// The property that matters is the mapping: node `n` of a clique is the clique's own
// multidimensional index with the tree's dimension set to `n`. It is asserted through the node
// state's own index conversion rather than by rebuilding the arithmetic beside it, so a view that
// agrees with a second copy of its own mistake still fails.
//
// The rest is what a walk needs from it: a read sees what the node state holds, a read after a
// write sees that write, a write reaches the node state -- in place or through the deferred list
// -- and no write leaves the clique it was made through.
//
// One body per storage, run against the dense and the sparse node state, because the two reach a
// cell differently: one through a handle, one through the window it keeps until it owns its cells.
// Nothing here builds a tree. A phylogeny and a node state are values.
//

#include "constants.h"
#include "coretools/algorithms.h"
#include "phylogeny_generators.h"
#include "storages/TDenseStateArray.h"
#include "storages/z_storage/TStorageZDense.h"
#include "storages/z_storage/TStorageZMatrix.h"
#include "tree/TPhylogeny.h"
#include "tree/clique/TCliqueView.h"
#include "tree/node_state_density.h"
#include "tree/node_state_draw.h"
#include "tree/node_state_shape.h"
#include "gtest/gtest.h"

#include <cstddef>
#include <cstdint>
#include <random>
#include <string>
#include <type_traits>
#include <vector>

namespace {

// -------------------------------------------------------------------------
// The concept the two headers that read a clique are written against
// -------------------------------------------------------------------------

// A view is what the node-state draw writes through and what the node-state density reads. Both
// headers name a concept and not this type, so this is the whole of what they need from it.
static_assert(node_state_draw::CliqueColumn<TCliqueView<TStorageZDense>>);
static_assert(node_state_draw::CliqueColumn<TCliqueView<TStorageZMatrix>>);
static_assert(node_state_density::CliqueStates<TCliqueView<TStorageZDense>>);
static_assert(node_state_density::CliqueStates<TCliqueView<TStorageZMatrix>>);

// -------------------------------------------------------------------------
// The shapes the properties are asserted over
// -------------------------------------------------------------------------

/// The trees every pairing below is taken from. Not balanced on purpose: a chain has a single
/// leaf, which makes one dimension of the container one cell wide, and a star has a single
/// internal node. A balanced fixture satisfies an index property by accident.
const std::vector<TPhylogeny> &trees() {
	static const std::vector<TPhylogeny> built_trees = [] {
		std::mt19937_64 rng(20260908);
		std::vector<TPhylogeny> all;
		all.push_back(build_phylogeny(phylo::random_forest(rng, 12, 1)));
		all.push_back(build_phylogeny(phylo::random_forest(rng, 14, 3)));
		all.push_back(build_phylogeny(phylo::chain(6)));
		all.push_back(build_phylogeny(phylo::star(5)));
		return all;
	}();
	return built_trees;
}

/// Every clique of the tree that owns `dimension`, as the multidimensional index it runs at: a
/// leaf in every other dimension, and 0 in its own. The enumeration TTree builds its cliques over.
std::vector<IndexArray> cliques_of(const IndexArray &leaf_counts, size_t dimension) {
	IndexArray clique_counts = leaf_counts;
	clique_counts[dimension] = 1;

	const size_t n_cliques = coretools::containerProduct(clique_counts);
	std::vector<IndexArray> cliques;
	cliques.reserve(n_cliques);
	for (size_t i = 0; i < n_cliques; ++i) {
		cliques.push_back(coretools::getSubscriptsAsArray(i, clique_counts));
	}
	return cliques;
}

/// The cell node `node` of `clique` occupies, written out here rather than taken from the header
/// under test. This is the rule the view is supposed to carry, and a test that called the view's
/// own arithmetic would agree with it whatever it did.
IndexArray expected_cell(const IndexArray &clique, size_t dimension, size_t node) {
	IndexArray cell = clique;
	cell[dimension] = node;
	return cell;
}

/// Runs `body(node_state, topology, clique, dimension)` for every clique of every tree of every
/// pairing. Each clique gets a node state of its own, so a property about what a view leaves
/// behind is asserted against an empty container.
template<typename Storage, typename Body> void for_every_clique(Body body) {
	for (const auto &first : trees()) {
		for (const auto &second : trees()) {
			const IndexArray leaf_counts{first.n_leaves(), second.n_leaves()};
			for (size_t dimension = 0; dimension < NUMBER_OF_TREES; ++dimension) {
				const TPhylogeny &owner = (dimension == 0) ? first : second;
				for (const auto &clique : cliques_of(leaf_counts, dimension)) {
					Storage Z(node_state_dimensions(leaf_counts, dimension, owner));
					body(Z, owner, clique, dimension);
				}
			}
		}
	}
}

/// The stride TTree built its clique windows with: the product of the leaf counts of every
/// dimension after the tree's own.
size_t increment_of(const IndexArray &leaf_counts, size_t dimension) {
	size_t increment = 1;
	for (size_t d = dimension + 1; d < NUMBER_OF_TREES; ++d) { increment *= leaf_counts[d]; }
	return increment;
}

// -------------------------------------------------------------------------
// One body per storage, instantiated for both node states
// -------------------------------------------------------------------------

class NodeStateNames {
public:
	template<typename Storage> static std::string GetName(int /*index*/) {
		if constexpr (std::is_same_v<Storage, TStorageZDense>) {
			return "dense_node_state";
		} else {
			return "sparse_node_state";
		}
	}
};

template<typename Storage> class CliqueView : public ::testing::Test {};
using NodeStates = ::testing::Types<TStorageZDense, TStorageZMatrix>;
TYPED_TEST_SUITE(CliqueView, NodeStates, NodeStateNames);

// -------------------------------------------------------------------------
// The mapping
// -------------------------------------------------------------------------

TYPED_TEST(CliqueView, a_node_addresses_the_cliques_cell_with_the_trees_dimension_set_to_it) {
	for_every_clique<TypeParam>(
	    [](auto &Z, const TPhylogeny &topology, const IndexArray &clique, size_t dimension) {
		    const TCliqueView<TypeParam> view(Z, topology, clique, dimension);
		    ASSERT_EQ(view.size(), topology.n_nodes());

		    for (size_t node = 0; node < topology.n_nodes(); ++node) {
			    const IndexArray expected = expected_cell(clique, dimension, node);
			    EXPECT_EQ(view.cell_of(node), expected) << "node " << node;
			    // And back through the node state's own index conversion, so the cell the walk
			    // really reads is that one too, and not a second copy of the view's arithmetic.
			    EXPECT_EQ(Z.get_multi_dimensional_index(view.linear_index(node)), expected)
			        << "node " << node;
		    }
	    });
}

TYPED_TEST(CliqueView, a_view_addresses_the_cells_the_window_over_the_same_clique_addressed) {
	// The clique's cells are a strided run, and the window over that run is what a tree used to
	// walk. This is what says the change of address is not a change of chain. It goes with the
	// window itself.
	for (const auto &first : trees()) {
		for (const auto &second : trees()) {
			const IndexArray leaf_counts{first.n_leaves(), second.n_leaves()};
			for (size_t dimension = 0; dimension < NUMBER_OF_TREES; ++dimension) {
				const TPhylogeny &owner = (dimension == 0) ? first : second;
				for (const auto &clique : cliques_of(leaf_counts, dimension)) {
					TypeParam Z(node_state_dimensions(leaf_counts, dimension, owner));
					const TCliqueView<TypeParam> view(Z, owner, clique, dimension);

					auto window =
					    Z.open_window(expected_cell(clique, dimension, 0), owner.n_nodes(),
						              increment_of(leaf_counts, dimension));
					ASSERT_EQ(view.size(), window.size());
					for (size_t node = 0; node < owner.n_nodes(); ++node) {
						EXPECT_EQ(view.linear_index(node), window.linear_index(node))
						    << "node " << node;
					}
					// Nothing was written, and the window must not reach its storage from here.
					const std::vector<size_t> inserts = window.take_buffered_inserts();
					EXPECT_TRUE(inserts.empty());
				}
			}
		}
	}
}

// -------------------------------------------------------------------------
// Reading
// -------------------------------------------------------------------------

TYPED_TEST(CliqueView, a_view_reads_the_state_the_node_state_holds) {
	for_every_clique<TypeParam>(
	    [](auto &Z, const TPhylogeny &topology, const IndexArray &clique, size_t dimension) {
		    // Written straight into the node state, so the view is read against the container
		    // rather than against itself.
		    std::vector<bool> expected(topology.n_nodes(), false);
		    for (size_t node = 0; node < topology.n_nodes(); node += 3) {
			    Z.insert_one(
			        Z.get_linear_index_in_container_space(expected_cell(clique, dimension, node)));
			    expected[node] = true;
		    }

		    const TCliqueView<TypeParam> view(Z, topology, clique, dimension);
		    for (size_t node = 0; node < topology.n_nodes(); ++node) {
			    EXPECT_EQ(view.is_one(node), expected[node]) << "node " << node;
		    }
	    });
}

TYPED_TEST(CliqueView, a_view_shows_its_own_write_to_a_later_read) {
	// What a post-order walk needs: a parent reads the states its children were just given, and
	// gets them whether or not the node state could take the write in place.
	for_every_clique<TypeParam>(
	    [](auto &Z, const TPhylogeny &topology, const IndexArray &clique, size_t dimension) {
		    TCliqueView<TypeParam> view(Z, topology, clique, dimension);

		    for (const size_t node : topology.internal_nodes()) {
			    view.set_state(node, node % 2 == 0);
		    }
		    for (const size_t node : topology.internal_nodes()) {
			    EXPECT_EQ(view.is_one(node), node % 2 == 0) << "node " << node;
		    }
		    const std::vector<size_t> deferred = view.take_deferred_inserts();
		    (void)deferred;
	    });
}

// -------------------------------------------------------------------------
// Writing
// -------------------------------------------------------------------------

TYPED_TEST(CliqueView, a_write_leaves_the_node_state_where_a_direct_write_leaves_it) {
	for_every_clique<TypeParam>([](auto &Z, const TPhylogeny &topology, const IndexArray &clique,
	                               size_t dimension) {
		std::vector<uint8_t> expected(Z.total_size_of_container_space(), 0);

		std::vector<size_t> deferred;
		{
			TCliqueView<TypeParam> view(Z, topology, clique, dimension);
			for (const size_t node : topology.internal_nodes()) {
				const bool state = node % 3 != 0;
				view.set_state(node, state);
				const size_t linear =
				    Z.get_linear_index_in_container_space(expected_cell(clique, dimension, node));
				expected[linear] = static_cast<uint8_t>(state);
			}
			deferred = view.take_deferred_inserts();
		}
		// The one exit a view inside a parallel region may take: the inserts it could not make are
		// handed out and committed once the region ends.
		Z.insert_in_Z({deferred});

		EXPECT_EQ(whole_space_states<uint8_t>(Z), expected);
	});
}

TYPED_TEST(CliqueView, a_view_defers_nothing_when_it_wrote_no_new_one) {
	for_every_clique<TypeParam>(
	    [](auto &Z, const TPhylogeny &topology, const IndexArray &clique, size_t dimension) {
		    TCliqueView<TypeParam> view(Z, topology, clique, dimension);
		    EXPECT_TRUE(view.take_deferred_inserts().empty()) << "a view that wrote nothing";

		    // A cell the node state does not hold already reads as zero, so writing one there
		    // stores a cell to say what its absence says.
		    for (const size_t node : topology.internal_nodes()) { view.set_state(node, false); }
		    EXPECT_TRUE(view.take_deferred_inserts().empty()) << "a view that wrote only zeros";
	    });
}

TYPED_TEST(CliqueView, a_view_writes_no_cell_outside_its_own_clique) {
	for_every_clique<TypeParam>(
	    [](auto &Z, const TPhylogeny &topology, const IndexArray &clique, size_t dimension) {
		    std::vector<size_t> deferred;
		    {
			    TCliqueView<TypeParam> view(Z, topology, clique, dimension);
			    for (const size_t node : topology.internal_nodes()) { view.set_state(node, true); }
			    deferred = view.take_deferred_inserts();
		    }

		    // Every deferred insert names a cell of this clique, before it is even applied. This is
		    // what makes "no clique writes another clique's cells" checkable rather than asserted.
		    for (const size_t linear : deferred) {
			    const IndexArray cell = Z.get_multi_dimensional_index(linear);
			    EXPECT_EQ(cell, expected_cell(clique, dimension, cell[dimension]));
			    EXPECT_LT(cell[dimension], topology.n_nodes());
		    }

		    Z.insert_in_Z({deferred});
		    for (size_t linear = 0; linear < Z.total_size_of_container_space(); ++linear) {
			    if (!Z.is_one(linear)) { continue; }
			    const IndexArray cell = Z.get_multi_dimensional_index(linear);
			    EXPECT_EQ(cell, expected_cell(clique, dimension, cell[dimension]))
			        << "linear index " << linear;
			    EXPECT_FALSE(topology.is_leaf(cell[dimension])) << "linear index " << linear;
		    }
	    });
}

#ifndef NDEBUG
TYPED_TEST(CliqueView, a_view_refuses_a_write_to_a_leaf) {
	// The walk assigns internal nodes. A leaf's state is drawn with the field and the other tree's
	// leaf, as one block, and not here. The view is the only place left that can catch a write to
	// the wrong block. The assertion inverts when the walk covers leaves.
	for_every_clique<TypeParam>(
	    [](auto &Z, const TPhylogeny &topology, const IndexArray &clique, size_t dimension) {
		    TCliqueView<TypeParam> view(Z, topology, clique, dimension);
		    for (const size_t leaf : topology.leaves()) {
			    EXPECT_ANY_THROW(view.set_state(leaf, true)) << "leaf " << leaf;
			    // A write of the state the cell already carries is dropped, but not before the
			    // block it lands in is checked.
			    EXPECT_ANY_THROW(view.set_state(leaf, false)) << "leaf " << leaf;
		    }
		    EXPECT_TRUE(view.take_deferred_inserts().empty());
	    });
}
#endif

} // namespace
