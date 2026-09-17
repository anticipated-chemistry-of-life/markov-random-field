//
// The clique-parallel-pass shape every one of TTree's own passes spells out by hand.
//
// The traversal itself does two things past what TCliqueView already carries: it gathers a
// kernel's per-clique result into clique order, and it collects every clique's deferred inserts
// into one list per clique. Both are asserted against a baseline built without the traversal, so a
// bug that made the traversal agree with itself would still fail here.
//
// The property that matters most is the one TTree's own hand-written regions did not have until
// now: the gathered result does not depend on the thread count. Nothing here builds a tree or a
// chain -- a phylogeny and a node state are values, the way node_state_walk.h and
// node_state_density.h already rely on.
//

#include "constants.h"
#include "phylogeny_generators.h"
#include "storages/z_storage/TStorageZDense.h"
#include "storages/z_storage/TStorageZSparse.h"
#include "thread_count.h"
#include "tree/TPhylogeny.h"
#include "tree/clique/TCliqueSpace.h"
#include "tree/clique/TCliqueView.h"
#include "tree/clique/clique_traversal.h"
#include "tree/node_state_shape.h"
#include "gtest/gtest.h"

#include <random>
#include <vector>

namespace {

using threads::TThreadCount;

constexpr size_t SPECIES = 0;

/// A named kernel type, so the traversal's own concept can be asserted against it directly rather
/// than only through a lambda a test never names.
struct TSumKernel {
	double operator()(size_t clique, TCliqueView<TStorageZDense> &view) const {
		return static_cast<double>(clique) + static_cast<double>(view.size());
	}
};
static_assert(clique_traversal::CliqueKernel<TSumKernel, TCliqueView<TStorageZDense>>);

/// A small, deliberately unbalanced tree: several roots, no two subtrees the same depth. A
/// balanced fixture would satisfy a clique-ordering bug by accident.
TPhylogeny small_tree() {
	std::mt19937_64 rng(20260917);
	return build_phylogeny(phylo::random_forest(rng, 11, 2));
}

/// This tree's leaf count, and the other tree's, in the two-dimension shape every test below opens
/// a clique space and a node state over.
IndexArray leaf_counts(const TPhylogeny &topology, size_t other_tree_leaves) {
	return {topology.n_leaves(), other_tree_leaves};
}

// -------------------------------------------------------------------------
// Gathering a value-returning kernel
// -------------------------------------------------------------------------

TEST(CliqueTraversal, gathers_a_value_returning_kernels_results_in_clique_order) {
	const TPhylogeny topology = small_tree();
	const IndexArray counts   = leaf_counts(topology, /*other_tree_leaves=*/5);
	TStorageZDense Z(node_state_dimensions(counts, SPECIES, topology));
	const TCliqueSpace<> space(counts, SPECIES);
	const size_t n_cliques = space.n_cliques();
	ASSERT_GT(n_cliques, 1u);

	auto view_factory = [&](size_t clique) {
		return TCliqueView<TStorageZDense>(Z, topology, space.index_of(clique), SPECIES);
	};
	// A value distinct per clique, so a result landing at the wrong index is a wrong number and
	// not a coincidence that a symmetric fixture would hide.
	auto kernel = [](size_t clique, TCliqueView<TStorageZDense> &view) {
		return static_cast<double>(clique) + static_cast<double>(view.size()) / 100.0;
	};

	const auto [gathered, inserts] = clique_traversal::run(n_cliques, view_factory, kernel);

	std::vector<double> expected(n_cliques);
	for (size_t clique = 0; clique < n_cliques; ++clique) {
		expected[clique] =
		    static_cast<double>(clique) + static_cast<double>(topology.n_nodes()) / 100.0;
	}
	EXPECT_EQ(gathered, expected);
	EXPECT_EQ(inserts.size(), n_cliques);
}

// -------------------------------------------------------------------------
// A void kernel
// -------------------------------------------------------------------------

TEST(CliqueTraversal, a_void_kernel_gathers_no_result_only_inserts) {
	const TPhylogeny topology = small_tree();
	const IndexArray counts   = leaf_counts(topology, /*other_tree_leaves=*/4);
	TStorageZSparse Z(node_state_dimensions(counts, SPECIES, topology));
	const TCliqueSpace<> space(counts, SPECIES);
	const size_t n_cliques = space.n_cliques();

	auto view_factory = [&](size_t clique) {
		return TCliqueView<TStorageZSparse>(Z, topology, space.index_of(clique), SPECIES);
	};
	auto kernel = [&topology](size_t /*clique*/, TCliqueView<TStorageZSparse> &view) {
		for (const size_t node : topology.internal_nodes()) { view.set_state(node, true); }
	};

	// The kernel returns nothing, so `run`'s result is only the inserts -- passing it straight to
	// insert_in_Z, which takes exactly that type, is what proves it rather than a name for it.
	const auto inserts = clique_traversal::run(n_cliques, view_factory, kernel);
	EXPECT_EQ(inserts.size(), n_cliques);
	Z.insert_in_Z(inserts);

	std::vector<size_t> expected(Z.total_size_of_container_space(), 0);
	for (size_t clique = 0; clique < n_cliques; ++clique) {
		IndexArray cell = space.index_of(clique);
		for (const size_t node : topology.internal_nodes()) {
			cell[SPECIES]                                          = node;
			expected[Z.get_linear_index_in_container_space(cell)] = 1;
		}
	}
	EXPECT_EQ(Z.get_full_Z_binary_vector(), expected);
}

// -------------------------------------------------------------------------
// Deferred inserts
// -------------------------------------------------------------------------

TEST(CliqueTraversal, deferred_inserts_reproduce_what_writing_each_clique_directly_would) {
	const TPhylogeny topology = small_tree();
	const IndexArray counts   = leaf_counts(topology, /*other_tree_leaves=*/6);
	const TCliqueSpace<> space(counts, SPECIES);
	const size_t n_cliques = space.n_cliques();

	auto write_kernel = [&topology](size_t clique, TCliqueView<TStorageZSparse> &view) {
		for (const size_t node : topology.internal_nodes()) {
			view.set_state(node, (node + clique) % 3 != 0);
		}
	};

	// The baseline: every clique's view opened and written one at a time, outside the traversal,
	// and every insert applied in one batch -- the same order production applies them in.
	TStorageZSparse Z_direct(node_state_dimensions(counts, SPECIES, topology));
	std::vector<std::vector<size_t>> inserts_direct(n_cliques);
	for (size_t clique = 0; clique < n_cliques; ++clique) {
		TCliqueView<TStorageZSparse> view(Z_direct, topology, space.index_of(clique), SPECIES);
		write_kernel(clique, view);
		inserts_direct[clique] = view.take_deferred_inserts();
	}
	Z_direct.insert_in_Z(inserts_direct);

	// The same writes, through the traversal.
	TStorageZSparse Z_traversal(node_state_dimensions(counts, SPECIES, topology));
	auto view_factory = [&](size_t clique) {
		return TCliqueView<TStorageZSparse>(Z_traversal, topology, space.index_of(clique), SPECIES);
	};
	const auto inserts = clique_traversal::run(n_cliques, view_factory, write_kernel);
	EXPECT_EQ(inserts.size(), n_cliques);
	Z_traversal.insert_in_Z(inserts);

	EXPECT_EQ(Z_traversal.get_full_Z_binary_vector(), Z_direct.get_full_Z_binary_vector());
}

// -------------------------------------------------------------------------
// Thread-count independence
// -------------------------------------------------------------------------

/// The check TTree's own hand-written clique-parallel regions did not have: the gathered result is
/// the same whatever the thread count. Each clique's sum is order-sensitive within that clique --
/// it walks the clique's internal nodes and folds a per-node term into a running total -- but
/// nothing is summed *across* cliques here, so a correct traversal passes this by construction.
/// What it actually guards is the traversal's own bookkeeping: that clique `c`'s result lands at
/// `gathered[c]` and nowhere else, whichever thread computed it and in whatever order the threads
/// finished.
TEST(CliqueTraversal, gathers_the_same_result_whatever_the_thread_count) {
	std::mt19937_64 rng(20260917);
	const TPhylogeny topology = build_phylogeny(phylo::random_forest(rng, 9, 2));
	// A prime clique count, so a dynamic four-way split can never land evenly and a result at the
	// wrong index is not masked by every thread happening to get the same share.
	constexpr size_t N_CLIQUES = 37;
	const IndexArray counts    = leaf_counts(topology, N_CLIQUES);
	const TCliqueSpace<> space(counts, SPECIES);
	ASSERT_EQ(space.n_cliques(), N_CLIQUES);

	auto kernel = [&topology](size_t clique, TCliqueView<TStorageZDense> &view) {
		double sum = 0.0;
		for (const size_t node : topology.internal_nodes()) {
			const bool state = (node + clique) % 2 == 0;
			view.set_state(node, state);
			sum += state ? 0.75 : 0.25;
		}
		return sum;
	};

	TStorageZDense Z_one(node_state_dimensions(counts, SPECIES, topology));
	auto view_factory_one = [&](size_t clique) {
		return TCliqueView<TStorageZDense>(Z_one, topology, space.index_of(clique), SPECIES);
	};
	std::vector<double> at_one_thread;
	{
		TThreadCount guard(1);
		at_one_thread = clique_traversal::run(N_CLIQUES, view_factory_one, kernel).first;
	}

	// A fresh storage for the second run, so the first run's writes cannot leak into it.
	TStorageZDense Z_four(node_state_dimensions(counts, SPECIES, topology));
	auto view_factory_four = [&](size_t clique) {
		return TCliqueView<TStorageZDense>(Z_four, topology, space.index_of(clique), SPECIES);
	};
	std::vector<double> at_four_threads;
	{
		TThreadCount guard(4);
		at_four_threads = clique_traversal::run(N_CLIQUES, view_factory_four, kernel).first;
	}

	EXPECT_EQ(at_one_thread, at_four_threads);
}

} // namespace
