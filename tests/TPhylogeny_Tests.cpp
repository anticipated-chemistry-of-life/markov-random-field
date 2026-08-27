// What a tree file is allowed to say, and what comes back when it says it.
//
// Until TPhylogeny was pulled out of TTree none of this could be asked: building a tree meant
// standing up three live stattools parameters and a filename, so the loader -- including the four
// branches of its edge-list walk -- had no test at all.
//
// Three layers, per issue #15: a table of malformed edge lists, a golden test on a known tree, and
// a generator asserting the structural invariants over arbitrary forests. The generator is the one
// that earns "strongly tested": the invariants hold for *all* inputs, and a fixed corpus tests
// however many trees somebody thought to write down.

#include "tree/TPhylogeny.h"

#include "coretools/Main/TError.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

TEdge edge(std::string child, std::string parent, double length = 1.0) {
	return TEdge{std::move(child), std::move(parent), length};
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

// -------------------------------------------------------------------------
// Layer 1: malformed edge lists
// -------------------------------------------------------------------------

TEST(Phylogeny, rejects_a_node_that_is_its_own_parent) {
	expect_user_error([&] { static_cast<void>(build_phylogeny({edge("a", "a")})); });
}

TEST(Phylogeny, rejects_a_zero_or_negative_branch_length) {
	expect_user_error([&] { static_cast<void>(build_phylogeny({edge("a", "root", 0.0)})); });
	expect_user_error([&] { static_cast<void>(build_phylogeny({edge("a", "root", -1.0)})); });
}

TEST(Phylogeny, rejects_a_node_claimed_as_a_child_twice) {
	expect_user_error([&] {
		static_cast<void>(build_phylogeny({edge("child", "first"), edge("child", "second")}));
	});
}

TEST(Phylogeny, rejects_a_cycle) {
	// Two edges that between them leave every node with a parent: a -> b -> a.
	expect_user_error(
	    [&] { static_cast<void>(build_phylogeny({edge("a", "b"), edge("b", "a")})); });
}

TEST(Phylogeny, rejects_a_longer_cycle) {
	expect_user_error([&] {
		static_cast<void>(build_phylogeny({edge("a", "b"), edge("b", "c"), edge("c", "a")}));
	});
}

TEST(Phylogeny, rejects_an_empty_edge_list) {
	// No nodes at all means no root, which is the first post-condition to fail.
	expect_user_error([&] { static_cast<void>(build_phylogeny({})); });
}

TEST(Phylogeny, accepts_a_forest_with_several_roots) {
	// A tree may have more than one root; this is the shape the glossary calls out.
	const auto tree = build_phylogeny({edge("a", "r1"), edge("b", "r2")});
	EXPECT_EQ(tree.n_roots(), 2u);
	EXPECT_EQ(tree.n_leaves(), 2u);
}

// -------------------------------------------------------------------------
// Layer 2: golden test on a known tree
// -------------------------------------------------------------------------

// The fixture is a three-root forest: animal -> mammal -> {human -> denisovan, chimp},
// insect -> {mosquito, fly}, and fish -> {crab, baby_shark}.
TPhylogeny golden() {
	return read_phylogeny(std::string(ACOL_TEST_DATA_DIR) + "/loading_tree.tsv");
}

TEST(PhylogenyGolden, counts) {
	const auto tree = golden();
	EXPECT_EQ(tree.n_nodes(), 11u);
	EXPECT_EQ(tree.n_roots(), 3u);
	EXPECT_EQ(tree.n_leaves(), 6u);
	EXPECT_EQ(tree.n_internal_nodes(), 5u); // roots included
	EXPECT_EQ(tree.internal_nodes_without_roots().size(), 2u);
	EXPECT_EQ(tree.n_branches(), 8u); // n_nodes - n_roots
}

TEST(PhylogenyGolden, nodes_are_in_canonical_order) {
	// Leaves, then internal non-root nodes in post-order, then roots (ADR-0004). The file lists
	// the nodes as mammal, human, chimp, insect, mosquito, animal, denisovan, fly, fish, crab,
	// baby_shark; within the leaf block and the root block that relative order survives.
	const auto tree                         = golden();
	const std::vector<std::string> expected = {"chimp",  "mosquito", "denisovan",
	                                           "fly",    "crab",     "baby_shark", // leaves
	                                           "human",  "mammal", // internal, post-order
	                                           "insect", "animal",   "fish"}; // roots
	ASSERT_EQ(tree.n_nodes(), expected.size());
	for (size_t i = 0; i < expected.size(); ++i) {
		EXPECT_EQ(tree.id_of(i), expected[i]) << "at " << i;
	}
}

TEST(PhylogenyGolden, relationships) {
	const auto tree     = golden();
	const size_t mammal = tree.index_of("mammal");
	const size_t animal = tree.index_of("animal");
	const size_t human  = tree.index_of("human");

	EXPECT_TRUE(tree.is_root(animal));
	EXPECT_FALSE(tree.is_root(mammal));
	EXPECT_EQ(tree.parent_of(mammal), animal);
	EXPECT_EQ(tree.parent_of(animal), TPhylogeny::NO_PARENT);

	// A node adopted later keeps the children it already had, in the order they arrived.
	const auto kids = tree.children_of(mammal);
	ASSERT_EQ(kids.size(), 2u);
	EXPECT_EQ(tree.id_of(kids[0]), "human");
	EXPECT_EQ(tree.id_of(kids[1]), "chimp");

	EXPECT_TRUE(tree.is_leaf(tree.index_of("denisovan")));
	EXPECT_FALSE(tree.is_leaf(human));
}

TEST(PhylogenyGolden, branch_lengths_follow_branch_order) {
	// Branch order is node order over the non-root nodes, so it is the leaf block followed by the
	// internal non-root block: chimp, mosquito, denisovan, fly, crab, baby_shark, human, mammal.
	// mammal's length arrives on the fourth line, long after mammal was created as a root.
	const auto tree                    = golden();
	const std::vector<double> expected = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.1, 0.8};
	ASSERT_EQ(tree.branch_lengths().size(), expected.size());
	for (size_t i = 0; i < expected.size(); ++i) {
		EXPECT_DOUBLE_EQ(tree.branch_lengths()[i], expected[i]) << "at branch " << i;
	}
	EXPECT_DOUBLE_EQ(tree.branch_length_of(tree.index_of("mammal")), 0.8);
	EXPECT_DOUBLE_EQ(tree.branch_length_of(tree.index_of("baby_shark")), 0.7);
}

TEST(PhylogenyGolden, index_spaces) {
	const auto tree = golden();
	EXPECT_EQ(tree.leaf_index(tree.index_of("chimp")), 0u);
	EXPECT_EQ(tree.leaf_index(tree.index_of("baby_shark")), 5u);
	EXPECT_EQ(tree.internal_index(tree.index_of("human")), 0u);
	EXPECT_EQ(tree.internal_index(tree.index_of("fish")), 4u);
	EXPECT_EQ(tree.branch_index(tree.index_of("mammal")), 7u);

	// A node has no index in a space it does not belong to.
	EXPECT_EQ(tree.internal_index(tree.index_of("chimp")), TPhylogeny::NOT_IN_SPACE);
	EXPECT_EQ(tree.leaf_index(tree.index_of("mammal")), TPhylogeny::NOT_IN_SPACE);
	EXPECT_EQ(tree.branch_index(tree.index_of("animal")), TPhylogeny::NOT_IN_SPACE);
}

TEST(PhylogenyGolden, unknown_node_is_an_error) {
	expect_user_error([&] { (void)golden().index_of("tyrannosaurus"); });
	EXPECT_FALSE(golden().contains("tyrannosaurus"));
}

// -------------------------------------------------------------------------
// Layer 2b: file shape
// -------------------------------------------------------------------------

/// A file that only exists for the duration of one test, and never inside the source tree: a
/// failing test must not leave a stray fixture behind for the next run to trip over.
class TTempFile {
private:
	std::filesystem::path _path;

public:
	TTempFile(const std::string &name, const std::string &content) {
		// A per-process salt and a counter, because two test binaries running side by side would
		// otherwise write the same path and read each other's content.
		static const auto salt = std::to_string(std::random_device{}());
		static std::atomic<unsigned> counter{0};
		_path = std::filesystem::temp_directory_path() /
		        ("acol_" + salt + "_" + std::to_string(counter++) + "_" + name);
		std::ofstream out(_path);
		// Without this, a test that could not write its own fixture reads an absent file, the
		// reader throws, and the test passes for entirely the wrong reason. gtest assertions are
		// not usable in a constructor, so this throws instead.
		if (!out.is_open()) { throw std::runtime_error("could not create temp file " + path()); }
		out << content;
	}
	~TTempFile() {
		std::error_code ignored;
		std::filesystem::remove(_path, ignored);
	}
	TTempFile(const TTempFile &)            = delete;
	TTempFile &operator=(const TTempFile &) = delete;

	[[nodiscard]] std::string path() const { return _path.string(); }
};

TEST(PhylogenyReader, rejects_a_file_without_three_columns) {
	const TTempFile file("acol_phylogeny_two_columns.tsv", "child\tparent\na\tb\n");
	expect_user_error([&] { (void)read_phylogeny(file.path()); });
}

// -------------------------------------------------------------------------
// Layer 3: invariants over arbitrary forests
// -------------------------------------------------------------------------

/// A random forest, as an edge list in arbitrary order. Every root is given a child up front, so
/// the "a root must have at least one child" post-condition is never the thing under test here.
std::vector<TEdge> random_forest(std::mt19937_64 &rng, size_t n_nodes, size_t n_roots) {
	std::vector<TEdge> edges;
	auto name = [](size_t i) { return "n" + std::to_string(i); };

	// Nodes [0, n_roots) are roots; the next n_roots nodes are one guaranteed child each.
	for (size_t i = 0; i < n_roots; ++i) { edges.push_back(edge(name(n_roots + i), name(i))); }
	// Everything after that hangs off a uniformly chosen node that already exists.
	for (size_t i = 2 * n_roots; i < n_nodes; ++i) {
		std::uniform_int_distribution<size_t> pick(0, i - 1);
		edges.push_back(edge(name(i), name(pick(rng))));
	}
	std::shuffle(edges.begin(), edges.end(), rng);
	return edges;
}

/// The order the nodes were created in, replicating the rule `TPhylogenyBuilder::add_edge`
/// follows: a new parent is appended before a new child, and a node adopted later keeps the
/// position it already had. Canonical order permutes the nodes, but preserves this relative order
/// *within* the leaf block and the root block -- which is what lets a user reorder the lines of
/// their tree file and get a correspondingly reordered output rather than an unrelated one.
std::vector<std::string> file_insertion_order(const std::vector<TEdge> &edges) {
	std::vector<std::string> order;
	std::unordered_set<std::string> seen;
	auto intern = [&](const std::string &id) {
		if (seen.insert(id).second) { order.push_back(id); }
	};
	for (const auto &e : edges) {
		// Parent first covers both cases: when the child is new too it is genuinely created second,
		// and when the child was created by an earlier row it already holds the lower position.
		intern(e.parent);
		intern(e.child);
	}
	return order;
}

void check_file_order_within_blocks(const TPhylogeny &tree, const std::vector<TEdge> &edges) {
	const auto order = file_insertion_order(edges);
	std::unordered_map<std::string, size_t> rank;
	for (size_t i = 0; i < order.size(); ++i) { rank[order[i]] = i; }
	ASSERT_EQ(rank.size(), tree.n_nodes());

	auto ascending_in_file = [&](size_t from, size_t to, const char *block) {
		for (size_t i = from + 1; i < to; ++i) {
			EXPECT_LT(rank.at(tree.id_of(i - 1)), rank.at(tree.id_of(i)))
			    << block << " block: '" << tree.id_of(i - 1) << "' and '" << tree.id_of(i)
			    << "' are out of file order";
		}
	};
	ascending_in_file(0, tree.n_leaves(), "leaf");
	ascending_in_file(tree.n_nodes() - tree.n_roots(), tree.n_nodes(), "root");
}

void check_invariants(const TPhylogeny &tree) {
	const size_t n = tree.n_nodes();

	EXPECT_EQ(tree.n_leaves() + tree.n_internal_nodes(), n) << "leaves and internals partition";
	EXPECT_EQ(tree.n_branches(), n - tree.n_roots()) << "one branch per non-root node";
	EXPECT_EQ(tree.branch_lengths().size(), tree.n_branches());
	EXPECT_EQ(tree.internal_nodes_without_roots().size(), tree.n_internal_nodes() - tree.n_roots());

	for (size_t i = 0; i < n; ++i) {
		// id -> index -> id round trips.
		EXPECT_EQ(tree.index_of(tree.id_of(i)), i);

		// A node is a root exactly when it has no parent, and a leaf exactly when it has no child.
		EXPECT_EQ(tree.is_root(i), tree.parent_of(i) == TPhylogeny::NO_PARENT);
		EXPECT_EQ(tree.is_leaf(i), tree.children_of(i).empty());

		// Every child names this node as its parent, and vice versa.
		for (const size_t child : tree.children_of(i)) { EXPECT_EQ(tree.parent_of(child), i); }
		if (!tree.is_root(i)) {
			const auto siblings = tree.children_of(tree.parent_of(i));
			EXPECT_NE(std::find(siblings.begin(), siblings.end(), i), siblings.end())
			    << "node " << i << " is not among its parent's children";
		}

		// A node has an index in a space exactly when it belongs to it.
		EXPECT_EQ(tree.leaf_index(i) != TPhylogeny::NOT_IN_SPACE, tree.is_leaf(i));
		EXPECT_EQ(tree.internal_index(i) != TPhylogeny::NOT_IN_SPACE, !tree.is_leaf(i));
		EXPECT_EQ(tree.branch_index(i) != TPhylogeny::NOT_IN_SPACE, !tree.is_root(i));
	}

	// -- the canonical ordering (ADR-0004) ---------------------------------
	// Everything below is a statement about the *guarantee*, not about the permutation that
	// produces it: what a caller may rely on when it sees a bare node index.

	const size_t n_leaves   = tree.n_leaves();
	const size_t first_root = n - tree.n_roots();
	ASSERT_LE(n_leaves, first_root) << "the leaf block and the root block overlap";

	for (size_t i = 0; i < n; ++i) {
		EXPECT_EQ(tree.is_leaf(i), i < n_leaves) << "leaf block at " << i;
		EXPECT_EQ(tree.is_root(i), i >= first_root) << "root block at " << i;

		// Non-root nodes are exactly the first two blocks, so branch index is node index.
		if (i < first_root) { EXPECT_EQ(tree.branch_index(i), i) << "branch index at " << i; }
		// A leaf's index in leaf space is its node index; an internal node's is a subtraction.
		if (i < n_leaves) {
			EXPECT_EQ(tree.leaf_index(i), i) << "leaf index at " << i;
		} else {
			EXPECT_EQ(tree.internal_index(i), i - n_leaves) << "internal index at " << i;
		}

		// Children precede parents, which is what makes a bottom-up walk a forward loop.
		if (!tree.is_root(i)) { EXPECT_GT(tree.parent_of(i), i) << "parent of " << i; }
	}

	// The middle block is internal *and* non-root: every node in it has both a parent and a child.
	for (size_t i = n_leaves; i < first_root; ++i) {
		EXPECT_FALSE(tree.is_root(i)) << "middle block node " << i << " has no parent";
		EXPECT_FALSE(tree.children_of(i).empty()) << "middle block node " << i << " has no child";
	}

	// The three blocks partition the nodes: none left over, and none counted twice.
	EXPECT_EQ(n_leaves + tree.internal_nodes_without_roots().size() + tree.n_roots(), n);

	// The category lists agree with the maps back into them.
	for (size_t k = 0; k < tree.leaves().size(); ++k) {
		EXPECT_EQ(tree.leaf_index(tree.leaves()[k]), k);
	}
	for (size_t k = 0; k < tree.internal_nodes().size(); ++k) {
		EXPECT_EQ(tree.internal_index(tree.internal_nodes()[k]), k);
	}
	for (size_t k = 0; k < tree.branches().size(); ++k) {
		EXPECT_EQ(tree.branch_index(tree.branches()[k]), k);
	}

	// Every node is reachable from a root.
	std::vector<bool> seen(n, false);
	std::vector<size_t> stack(tree.roots().begin(), tree.roots().end());
	for (const size_t r : tree.roots()) { seen[r] = true; }
	while (!stack.empty()) {
		const size_t node = stack.back();
		stack.pop_back();
		for (const size_t child : tree.children_of(node)) {
			ASSERT_FALSE(seen[child]) << "node " << child << " reached twice: not a forest";
			seen[child] = true;
			stack.push_back(child);
		}
	}
	EXPECT_EQ(std::count(seen.begin(), seen.end(), true), static_cast<long>(n));
}

TEST(PhylogenyInvariants, hold_for_arbitrary_forests) {
	std::mt19937_64 rng(20260826);
	for (size_t trial = 0; trial < 200; ++trial) {
		std::uniform_int_distribution<size_t> n_roots_dist(1, 4);
		const size_t n_roots = n_roots_dist(rng);
		std::uniform_int_distribution<size_t> n_nodes_dist(2 * n_roots, 2 * n_roots + 40);
		const size_t n_nodes = n_nodes_dist(rng);

		const auto edges = random_forest(rng, n_nodes, n_roots);
		const auto tree  = build_phylogeny(edges);
		ASSERT_EQ(tree.n_nodes(), n_nodes) << "trial " << trial;
		EXPECT_EQ(tree.n_roots(), n_roots) << "trial " << trial;
		check_invariants(tree);
		check_file_order_within_blocks(tree, edges);
		if (::testing::Test::HasFailure()) { FAIL() << "invariants broke on trial " << trial; }
	}
}

TEST(PhylogenyInvariants, hold_for_a_deep_unbalanced_chain) {
	// A chain is the shape a balanced fixture never produces.
	std::vector<TEdge> edges;
	for (size_t i = 1; i < 60; ++i) {
		edges.push_back(edge("n" + std::to_string(i), "n" + std::to_string(i - 1)));
	}
	const auto tree = build_phylogeny(edges);
	EXPECT_EQ(tree.n_roots(), 1u);
	EXPECT_EQ(tree.n_leaves(), 1u);
	check_invariants(tree);
	check_file_order_within_blocks(tree, edges);
}

TEST(PhylogenyInvariants, hold_for_a_wide_shallow_star) {
	std::vector<TEdge> edges;
	for (size_t i = 1; i < 200; ++i) { edges.push_back(edge("n" + std::to_string(i), "root")); }
	const auto tree = build_phylogeny(edges);
	EXPECT_EQ(tree.n_roots(), 1u);
	EXPECT_EQ(tree.n_leaves(), 199u);
	check_invariants(tree);
	check_file_order_within_blocks(tree, edges);
}

} // namespace
