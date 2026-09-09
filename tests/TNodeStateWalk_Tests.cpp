//
// The Gibbs walk over one clique's node state.
//
// The property that matters is that it redraws each node from its full conditional: a root from the
// stationary distribution and its children, every other internal node from its parent and its
// children, and a leaf from its parent and from the link. The first three terms are the density
// tree/node_state_density.h scores, and the fourth is the seam the caller binds. The suite
// therefore takes the threshold each node is drawn against from the density and the link together,
// rather than restating the arithmetic the walk already carries.
//
// The walk covers every node, leaves included: a tree field is the leaf block of that tree's node
// state, and each tree draws its own (ADR-0005). So this suite asserts that a leaf moves, and that
// it moves with what the link says about it.
//
// A phylogeny, a transition grid and a link are values, so nothing here builds a tree, a field or a
// chain.
//

#include "clique_columns.h"
#include "phylogeny_generators.h"
#include "tree/TPhylogeny.h"
#include "tree/branch/TBinGrid.h"
#include "tree/branch/TTransitionGrid.h"
#include "tree/node_state_density.h"
#include "tree/node_state_walk.h"
#include "written_uniforms.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <random>
#include <utility>
#include <vector>

namespace {

using clique::TBinPerBranch;
using clique::TColumn;
using clique::TNoLink;
using clique::TOneBin;
using clique::TWrittenLink;
using clique::uniforms_for;
using node_state_density::log_density_of_clique;
using node_state_walk::update_clique;
using phylo::edge;

constexpr size_t N_BINS = 5;

/// A stream that records which cells were asked for a uniform. A node the walk visits draws
/// exactly one, so the record is the visit list. It sits here and not in the column, because what
/// a walk asked of its stream is the stream's own business.
class TRecordingUniforms {
private:
	uniforms::TWrittenUniforms _values;
	mutable std::vector<size_t> _asked;

public:
	explicit TRecordingUniforms(uniforms::TWrittenUniforms values) : _values(std::move(values)) {}

	[[nodiscard]] double at(size_t linear_index) const {
		_asked.push_back(linear_index);
		return _values.at(linear_index);
	}

	/// The linear indices a walk asked for, in the order it asked.
	[[nodiscard]] const std::vector<size_t> &asked() const { return _asked; }
};

/// How far either side of a threshold a uniform is placed. The density sums a term per node and
/// the walk sums the two or three terms that name one node, so the two derive one threshold by two
/// routes and agree to a handful of bits rather than exactly. This is wide enough to swallow that
/// and still ten million times tighter than any difference a wrong term would make.
constexpr double CLEAR_OF_THE_THRESHOLD = 1e-9;

/// The log density of a whole configuration, the link included: one term per node, plus what the
/// link says about each leaf's state.
///
/// This is what the walk's conditionals have to come from. The link enters exactly once per leaf,
/// which is what makes an internal node's conditional free of it.
template<typename Bins, typename Link>
double log_joint(const TPhylogeny &topology, const TTransitionGrid &process, const Bins &bin_of,
                 const TColumn &column, const Link &link) {
	double sum = log_density_of_clique(topology, process, column, bin_of);
	for (const size_t leaf : topology.leaves()) {
		sum += std::log(link.prob_of_leaf_states(leaf)[column.is_one(leaf) ? 1 : 0]);
	}
	return sum;
}

/// The probability the joint gives node `node` of being a one, with every other node held where
/// `column` holds it. This is the full conditional the walk draws from, derived from the density
/// and the link and not from the walk.
template<typename Bins, typename Link>
double conditional_probability_of_one(const TPhylogeny &topology, const TTransitionGrid &process,
                                      const Bins &bin_of, const TColumn &column, size_t node,
                                      const Link &link) {
	TColumn scored(topology.n_nodes());
	scored.set_from_mask(column.mask());

	scored.set_state(node, false);
	const double log_zero = log_joint(topology, process, bin_of, scored, link);
	scored.set_state(node, true);
	const double log_one = log_joint(topology, process, bin_of, scored, link);

	return 1.0 / (1.0 + std::exp(log_zero - log_one));
}

/// A link that depends on the clique and on the leaf, so that a clique reading another clique's
/// link, or a leaf reading another leaf's, moves the chain.
TWrittenLink link_of_clique(const TPhylogeny &topology, size_t clique) {
	TWrittenLink link(topology.n_leaves());
	for (const size_t leaf : topology.leaves()) {
		const double at_one = 0.2 + 0.1 * static_cast<double>((clique + 3 * leaf) % 6U);
		link.set(leaf, {1.0 - at_one, at_one});
	}
	return link;
}

/// Walks a run of cliques over one stream, in the order given and on the number of threads given.
/// Clique `c` occupies the cells at offset `c`, so no two cliques share a cell or a uniform. The
/// answer is one configuration per clique, in clique order whatever order they were walked in.
std::vector<size_t> walk_the_cliques(const TPhylogeny &topology, const TTransitionGrid &process,
                                     const TOneBin &bins,
                                     const uniforms::TWrittenUniforms &uniforms,
                                     const std::vector<size_t> &order, int n_threads) {
	const size_t n_cliques = order.size();
	std::vector<size_t> masks(n_cliques, 0);

#pragma omp parallel for num_threads(n_threads) schedule(dynamic) default(none)                    \
    shared(topology, process, bins, uniforms, order, masks, n_cliques)
	for (size_t k = 0; k < n_cliques; ++k) {
		const size_t clique = order[k];
		TColumn column(topology.n_nodes(), clique, n_cliques);
		update_clique(topology, process, bins, uniforms, link_of_clique(topology, clique), column);
		masks[clique] = column.mask();
	}
	return masks;
}

// -------------------------------------------------------------------------
// Which nodes the walk covers
// -------------------------------------------------------------------------

/// Every node is visited exactly once, leaves included. A tree draws its own tree field, which is
/// the leaf block of this very run of cells (ADR-0005).
TEST(NodeStateWalk, visits_every_node_once) {
	const TBinGrid grid(N_BINS);
	const TTransitionGrid process(0.4, 1.0, grid);
	const TPhylogeny topology =
	    build_phylogeny({edge("a", "inner"), edge("b", "inner"), edge("inner", "root"),
	                     edge("c", "root"), edge("d", "second_root")});
	ASSERT_GT(topology.n_leaves(), 0u);
	ASSERT_EQ(topology.n_roots(), 2u);

	constexpr size_t OFFSET = 3;
	constexpr size_t STRIDE = 2;
	TColumn column(topology.n_nodes(), OFFSET, STRIDE);
	const TRecordingUniforms uniforms(uniforms_for(topology, OFFSET, STRIDE, 0.5));
	update_clique(topology, process, TOneBin{1}, uniforms, TNoLink{}, column);

	std::vector<size_t> asked = uniforms.asked();
	std::sort(asked.begin(), asked.end());
	std::vector<size_t> expected;
	for (size_t node = 0; node < topology.n_nodes(); ++node) {
		expected.push_back(OFFSET + node * STRIDE);
	}
	std::sort(expected.begin(), expected.end());
	EXPECT_EQ(asked, expected);
}

/// The link is asked about every leaf, once each, and about no other node. An internal node's
/// conditional does not name it.
TEST(NodeStateWalk, asks_the_link_about_every_leaf_and_no_other_node) {
	const TBinGrid grid(N_BINS);
	const TTransitionGrid process(0.4, 1.0, grid);
	const TPhylogeny topology =
	    build_phylogeny({edge("a", "inner"), edge("b", "inner"), edge("inner", "root"),
	                     edge("c", "root")});

	/// A link that records which nodes it was asked about.
	class TRecordingLink {
	public:
		mutable std::vector<size_t> asked;
		[[nodiscard]] std::array<double, 2> prob_of_leaf_states(size_t leaf) const {
			asked.push_back(leaf);
			return {1.0, 1.0};
		}
	};

	TColumn column(topology.n_nodes());
	const TRecordingLink link;
	update_clique(topology, process, TOneBin{1},
	              uniforms::TWrittenUniforms(topology.n_nodes(), 0.5), link, column);

	std::vector<size_t> asked = link.asked;
	std::sort(asked.begin(), asked.end());
	const std::vector<size_t> leaves(topology.leaves().begin(), topology.leaves().end());
	EXPECT_EQ(asked, leaves);
}

// -------------------------------------------------------------------------
// It draws the conditional the density and the link state
// -------------------------------------------------------------------------

/// A root is drawn from the stationary distribution and from its children. The threshold comes
/// from the density, so this pins the two terms together rather than restating either.
TEST(NodeStateWalk, draws_a_root_from_the_stationary_distribution_and_its_children) {
	const TBinGrid grid(N_BINS);
	const TTransitionGrid process(0.25, 1.1, grid);
	// One root over two leaves. The leaves are drawn first and the root reads what they were
	// given, so a link that says nothing holds them where this test puts them.
	const TPhylogeny topology = build_phylogeny({edge("a", "root"), edge("b", "root")});
	const size_t root         = topology.index_of("root");
	const TOneBin bins{2};

	for (const size_t leaves_mask : {0b00U, 0b01U, 0b11U}) {
		SCOPED_TRACE("leaves " + std::to_string(leaves_mask));
		TColumn held(topology.n_nodes());
		held.set_from_mask(leaves_mask);
		const double threshold =
		    conditional_probability_of_one(topology, process, bins, held, root, TNoLink{});
		ASSERT_GT(threshold, 0.0);
		ASSERT_LT(threshold, 1.0);

		for (const auto &[uniform, expected] :
		     {std::pair{threshold * (1.0 - CLEAR_OF_THE_THRESHOLD), true},
		      std::pair{threshold * (1.0 + CLEAR_OF_THE_THRESHOLD), false}}) {
			TColumn column(topology.n_nodes());
			column.set_from_mask(leaves_mask);
			uniforms::TWrittenUniforms uniforms(topology.n_nodes(), 0.5);
			// Each leaf is held where it started, so the root reads the states this test set.
			for (const size_t leaf : topology.leaves()) {
				uniforms.set(column.linear_index(leaf),
				             column.is_one(leaf) ? 0.0 : std::nextafter(1.0, 0.0));
			}
			uniforms.set(column.linear_index(root), uniform);
			update_clique(topology, process, bins, uniforms, TNoLink{}, column);
			EXPECT_EQ(column.is_one(root), expected);
		}
	}
}

/// A non-root internal node is drawn from its parent and from its children. The parent's state is
/// the one the update started from, because the parent comes after this node in canonical order.
TEST(NodeStateWalk, draws_a_non_root_from_its_parent_and_its_children) {
	const TBinGrid grid(N_BINS);
	const TTransitionGrid process(0.35, 0.8, grid);
	const TPhylogeny topology =
	    build_phylogeny({edge("a", "inner"), edge("b", "inner"), edge("inner", "root")});
	const size_t inner = topology.index_of("inner");
	const size_t root  = topology.index_of("root");
	const TOneBin bins{3};

	for (const bool root_state : {false, true}) {
		for (const size_t leaves_mask : {0b00U, 0b01U, 0b11U}) {
			SCOPED_TRACE("root " + std::to_string(static_cast<int>(root_state)) + ", leaves " +
			             std::to_string(leaves_mask));
			TColumn held(topology.n_nodes());
			held.set_from_mask(leaves_mask);
			held.set_state(root, root_state);
			const double threshold =
			    conditional_probability_of_one(topology, process, bins, held, inner, TNoLink{});
			ASSERT_GT(threshold, 0.0);
			ASSERT_LT(threshold, 1.0);

			for (const auto &[uniform, expected] :
			     {std::pair{threshold * (1.0 - CLEAR_OF_THE_THRESHOLD), true},
			      std::pair{threshold * (1.0 + CLEAR_OF_THE_THRESHOLD), false}}) {
				TColumn column(topology.n_nodes());
				column.set_from_mask(leaves_mask);
				column.set_state(root, root_state);
				uniforms::TWrittenUniforms uniforms(topology.n_nodes(), 0.5);
				// The leaves and the root are held where they started, so the inner node reads
				// the states this update began with.
				for (const size_t leaf : topology.leaves()) {
					uniforms.set(column.linear_index(leaf),
					             column.is_one(leaf) ? 0.0 : std::nextafter(1.0, 0.0));
				}
				uniforms.set(column.linear_index(root),
				             root_state ? 0.0 : std::nextafter(1.0, 0.0));
				uniforms.set(column.linear_index(inner), uniform);
				update_clique(topology, process, bins, uniforms, TNoLink{}, column);
				ASSERT_EQ(column.is_one(root), root_state);
				EXPECT_EQ(column.is_one(inner), expected);
			}
		}
	}
}

/// A leaf is drawn from its parent and from the link, and from nothing else. A leaf has no
/// children, so those two terms are its whole conditional.
TEST(NodeStateWalk, draws_a_leaf_from_its_parent_and_the_link) {
	const TBinGrid grid(N_BINS);
	const TTransitionGrid process(0.3, 0.9, grid);
	const TPhylogeny topology =
	    build_phylogeny({edge("a", "inner"), edge("b", "inner"), edge("inner", "root")});
	const size_t leaf  = topology.index_of("a");
	const size_t inner = topology.index_of("inner");
	const TOneBin bins{2};

	for (const bool parent_state : {false, true}) {
		for (const auto &below : {std::array<double, 2>{1.0, 1.0},
		                          std::array<double, 2>{0.9, 0.1},
		                          std::array<double, 2>{0.05, 0.95}}) {
			SCOPED_TRACE("parent " + std::to_string(static_cast<int>(parent_state)) + ", link " +
			             std::to_string(below[1]));
			TWrittenLink link(topology.n_leaves());
			link.set(leaf, below);

			TColumn held(topology.n_nodes());
			held.set_state(inner, parent_state);
			const double threshold =
			    conditional_probability_of_one(topology, process, bins, held, leaf, link);
			ASSERT_GT(threshold, 0.0);
			ASSERT_LT(threshold, 1.0);

			for (const auto &[uniform, expected] :
			     {std::pair{threshold * (1.0 - CLEAR_OF_THE_THRESHOLD), true},
			      std::pair{threshold * (1.0 + CLEAR_OF_THE_THRESHOLD), false}}) {
				TColumn column(topology.n_nodes());
				column.set_state(inner, parent_state);
				uniforms::TWrittenUniforms uniforms(topology.n_nodes(), 0.5);
				uniforms.set(column.linear_index(leaf), uniform);
				update_clique(topology, process, bins, uniforms, link, column);
				EXPECT_EQ(column.is_one(leaf), expected);
			}
		}
	}
}

/// A link that all but rules a leaf state out keeps it out, whatever the parent holds and wherever
/// the uniform falls inside the interval.
///
/// The link cannot say a state is impossible: the error probability is strictly inside (0, 0.5), so
/// the link table reaches neither 0 nor 1. What it can do is put the odds far enough out that the
/// parent term cannot move them back, and this is the leaf term earning its place in the
/// conditional. The two ends of the interval are left out on purpose -- a uniform of exactly 0
/// takes state 1 against any probability above 0, which is the inverse-CDF draw doing its job
/// rather than the link failing to be heard.
TEST(NodeStateWalk, a_link_that_all_but_rules_a_state_out_keeps_it_out) {
	const TBinGrid grid(N_BINS);
	const TTransitionGrid process(0.4, 1.0, grid);
	const TPhylogeny topology =
	    build_phylogeny({edge("a", "inner"), edge("b", "inner"), edge("inner", "root")});
	const size_t leaf  = topology.index_of("a");
	const size_t inner = topology.index_of("inner");

	for (const bool forced : {false, true}) {
		TWrittenLink link(topology.n_leaves());
		link.set(leaf, forced ? std::array<double, 2>{1e-12, 1.0}
		                      : std::array<double, 2>{1.0, 1e-12});
		for (const bool parent_state : {false, true}) {
			for (const double uniform : {0.01, 0.25, 0.5, 0.75, 0.99}) {
				TColumn column(topology.n_nodes());
				column.set_state(inner, parent_state);
				uniforms::TWrittenUniforms uniforms(topology.n_nodes(), 0.5);
				uniforms.set(column.linear_index(leaf), uniform);
				update_clique(topology, process, TOneBin{1}, uniforms, link, column);
				EXPECT_EQ(column.is_one(leaf), forced)
				    << "parent " << parent_state << ", uniform " << uniform;
			}
		}
	}
}

/// Every node, whatever its kind and whatever the configuration around it, is drawn against the
/// threshold the density and the link name together. One assertion over the whole walk, so a node
/// scored against the wrong neighbours has nowhere to hide.
TEST(NodeStateWalk, draws_every_node_against_the_joint_threshold) {
	const TBinGrid grid(N_BINS);
	const TTransitionGrid process(0.3, 1.2, grid);
	const TPhylogeny topology = build_phylogeny(
	    {edge("a", "inner"), edge("b", "inner"), edge("inner", "root"), edge("c", "root")});
	// A bin per branch, so a node scored against a neighbour's branch reads a different matrix.
	const TBinPerBranch bins{N_BINS};
	// A link that differs from leaf to leaf, so a leaf reading another leaf's link moves.
	const TWrittenLink link = link_of_clique(topology, 1);

	std::mt19937_64 rng(20260909);
	std::uniform_int_distribution<size_t> any_mask(0, (size_t{1} << topology.n_nodes()) - 1);

	for (size_t replicate = 0; replicate < 40; ++replicate) {
		const size_t start = any_mask(rng);
		// The walk as it runs, and the states it leaves behind.
		TColumn walked(topology.n_nodes());
		walked.set_from_mask(start);
		uniforms::TWrittenUniforms uniforms(topology.n_nodes());
		uniforms.fill_from(rng);
		update_clique(topology, process, bins, uniforms, link, walked);

		// The same walk, node by node, with each threshold taken from the joint over the states
		// standing at that point. The two agree only if the walk scores each node against exactly
		// its parent, its children and -- for a leaf -- its own link.
		TColumn scored(topology.n_nodes());
		scored.set_from_mask(start);
		for (size_t node = 0; node < topology.n_nodes(); ++node) {
			const double threshold =
			    conditional_probability_of_one(topology, process, bins, scored, node, link);
			scored.set_state(node, uniforms.at(scored.linear_index(node)) < threshold);
		}

		for (size_t node = 0; node < topology.n_nodes(); ++node) {
			EXPECT_EQ(walked.is_one(node), scored.is_one(node))
			    << "node " << topology.id_of(node) << ", replicate " << replicate;
		}
	}
}

// -------------------------------------------------------------------------
// What the states may not depend on
// -------------------------------------------------------------------------

/// A cell draws the uniform its own position names (ADR-0007), so moving one cell's uniform moves
/// that node and nothing above it that does not read it.
TEST(NodeStateWalk, a_cell_draws_the_uniform_its_own_linear_index_names) {
	const TBinGrid grid(N_BINS);
	const TTransitionGrid process(0.4, 1.0, grid);
	// Two roots, so a root can be moved without any other node reading it. The leaves below them
	// are drawn first, and read the states this update started from.
	const TPhylogeny topology = build_phylogeny({edge("a", "r1"), edge("b", "r2")});
	constexpr size_t OFFSET   = 5;
	constexpr size_t STRIDE   = 3;
	const size_t moved        = *topology.roots().begin();

	TColumn before(topology.n_nodes(), OFFSET, STRIDE);
	uniforms::TWrittenUniforms uniforms = uniforms_for(topology, OFFSET, STRIDE, 0.5);
	update_clique(topology, process, TOneBin{1}, uniforms, TNoLink{}, before);

	TColumn after(topology.n_nodes(), OFFSET, STRIDE);
	uniforms.set(before.linear_index(moved), before.is_one(moved) ? std::nextafter(1.0, 0.0) : 0.0);
	update_clique(topology, process, TOneBin{1}, uniforms, TNoLink{}, after);

	EXPECT_NE(after.is_one(moved), before.is_one(moved));
	for (size_t node = 0; node < topology.n_nodes(); ++node) {
		if (node == moved) { continue; }
		EXPECT_EQ(after.is_one(node), before.is_one(node))
		    << "node " << topology.id_of(node) << " moved with another cell's uniform";
	}
}

/// The states a run of cliques is given do not move with the order the cliques are walked in, nor
/// with the number of threads that walk them. Each clique reads its own cells, its own uniforms and
/// its own link, so the answer is a property of the cell and not of the schedule (ADR-0007).
TEST(NodeStateWalk, the_states_do_not_move_with_traversal_order_or_thread_count) {
	const TBinGrid grid(N_BINS);
	const TTransitionGrid process(0.45, 1.3, grid);
	const TPhylogeny topology = build_phylogeny(phylo::chain(7));
	const TOneBin bins{2};

	// The stride between two cells of one clique is the number of cliques, so the run of cells is
	// interleaved exactly the way a node state interleaves them.
	constexpr size_t N_CLIQUES = 32;
	uniforms::TWrittenUniforms uniforms(N_CLIQUES * topology.n_nodes());
	std::mt19937_64 rng(20260910);
	uniforms.fill_from(rng);

	std::vector<size_t> ascending(N_CLIQUES);
	std::iota(ascending.begin(), ascending.end(), size_t{0});
	std::vector<size_t> shuffled = ascending;
	std::shuffle(shuffled.begin(), shuffled.end(), rng);

	const std::vector<size_t> expected =
	    walk_the_cliques(topology, process, bins, uniforms, ascending, 1);
	EXPECT_EQ(walk_the_cliques(topology, process, bins, uniforms, shuffled, 1), expected);
	EXPECT_EQ(walk_the_cliques(topology, process, bins, uniforms, ascending, 4), expected);
	EXPECT_EQ(walk_the_cliques(topology, process, bins, uniforms, shuffled, 4), expected);
}

} // namespace
