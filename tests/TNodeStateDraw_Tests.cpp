//
// The top-down draw of one clique's node state.
//
// The property that matters is that it draws from the density tree/node_state_density.h scores: a
// root from the stationary distribution, every other node from its parent, and every node once.
// The suite therefore counts configurations over small trees and compares the frequencies with
// that density, rather than restating the arithmetic the draw already carries.
//
// The two things the walk it replaces got wrong are asserted outright: it stopped one level above
// the leaves, and it read a parent state nothing had written.
//
// A phylogeny and a transition grid are values, so nothing here builds a tree or a chain.
//

#include "phylogeny_generators.h"
#include "tree/TPhylogeny.h"
#include "tree/branch/TBinGrid.h"
#include "tree/branch/TTransitionGrid.h"
#include "tree/node_state_density.h"
#include "tree/node_state_draw.h"
#include "written_uniforms.h"
#include "gtest/gtest.h"

#include <cmath>
#include <cstddef>
#include <random>
#include <vector>

namespace {

using node_state_density::log_density_of_clique;
using node_state_draw::draw_clique;
using phylo::edge;

constexpr size_t N_BINS = 5;

/// One clique's column, as a test writes it. The linear index of a node is deliberately not the
/// node index: a clique's cells are a strided run of the node state, and the uniform a node draws
/// is named by that stride and not by the node.
class TColumn {
private:
	std::vector<bool> _states;
	size_t _offset;
	size_t _stride;

public:
	TColumn(size_t n_nodes, size_t offset, size_t stride)
	    : _states(n_nodes, false), _offset(offset), _stride(stride) {}

	[[nodiscard]] bool is_one(size_t node) const { return static_cast<bool>(_states[node]); }
	[[nodiscard]] size_t linear_index(size_t node) const { return _offset + node * _stride; }
	void set_state(size_t node, bool state) { _states[node] = state; }

	/// The configuration as one integer, bit `node` per node. What the density is enumerated by.
	[[nodiscard]] size_t mask() const {
		size_t mask = 0;
		for (size_t node = 0; node < _states.size(); ++node) {
			if (_states[node]) { mask |= size_t{1} << node; }
		}
		return mask;
	}
};

/// Every branch in the same bin, which is what a test that is not about branch lengths wants.
struct TOneBin {
	size_t bin = 0;
	size_t operator()(size_t /*node*/) const { return bin; }
};

/// A column wide enough for the uniforms a stride of `stride` asks for.
uniforms::TWrittenUniforms uniforms_for(const TPhylogeny &topology, size_t offset, size_t stride,
                                        double value) {
	return uniforms::TWrittenUniforms(offset + topology.n_nodes() * stride, value);
}

// -------------------------------------------------------------------------
// Every node, leaves included
// -------------------------------------------------------------------------

/// The whole point of the rebuilt draw. The walk it replaces stopped at the internal nodes, so a
/// tree field -- the leaf block of this very column -- was never drawn.
TEST(NodeStateDraw, draws_every_node_including_the_leaves) {
	const TBinGrid grid(N_BINS);
	const TTransitionGrid process(0.5, 1.0, grid);
	const TPhylogeny topology = build_phylogeny(
	    {edge("a", "inner"), edge("b", "inner"), edge("inner", "root"), edge("c", "root")});
	ASSERT_GT(topology.n_leaves(), 0u);

	// Every probability in this process is strictly inside (0, 1), so a uniform of 0 gives state 1
	// at every node and a uniform just below 1 gives state 0 at every node. Which node is which
	// kind therefore cannot decide the answer -- only whether the node was reached at all.
	TColumn all_ones(topology.n_nodes(), 3, 2);
	draw_clique(topology, process, TOneBin{1}, uniforms_for(topology, 3, 2, 0.0), all_ones);

	TColumn all_zeros(topology.n_nodes(), 3, 2);
	draw_clique(topology, process, TOneBin{1},
	            uniforms_for(topology, 3, 2, std::nextafter(1.0, 0.0)), all_zeros);

	for (size_t node = 0; node < topology.n_nodes(); ++node) {
		SCOPED_TRACE("node " + topology.id_of(node));
		EXPECT_TRUE(all_ones.is_one(node));
		EXPECT_FALSE(all_zeros.is_one(node));
	}
	for (const size_t leaf : topology.leaves()) {
		EXPECT_TRUE(all_ones.is_one(leaf)) << "leaf " << topology.id_of(leaf) << " was not drawn";
	}
}

// -------------------------------------------------------------------------
// It draws what the density scores
// -------------------------------------------------------------------------

/// A root is drawn from the stationary distribution, and from nothing else.
TEST(NodeStateDraw, draws_a_root_from_the_stationary_distribution) {
	const TBinGrid grid(N_BINS);
	constexpr double ALPHA = 0.25;
	const TTransitionGrid process(ALPHA, 1.0, grid);
	const TPhylogeny topology = build_phylogeny({edge("leaf", "root")});
	const size_t root         = topology.index_of("root");

	for (const auto &[uniform, expected] :
	     {std::pair{std::nextafter(ALPHA, 0.0), true}, std::pair{ALPHA, false}}) {
		TColumn column(topology.n_nodes(), 0, 1);
		uniforms::TWrittenUniforms uniforms(topology.n_nodes(), 0.5);
		uniforms.set(column.linear_index(root), uniform);
		draw_clique(topology, process, TOneBin{2}, uniforms, column);
		EXPECT_EQ(column.is_one(root), expected);
	}
}

/// Every other node is drawn from the state its parent was just given. The draw this replaces read
/// every parent as 0 instead, which is the defect this pins.
TEST(NodeStateDraw, draws_a_node_from_the_state_its_parent_was_given) {
	const TBinGrid grid(N_BINS);
	const TTransitionGrid process(0.25, 1.1, grid);
	const TPhylogeny topology = build_phylogeny({edge("leaf", "root")});
	const size_t root         = topology.index_of("root");
	const size_t leaf         = topology.index_of("leaf");
	constexpr size_t BIN      = 2;

	for (const bool root_state : {false, true}) {
		SCOPED_TRACE("root state " + std::to_string(static_cast<int>(root_state)));
		const double threshold = process.probability(BIN, root_state, /*to =*/true);
		ASSERT_GT(threshold, 0.0);
		ASSERT_LT(threshold, 1.0);

		for (const auto &[uniform, expected] :
		     {std::pair{std::nextafter(threshold, 0.0), true}, std::pair{threshold, false}}) {
			TColumn column(topology.n_nodes(), 0, 1);
			uniforms::TWrittenUniforms uniforms(topology.n_nodes(), 0.5);
			// The root's own draw, chosen so it lands on the state under test.
			uniforms.set(column.linear_index(root), root_state ? 0.0 : std::nextafter(1.0, 0.0));
			uniforms.set(column.linear_index(leaf), uniform);
			draw_clique(topology, process, TOneBin{BIN}, uniforms, column);
			ASSERT_EQ(column.is_one(root), root_state);
			EXPECT_EQ(column.is_one(leaf), expected);
		}
	}
}

/// The frequencies it draws are the density tree/node_state_density.h scores. The two derive the
/// same conditional terms in opposite directions, so agreement is a cross-check and not a
/// restatement: a branch scored twice, or a root drawn from the wrong distribution, breaks it.
TEST(NodeStateDraw, draws_the_configurations_the_clique_density_scores) {
	const TBinGrid grid(N_BINS);
	const TTransitionGrid process(0.35, 0.8, grid);
	const TPhylogeny topology =
	    build_phylogeny({edge("a", "inner"), edge("b", "inner"), edge("inner", "root")});
	const TOneBin bins{2};

	constexpr size_t N_DRAWS      = 200000;
	const size_t n_configurations = size_t{1} << topology.n_nodes();
	std::vector<size_t> counted(n_configurations, 0);

	std::mt19937_64 rng(20260907);
	for (size_t draw = 0; draw < N_DRAWS; ++draw) {
		TColumn column(topology.n_nodes(), 0, 1);
		uniforms::TWrittenUniforms uniforms(topology.n_nodes());
		uniforms.fill_from(rng);
		draw_clique(topology, process, bins, uniforms, column);
		++counted[column.mask()];
	}

	double total_scored = 0.0;
	for (size_t mask = 0; mask < n_configurations; ++mask) {
		TColumn column(topology.n_nodes(), 0, 1);
		for (size_t node = 0; node < topology.n_nodes(); ++node) {
			column.set_state(node, ((mask >> node) & 1U) != 0U);
		}
		const double scored = std::exp(log_density_of_clique(topology, process, column, bins));
		total_scored += scored;
		const double drawn = static_cast<double>(counted[mask]) / static_cast<double>(N_DRAWS);
		EXPECT_NEAR(drawn, scored, 0.005) << "configuration " << mask;
	}
	// The enumeration is complete, so nothing was drawn that the density does not score.
	EXPECT_NEAR(total_scored, 1.0, 1e-12);
}

/// A tree with several roots draws each of them independently, and the walk needs no special case:
/// canonical order puts every root above every node below it (ADR-0004).
TEST(NodeStateDraw, draws_every_root_of_a_tree_that_has_several) {
	const TBinGrid grid(N_BINS);
	constexpr double ALPHA = 0.25;
	const TTransitionGrid process(ALPHA, 1.0, grid);
	const TPhylogeny topology =
	    build_phylogeny({edge("a", "r1"), edge("b", "r2"), edge("c", "r3")});
	ASSERT_EQ(topology.n_roots(), 3u);

	TColumn column(topology.n_nodes(), 0, 1);
	uniforms::TWrittenUniforms uniforms(topology.n_nodes(), 0.5);
	// Each root reads its own uniform, so the three come out different.
	const auto roots = topology.roots();
	uniforms.set(column.linear_index(*roots.begin()), std::nextafter(ALPHA, 0.0));
	uniforms.set(column.linear_index(*(roots.begin() + 1)), ALPHA);
	uniforms.set(column.linear_index(*(roots.begin() + 2)), std::nextafter(ALPHA, 0.0));
	draw_clique(topology, process, TOneBin{2}, uniforms, column);

	EXPECT_TRUE(column.is_one(*roots.begin()));
	EXPECT_FALSE(column.is_one(*(roots.begin() + 1)));
	EXPECT_TRUE(column.is_one(*(roots.begin() + 2)));
}

/// A process that never switches copies its root all the way down, leaves included.
TEST(NodeStateDraw, a_frozen_process_copies_every_root_down_to_the_leaves) {
	const TBinGrid grid(N_BINS);
	// Slow enough that no branch on the grid switches within the tolerance below.
	const TTransitionGrid process(0.5, 1e-9, grid);
	const TPhylogeny topology = build_phylogeny(phylo::chain(6));

	std::mt19937_64 rng(20260908);
	for (size_t replicate = 0; replicate < 20; ++replicate) {
		TColumn column(topology.n_nodes(), 0, 1);
		uniforms::TWrittenUniforms uniforms(topology.n_nodes());
		uniforms.fill_from(rng);
		draw_clique(topology, process, TOneBin{N_BINS - 1}, uniforms, column);

		for (const size_t root : topology.roots()) {
			for (size_t node = 0; node < topology.n_nodes(); ++node) {
				EXPECT_EQ(column.is_one(node), column.is_one(root));
			}
		}
	}
}

// -------------------------------------------------------------------------
// Which uniform a cell reads
// -------------------------------------------------------------------------

/// A cell draws the uniform its own position names (ADR-0007). A leaf is what this is asserted on,
/// because nothing below it reads the state that changes.
TEST(NodeStateDraw, a_cell_draws_the_uniform_its_own_linear_index_names) {
	const TBinGrid grid(N_BINS);
	const TTransitionGrid process(0.4, 1.0, grid);
	const TPhylogeny topology = build_phylogeny(phylo::star(4));
	constexpr size_t OFFSET   = 5;
	constexpr size_t STRIDE   = 3;
	const size_t leaf         = *topology.leaves().begin();

	TColumn before(topology.n_nodes(), OFFSET, STRIDE);
	uniforms::TWrittenUniforms uniforms = uniforms_for(topology, OFFSET, STRIDE, 0.5);
	draw_clique(topology, process, TOneBin{1}, uniforms, before);

	// The uniform of the cell this leaf occupies, moved to the other end of the unit interval.
	TColumn after(topology.n_nodes(), OFFSET, STRIDE);
	uniforms.set(before.linear_index(leaf), before.is_one(leaf) ? std::nextafter(1.0, 0.0) : 0.0);
	draw_clique(topology, process, TOneBin{1}, uniforms, after);

	EXPECT_NE(after.is_one(leaf), before.is_one(leaf));
	for (size_t node = 0; node < topology.n_nodes(); ++node) {
		if (node == leaf) { continue; }
		EXPECT_EQ(after.is_one(node), before.is_one(node))
		    << "node " << topology.id_of(node) << " moved with another cell's uniform";
	}
}

} // namespace
