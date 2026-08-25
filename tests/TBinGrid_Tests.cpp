// Properties of the bin <-> branch-length correspondence.
//
// These are the C++ half of the checks that model_validation/tests/test_independent.py runs against
// the Python reference. Both implementations have to satisfy them, and until TBinGrid was pulled
// out of TTree only the Python one could be asked.

#include "process/TBinGrid.h"
#include "gtest/gtest.h"

#include <numeric>
#include <stdexcept>
#include <vector>

namespace {

// Matches ProgramOptions::BRANCH_LENGTHS_BINS, so the numbers here are the ones production runs on.
constexpr size_t N_BINS = 10;

/// A picker that hands out a fixed script of indices and refuses to be asked for more.
class TScriptedPicker {
private:
	std::vector<size_t> _indices;
	size_t _next = 0;

public:
	explicit TScriptedPicker(std::vector<size_t> indices) : _indices(std::move(indices)) {}

	size_t operator()(size_t n) {
		EXPECT_LT(_next, _indices.size()) << "the repair asked for more draws than were scripted";
		const size_t index = _indices[_next++];
		EXPECT_LT(index, n) << "scripted index is out of range for this many branches";
		return index;
	}

	[[nodiscard]] size_t draws_taken() const { return _next; }
};

/// Cycles through every branch in turn: always terminates, never needs a generator.
struct TRoundRobinPicker {
	size_t next = 0;
	size_t operator()(size_t n) { return next++ % n; }
};

/// Fails the test if the repair asks for a draw at all.
struct TNeverPicker {
	size_t operator()(size_t n) {
		ADD_FAILURE() << "the repair drew an index when it should have had nothing to do";
		return n - 1;
	}
};

size_t sum(const std::vector<size_t> &bins) {
	return std::accumulate(bins.begin(), bins.end(), size_t{0});
}

} // namespace

// --------------------------------------------------------------------------
// The grid itself
// --------------------------------------------------------------------------

TEST(BinGrid, delta_matches_the_grid_definition) {
	const TBinGrid grid(N_BINS);
	EXPECT_DOUBLE_EQ(grid.delta(), 2.0 / (double)(N_BINS + 1));
	EXPECT_EQ(grid.n_bins(), N_BINS);
	EXPECT_EQ(TBinGrid::min_bin(), 0u);
	EXPECT_EQ(grid.max_bin(), N_BINS - 1);
}

TEST(BinGrid, grid_centres_have_mean_one_at_budget) {
	// The budget is exactly the constraint "mean grid branch length is 1".
	const TBinGrid grid(N_BINS);
	const size_t n_branches  = 254;
	const double mean_bin    = (double)grid.budget(n_branches) / (double)n_branches;
	const double mean_length = grid.delta() * (mean_bin + 0.5);
	EXPECT_NEAR(mean_length, 1.0, 1e-12);
}

TEST(BinGrid, bin_of_its_own_grid_centre) {
	const TBinGrid grid(N_BINS);
	for (size_t k = 0; k < N_BINS; ++k) {
		EXPECT_EQ(grid.bin_from_length(grid.grid_branch_length(k)), k) << "bin " << k;
	}
}

TEST(BinGrid, non_positive_lengths_land_in_the_lowest_bin) {
	const TBinGrid grid(N_BINS);
	EXPECT_EQ(grid.bin_from_length(0.0), 0u);
	EXPECT_EQ(grid.bin_from_length(-1.0), 0u);
}

TEST(BinGrid, lengths_past_the_grid_clamp_to_the_top) {
	const TBinGrid grid(N_BINS);
	EXPECT_EQ(grid.bin_from_length(1000.0), N_BINS - 1);
}

TEST(BinGrid, a_grid_needs_at_least_one_bin) { EXPECT_THROW((TBinGrid(0)), std::invalid_argument); }

// --------------------------------------------------------------------------
// The budget
// --------------------------------------------------------------------------

TEST(BinGrid, repair_reaches_the_budget) {
	const TBinGrid grid(N_BINS);
	for (const size_t n_branches : {size_t{2}, size_t{7}, size_t{64}, size_t{254}}) {
		// Start every branch at the top, which is as far off budget as it gets.
		std::vector<size_t> bins(n_branches, N_BINS - 1);
		TRoundRobinPicker picker;
		grid.repair_to_budget(bins, std::ref(picker));

		EXPECT_EQ(sum(bins), grid.budget(n_branches)) << "n_branches " << n_branches;
		for (const size_t bin : bins) { EXPECT_LE(bin, grid.max_bin()); }
		EXPECT_TRUE(grid.is_on_budget(bins));
	}
}

TEST(BinGrid, repair_is_a_no_op_when_already_on_budget) {
	const TBinGrid grid(N_BINS);
	std::vector<size_t> bins(100, N_BINS / 2);
	const std::vector<size_t> before = bins;

	TNeverPicker picker;
	grid.repair_to_budget(bins, picker);

	EXPECT_EQ(bins, before);
}

TEST(BinGrid, repair_redraws_when_it_picks_a_bin_at_a_boundary) {
	// Two branches, budget 10. Bin 0 is pinned at the top and cannot rise, so a draw that picks it
	// has to be spent and retried rather than silently applied to its neighbour.
	const TBinGrid grid(N_BINS);
	std::vector<size_t> bins = {N_BINS - 1, 0};
	ASSERT_EQ(sum(bins), 9u);
	ASSERT_EQ(grid.budget(bins.size()), 10u);

	TScriptedPicker picker({0, 0, 1});
	grid.repair_to_budget(bins, std::ref(picker));

	EXPECT_EQ(picker.draws_taken(), 3u) << "the two draws landing on the full bin were not spent";
	EXPECT_EQ(bins, (std::vector<size_t>{N_BINS - 1, 1}));
	EXPECT_TRUE(grid.is_on_budget(bins));
}

// --------------------------------------------------------------------------
// Reading lengths back in
// --------------------------------------------------------------------------

TEST(BinGrid, grid_centres_round_trip_through_the_read_path) {
	// Writing grid centres recovers the bins that produced them. Only true because bins on budget
	// have mean grid length exactly 1, so the normalise-to-mean-1 step is a no-op -- and so is the
	// repair, which TNeverPicker asserts rather than assumes.
	const TBinGrid grid(N_BINS);
	// Every bin appears at least once, and the twelve of them land on budget(12) = 60.
	const std::vector<size_t> bins = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 6};
	ASSERT_TRUE(grid.is_on_budget(bins));

	std::vector<double> lengths;
	lengths.reserve(bins.size());
	for (const size_t bin : bins) { lengths.push_back(grid.grid_branch_length(bin)); }

	TNeverPicker picker;
	EXPECT_EQ(grid.bins_from_lengths(lengths, picker), bins);
}

TEST(BinGrid, flat_lengths_bin_to_the_middle) {
	// A flat tree file starts every branch at the budget's mean bin.
	const TBinGrid grid(N_BINS);
	const std::vector<double> lengths(254, 0.2);

	TNeverPicker picker;
	const auto bins = grid.bins_from_lengths(lengths, picker);

	for (const size_t bin : bins) { EXPECT_EQ(bin, N_BINS / 2); }
	EXPECT_TRUE(grid.is_on_budget(bins));
}

TEST(BinGrid, a_non_positive_length_is_rejected) {
	// Roots carry no branch; they must be filtered out before they get here.
	const TBinGrid grid(N_BINS);
	TNeverPicker picker;
	EXPECT_THROW((void)grid.bins_from_lengths({1.0, 0.0, 2.0}, picker), std::invalid_argument);
	EXPECT_THROW((void)grid.bins_from_lengths({}, picker), std::invalid_argument);
}

// --------------------------------------------------------------------------
// The budget-conserving step
// --------------------------------------------------------------------------

TEST(BinGrid, every_step_conserves_the_budget_and_stays_on_the_grid) {
	const TBinGrid grid(N_BINS);
	for (size_t first = 0; first < N_BINS; ++first) {
		for (size_t second = 0; second < N_BINS; ++second) {
			for (const bool decrease : {false, true}) {
				const int direction = grid.step_direction(first, second, decrease);
				ASSERT_GE(direction, -1);
				ASSERT_LE(direction, 1);

				const int moved_first  = (int)first + direction;
				const int moved_second = (int)second - direction;

				EXPECT_EQ(moved_first + moved_second, (int)(first + second))
				    << "budget moved at (" << first << ", " << second << ")";
				EXPECT_GE(moved_first, 0);
				EXPECT_GE(moved_second, 0);
				EXPECT_LE(moved_first, (int)grid.max_bin());
				EXPECT_LE(moved_second, (int)grid.max_bin());
			}
		}
	}
}

TEST(BinGrid, a_step_stalls_only_when_both_bins_sit_on_the_same_boundary) {
	const TBinGrid grid(N_BINS);
	for (size_t first = 0; first < N_BINS; ++first) {
		for (size_t second = 0; second < N_BINS; ++second) {
			const bool both_at_bottom = first == 0 && second == 0;
			const bool both_at_top    = first == N_BINS - 1 && second == N_BINS - 1;
			const bool stalled        = grid.step_direction(first, second, false) == 0;

			EXPECT_EQ(stalled, both_at_bottom || both_at_top)
			    << "at (" << first << ", " << second << ")";
		}
	}
}

TEST(BinGrid, the_coin_decides_exactly_when_the_step_is_free) {
	const TBinGrid grid(N_BINS);
	for (size_t first = 0; first < N_BINS; ++first) {
		for (size_t second = 0; second < N_BINS; ++second) {
			const bool coin_matters = grid.step_direction(first, second, false) !=
			                          grid.step_direction(first, second, true);
			EXPECT_EQ(grid.step_is_free(first, second), coin_matters)
			    << "at (" << first << ", " << second << ")";
		}
	}
}

TEST(BinGrid, a_forced_step_moves_away_from_the_boundary) {
	const TBinGrid grid(N_BINS);
	const size_t top = N_BINS - 1;

	EXPECT_EQ(grid.step_direction(0, 5, false), 1);    // first pinned at the bottom
	EXPECT_EQ(grid.step_direction(5, 0, false), -1);   // second pinned at the bottom
	EXPECT_EQ(grid.step_direction(top, 5, false), -1); // first pinned at the top
	EXPECT_EQ(grid.step_direction(5, top, false), 1);  // second pinned at the top

	// One bin at each boundary is legal in both directions; the bottom is checked first.
	EXPECT_EQ(grid.step_direction(0, top, false), 1);
	EXPECT_EQ(grid.step_direction(top, 0, false), -1);
}
