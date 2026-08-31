#include "constants.h"
#include "storages/y_storage/TStorageYDense.h"
#include "storages/y_storage/TStorageYMatrix.h"
#include "storages/z_storage/TStorageZDense.h"
#include "gtest/gtest.h"
#include <cstddef>
#include <cstdint>
#include <vector>

//-----------------------------------
// TStorageZDense -- the dense internal state, which is the shared state array plus the bulk paths
// `Z` is asked for by name. What is asserted of the array here therefore holds of the dense field's
// state too; the field's own suite below only has to check what the counter adds.
//
// Tests use a single-row layout {1, N} unless the shape matters, so that the linear index equals
// the column.
//-----------------------------------

TEST(ZStorageDense_Tests, a_sized_array_is_all_zeros) {
	const TStorageZDense Z({2, 3});
	EXPECT_EQ(Z.total_size_of_container_space(), 6u);
	EXPECT_TRUE(Z.empty());
	for (size_t i = 0; i < 6; ++i) { EXPECT_FALSE(Z.is_one(i)); }
}

TEST(ZStorageDense_Tests, set_state_round_trips) {
	TStorageZDense Z({1, 5});
	Z.set_state(1, true);
	Z.set_state(3, true);
	EXPECT_EQ(Z.is_one(0), false);
	EXPECT_EQ(Z.is_one(1), true);
	EXPECT_EQ(Z.is_one(2), false);
	EXPECT_EQ(Z.is_one(3), true);
	EXPECT_EQ(Z.is_one(4), false);

	// setting the same state twice keeps it, and setting it back clears it
	Z.set_state(1, true);
	EXPECT_TRUE(Z.is_one(1));
	Z.set_state(1, false);
	EXPECT_FALSE(Z.is_one(1));
}

TEST(ZStorageDense_Tests, insert_one_and_insert_zero_round_trip) {
	TStorageZDense Z({1, 4});
	Z.insert_one(2);
	EXPECT_TRUE(Z.is_one(2));
	Z.insert_zero(2);
	EXPECT_FALSE(Z.is_one(2));
}

// The one place the two implementations are deliberately allowed to answer differently: sparse
// reports whether anything is stored, and an insert_zero stores something. Dense has no such
// state to report and answers whether any cell is one, which is what the caller is asking.
TEST(ZStorageDense_Tests, empty_tracks_whether_any_cell_is_one) {
	TStorageZDense Z({1, 4});
	EXPECT_TRUE(Z.empty());
	Z.insert_one(2);
	EXPECT_FALSE(Z.empty());
	Z.set_state(2, false);
	EXPECT_TRUE(Z.empty());
	Z.insert_zero(2); // stored-as-zero in the sparse matrix; still no cell is one here
	EXPECT_TRUE(Z.empty());
}

TEST(ZStorageDense_Tests, initialize_dimensions_resizes_and_clears) {
	TStorageZDense Z({1, 4});
	Z.insert_one(3);
	Z.initialize_dimensions({2, 3});
	EXPECT_EQ(Z.total_size_of_container_space(), 6u);
	EXPECT_TRUE(Z.empty());
}

TEST(ZStorageDense_Tests, index_conversion_round_trips_over_the_whole_space) {
	const TStorageZDense Z({3, 4});
	for (size_t linear = 0; linear < Z.total_size_of_container_space(); ++linear) {
		const auto multidim = Z.get_multi_dimensional_index(linear);
		EXPECT_EQ(Z.get_linear_index_in_container_space(multidim), linear);
	}
	// row-major, so the linear index of (row, col) is row * n_cols + col
	EXPECT_EQ(Z.get_multi_dimensional_index(5), (IndexArray{1, 1}));
	EXPECT_EQ(Z.get_linear_index_in_container_space(IndexArray{2, 3}), 11u);
}

TEST(ZStorageDense_Tests, remove_zeros_leaves_every_cell_where_it_is) {
	TStorageZDense Z({1, 5});
	Z.insert_one(1);
	Z.insert_zero(2);
	Z.remove_zeros();
	EXPECT_EQ(Z.total_size_of_container_space(), 5u);
	EXPECT_TRUE(Z.is_one(1));
	EXPECT_FALSE(Z.is_one(2));
}

TEST(ZStorageDense_Tests, insert_one_outside_the_container_space_throws) {
	TStorageZDense Z({2, 3}); // total size 6
	EXPECT_ANY_THROW(Z.insert_one(6));
	EXPECT_NO_THROW(Z.insert_one(5));
}

TEST(ZStorageDense_Tests, insert_zero_outside_the_container_space_throws) {
	TStorageZDense Z({2, 3}); // total size 6
	EXPECT_ANY_THROW(Z.insert_zero(6));
	EXPECT_NO_THROW(Z.insert_zero(5));
}

// increment == 1: the variable dimension is the last one -- for the sparse implementation a row
// walk, for this one a straight run of linear indices.
TEST(ZStorageDense_Tests, fill_current_state_along_the_last_dimension) {
	TStorageZDense Z({1, 6});
	Z.insert_one(1);
	Z.insert_one(3);

	std::vector<uint8_t> state;
	std::vector<uint8_t> exists;
	std::vector<size_t> linear;
	Z.fill_current_state(IndexArray{0, 0}, /*K=*/6, /*increment=*/1, state, exists, linear);

	EXPECT_EQ(linear, (std::vector<size_t>{0, 1, 2, 3, 4, 5}));
	EXPECT_EQ(state, (std::vector<uint8_t>{0, 1, 0, 1, 0, 0}));
	// dense stores the whole container space, so every cell of the clique exists
	EXPECT_EQ(exists, (std::vector<uint8_t>{1, 1, 1, 1, 1, 1}));
}

// increment > 1: the variable dimension is the first one -- a column of the matrix.
TEST(ZStorageDense_Tests, fill_current_state_along_the_first_dimension) {
	TStorageZDense Z({3, 2}); // 3 rows, 2 cols -> nCols == 2 == increment
	Z.insert_one(2);          // (row 1, col 0)
	Z.insert_one(4);          // (row 2, col 0)
	Z.insert_one(3);          // (row 1, col 1) -> not on the col-0 walk

	std::vector<uint8_t> state;
	std::vector<uint8_t> exists;
	std::vector<size_t> linear;
	Z.fill_current_state(IndexArray{0, 0}, /*K=*/3, /*increment=*/2, state, exists, linear);

	EXPECT_EQ(linear, (std::vector<size_t>{0, 2, 4}));
	EXPECT_EQ(state, (std::vector<uint8_t>{0, 1, 1}));
	EXPECT_EQ(exists, (std::vector<uint8_t>{1, 1, 1}));
}

TEST(ZStorageDense_Tests, fill_current_state_honours_the_start_index) {
	TStorageZDense Z({1, 6});
	Z.insert_one(2);
	Z.insert_one(4);
	Z.insert_one(5); // past the end of the window -> must not be reported

	std::vector<uint8_t> state;
	std::vector<uint8_t> exists;
	std::vector<size_t> linear;
	// window covers columns [2, 5)
	Z.fill_current_state(IndexArray{0, 2}, /*K=*/3, /*increment=*/1, state, exists, linear);

	EXPECT_EQ(linear, (std::vector<size_t>{2, 3, 4}));
	EXPECT_EQ(state, (std::vector<uint8_t>{1, 0, 1}));
	EXPECT_EQ(exists, (std::vector<uint8_t>{1, 1, 1}));
}

//-----------------------------------
// TStorageYDense -- the dense field: the state array above, with a 16-bit counter beside it.
//-----------------------------------

TEST(YStorageDense_Tests, state_round_trips_and_out_of_range_inserts_throw) {
	TStorageYDense Y(1000, {2, 3}); // total size 6
	EXPECT_EQ(Y.total_size_of_container_space(), 6u);
	EXPECT_TRUE(Y.empty());

	Y.insert_one(5);
	EXPECT_TRUE(Y.is_one(5));
	EXPECT_FALSE(Y.empty());
	EXPECT_EQ(Y.get_multi_dimensional_index(5), (IndexArray{1, 2}));

	Y.set_state(5, false);
	EXPECT_FALSE(Y.is_one(5));
	EXPECT_TRUE(Y.empty());

	EXPECT_ANY_THROW(Y.insert_one(6));
	EXPECT_ANY_THROW(Y.insert_zero(6));
}

TEST(YStorageDense_Tests, thinning_factor_uses_the_full_16_bit_counter) {
	// 65534 == 2 * 32767, the capacity of the sparse field's 15-bit counter. The dense counter has
	// the sixteenth bit too, so the same chain needs no thinning at all where sparse needs half.
	constexpr size_t n_iterations = 65534;
	const TStorageYDense dense(n_iterations, {1, 4});
	const TStorageYMatrix sparse(n_iterations, IndexArray{1, 4});

	EXPECT_EQ(dense.get_thinning_factor(), 1u);
	EXPECT_EQ(sparse.get_thinning_factor(), 2u);

	// The counts are what a chain produced, not what its length predicts, so before one runs there
	// is nothing to report. Running it is what the two thinning factors then differ over.
	EXPECT_EQ(dense.get_total_counts(), 0u);
	EXPECT_EQ(sparse.get_total_counts(), 0u);

	TStorageYDense dense_run(n_iterations, {1, 4});
	TStorageYMatrix sparse_run(n_iterations, IndexArray{1, 4});
	for (size_t iteration = 0; iteration < n_iterations; ++iteration) {
		dense_run.add_to_counter(iteration);
		sparse_run.add_to_counter(iteration);
	}
	EXPECT_EQ(dense_run.get_total_counts(), 65534u);  // every iteration
	EXPECT_EQ(sparse_run.get_total_counts(), 32767u); // one in two
}

TEST(YStorageDense_Tests, counter_accumulates_for_ones_only) {
	TStorageYDense Y(1000, {1, 5}); // 1000 iterations -> thinning factor 1
	ASSERT_EQ(Y.get_thinning_factor(), 1u);
	Y.insert_one(1);
	Y.insert_zero(2);

	Y.add_to_counter(0);
	Y.add_to_counter(1);
	Y.add_to_counter(2);
	EXPECT_EQ(Y.get_counter(1), 3);
	EXPECT_EQ(Y.get_counter(2), 0);

	// a cell that stops being a one stops being counted, and keeps what it has
	Y.set_state(1, false);
	Y.add_to_counter(3);
	EXPECT_EQ(Y.get_counter(1), 3);
}

TEST(YStorageDense_Tests, counter_accumulates_once_per_thinning_factor) {
	// 196605 == 3 * 65535, so one iteration in three is counted.
	TStorageYDense Y(196605, {1, 2});
	ASSERT_EQ(Y.get_thinning_factor(), 3u);
	Y.insert_one(0);

	for (size_t iteration = 0; iteration < 9; ++iteration) { Y.add_to_counter(iteration); }
	EXPECT_EQ(Y.get_counter(0), 3); // iterations 0, 3 and 6
	// the denominator counted the same three, which is what keeps the fraction a probability
	EXPECT_EQ(Y.get_total_counts(), 3u);
	EXPECT_DOUBLE_EQ(Y.get_fraction_of_ones(0), 1.0);
}

TEST(YStorageDense_Tests, counter_past_the_16_bit_maximum_throws) {
	TStorageYDense Y(TStorageYDense::MAX_COUNTER, {1, 1});
	ASSERT_EQ(Y.get_thinning_factor(), 1u);
	Y.insert_one(0);

	// the chain it was sized for fills the counter exactly
	for (size_t iteration = 0; iteration < TStorageYDense::MAX_COUNTER; ++iteration) {
		Y.add_to_counter(iteration);
	}
	EXPECT_EQ(Y.get_counter(0), TStorageYDense::MAX_COUNTER);
	EXPECT_ANY_THROW(Y.add_to_counter(TStorageYDense::MAX_COUNTER));
}

TEST(YStorageDense_Tests, reset_counts_zeroes_every_counter_and_keeps_the_states) {
	TStorageYDense Y(1000, {1, 5});
	Y.insert_one(0);
	Y.insert_one(2);
	Y.add_to_counter(0);
	Y.add_to_counter(1);
	ASSERT_EQ(Y.get_counter(0), 2);

	Y.reset_counts();
	EXPECT_EQ(Y.get_counter(0), 0);
	EXPECT_EQ(Y.get_counter(2), 0);
	EXPECT_TRUE(Y.is_one(0));
	EXPECT_TRUE(Y.is_one(2));
}

TEST(YStorageDense_Tests, set_state_preserves_the_counter_but_insert_starts_it_over) {
	TStorageYDense Y(1000, {1, 5});
	Y.insert_one(2);
	Y.add_to_counter(0);
	ASSERT_EQ(Y.get_counter(2), 1);

	// flipping a cell in place is a move of the chain: what it was a one for still counts
	Y.set_state(2, false);
	EXPECT_EQ(Y.get_counter(2), 1);
	Y.set_state(2, true);
	EXPECT_EQ(Y.get_counter(2), 1);

	// inserting writes the cell anew, as it does in the sparse field
	Y.insert_one(2);
	EXPECT_EQ(Y.get_counter(2), 0);
}

TEST(YStorageDense_Tests, remove_zeros_drops_the_counter_of_a_cell_that_is_no_longer_one) {
	TStorageYDense Y(1000, {1, 5});
	Y.insert_one(1);
	Y.insert_one(3);
	Y.add_to_counter(0);
	Y.set_state(1, false);

	Y.remove_zeros();
	EXPECT_EQ(Y.get_counter(1), 0); // sparse erases the cell outright; this is the same answer
	EXPECT_EQ(Y.get_counter(3), 1);
	EXPECT_TRUE(Y.is_one(3));
	EXPECT_EQ(Y.total_size_of_container_space(), 5u);
}

TEST(YStorageDense_Tests, get_fraction_of_ones) {
	TStorageYDense Y(1000, {1, 5});
	Y.insert_one(3);
	for (size_t iteration = 0; iteration < 4; ++iteration) { Y.add_to_counter(iteration); }

	// as in TStorageY_Tests: four of the four counted iterations, spelled without dividing the
	// count by itself
	EXPECT_EQ(Y.get_total_counts(), 4u);
	EXPECT_DOUBLE_EQ(Y.get_fraction_of_ones(3), 1.0);
	// a cell that was never a one has fraction 0
	EXPECT_DOUBLE_EQ(Y.get_fraction_of_ones(0), 0.0);
}

TEST(YStorageDense_Tests, fill_current_state_reports_every_cell_as_existing) {
	TStorageYDense Y(1000, {1, 4});
	Y.insert_one(1);

	std::vector<uint8_t> state;
	std::vector<uint8_t> exists;
	std::vector<size_t> linear;
	Y.fill_current_state(IndexArray{0, 0}, /*K=*/4, /*increment=*/1, state, exists, linear);

	EXPECT_EQ(linear, (std::vector<size_t>{0, 1, 2, 3}));
	EXPECT_EQ(state, (std::vector<uint8_t>{0, 1, 0, 0}));
	EXPECT_EQ(exists, (std::vector<uint8_t>{1, 1, 1, 1}));
}
