#include "constants.h"
#include "storages/y_storage/TStorageY.h"
#include "storages/y_storage/TStorageYMatrix.h"
#include "gtest/gtest.h"
#include <cstdint>
#include <vector>

//-----------------------------------
// TStorageY (single 2-byte state + 15-bit counter; no longer stores a linear index)
//-----------------------------------

TEST(YStorage_Tests, flip_state) {
	TStorageY y;
	y.set_state(true);
	EXPECT_TRUE(y.is_one());
	y.switch_state();
	EXPECT_FALSE(y.is_one());
	y.switch_state();
	EXPECT_TRUE(y.is_one());
}

TEST(YStorage_Tests, switch_state_preserves_counter) {
	TStorageY y;
	y.set_state(true);
	y.set_counter(25);
	y.switch_state();
	EXPECT_FALSE(y.is_one());
	EXPECT_EQ(y.get_counter(), 25);
	y.switch_state();
	EXPECT_TRUE(y.is_one());
	EXPECT_EQ(y.get_counter(), 25);
}

TEST(YStorage_Tests, set_counter) {
	TStorageY y;
	y.set_counter(1458);
	EXPECT_EQ(y.get_counter(), 1458);
}

TEST(YStorage_Tests, set_count_zero) {
	TStorageY y;
	y.set_counter(0);
	EXPECT_EQ(y.get_counter(), 0);
}

TEST(YStorage_Tests, set_count_max) {
	TStorageY y;
	y.set_counter(TStorageY::MAX_COUNTER);
	EXPECT_EQ(y.get_counter(), TStorageY::MAX_COUNTER);
}

TEST(YStorage_Tests, set_counter_exceeds_max_throws) {
	TStorageY y;
	EXPECT_ANY_THROW(y.set_counter(TStorageY::MAX_COUNTER + 1));
}

TEST(YStorage_Tests, set_state_true) {
	TStorageY y;
	y.set_state(true);
	EXPECT_TRUE(y.is_one());
}

TEST(YStorage_Tests, set_state_false) {
	TStorageY y;
	y.set_state(false);
	EXPECT_FALSE(y.is_one());
}

TEST(YStorage_Tests, constructor_from_state) {
	TStorageY one{true};
	EXPECT_TRUE(one.is_one());
	EXPECT_EQ(one.get_counter(), 0);
	TStorageY zero{false};
	EXPECT_FALSE(zero.is_one());
	EXPECT_EQ(zero.get_counter(), 0);
}

TEST(YStorage_Tests, update_counter) {
	TStorageY y;
	y.set_state(true);
	y.set_counter(0);
	y.update_counter();
	EXPECT_EQ(y.get_counter(), 1);
}

TEST(YStorage_Tests, update_counter_false) {
	TStorageY y;
	y.set_state(false);
	y.set_counter(0);
	y.update_counter();
	EXPECT_EQ(y.get_counter(), 0);
}

TEST(YStorage_Tests, update_counter_to_two) {
	TStorageY y;
	y.set_counter(1);
	y.set_state(true);
	y.update_counter();
	EXPECT_EQ(y.get_counter(), 2);
}

TEST(YStorage_Tests, update_counter_incremental) {
	TStorageY y;
	y.set_counter(0);
	y.set_state(true);
	y.update_counter();
	EXPECT_EQ(y.get_counter(), 1);
}

TEST(YStorage_Tests, update_counter_incremental_2) {
	TStorageY y;
	y.set_state(true);
	y.update_counter();
	y.set_state(false);
	y.update_counter();
	EXPECT_EQ(y.get_counter(), 1);
}

TEST(YStorage_Tests, update_counter_incremental_3) {
	TStorageY y;
	y.set_state(true);
	y.update_counter();
	y.set_state(false);
	y.update_counter();
	y.set_state(true);
	y.update_counter();
	EXPECT_EQ(y.get_counter(), 2);
}

TEST(YStorage_Tests, test_constructor) {
	constexpr TStorageY y;
	EXPECT_EQ(y.get_counter(), 0);
	EXPECT_FALSE(y.is_one());
	EXPECT_TRUE(y.is_empty());
}

TEST(YStorage_Tests, reset_counter) {
	TStorageY y;
	y.set_counter(1458);
	y.reset_counter();
	EXPECT_EQ(y.get_counter(), 0);
}

TEST(YStorage_Tests, reset_counter_keeps_state) {
	TStorageY y;
	y.set_state(true);
	y.set_counter(1458);
	y.reset_counter();
	EXPECT_EQ(y.get_counter(), 0);
	EXPECT_TRUE(y.is_one());
}

TEST(YStorage_Tests, equality) {
	TStorageY a{true};
	a.set_counter(5);
	TStorageY b{true};
	b.set_counter(5);
	TStorageY c{false};
	c.set_counter(5);
	EXPECT_TRUE(a == b);
	EXPECT_FALSE(a == c);
	EXPECT_TRUE(a != c);
}

TEST(YStorage_Tests, is_empty) {
	TStorageY a;
	EXPECT_TRUE(a.is_empty()); // default: state false, counter 0
	a.set_state(true);
	EXPECT_FALSE(a.is_empty());
	a.set_state(false);
	EXPECT_TRUE(a.is_empty()); // state false + counter 0 == empty (== TStorageY{})
	TStorageY b;
	b.set_counter(1);
	EXPECT_FALSE(b.is_empty()); // counter alone makes it non-empty
}

//-----------------------------------
// TStorageYMatrix (sparse 2-D matrix; indices are linear indices in Y space)
// Tests use a single-row layout {1, N} so the linear index equals the column.
//-----------------------------------

TEST(YStorageMatrix_Tests, remove_zeros_removes_zero_state_elements) {
	TStorageYMatrix Y(1000, {1, 5});
	Y.insert_one(1);
	Y.insert_one(3);
	Y.insert_zero(2);
	EXPECT_EQ(Y.get_stored_entries().size(), 3);
	Y.remove_zeros();
	const auto entries = Y.get_stored_entries();
	ASSERT_EQ(entries.size(), 2);
	EXPECT_EQ(entries[0].first, 1u);
	EXPECT_EQ(entries[1].first, 3u);
	EXPECT_TRUE(Y.is_one(1));
	EXPECT_TRUE(Y.is_one(3));
}

TEST(YStorageMatrix_Tests, remove_zeros_all_zeros_empties_matrix) {
	TStorageYMatrix Y(1000, {1, 5});
	Y.insert_zero(0);
	Y.insert_zero(2);
	Y.insert_zero(4);
	EXPECT_EQ(Y.get_stored_entries().size(), 3);
	Y.remove_zeros();
	EXPECT_EQ(Y.get_stored_entries().size(), 0);
	EXPECT_TRUE(Y.empty());
}

TEST(YStorageMatrix_Tests, remove_zeros_all_ones_unchanged) {
	TStorageYMatrix Y(1000, {1, 5});
	Y.insert_one(0);
	Y.insert_one(2);
	Y.insert_one(4);
	EXPECT_EQ(Y.get_stored_entries().size(), 3);
	Y.remove_zeros();
	EXPECT_EQ(Y.get_stored_entries().size(), 3);
	EXPECT_EQ(Y.number_of_ones(), 3);
}

TEST(YStorageMatrix_Tests, set_to_zero_then_remove) {
	TStorageYMatrix Y(1000, {1, 5});
	Y.insert_one(1);
	Y.insert_one(2);
	Y.insert_one(3);
	Y.set_to_zero(2); // linear index 2
	Y.remove_zeros();
	const auto entries = Y.get_stored_entries();
	ASSERT_EQ(entries.size(), 2);
	EXPECT_EQ(entries[0].first, 1u);
	EXPECT_EQ(entries[1].first, 3u);
}

TEST(YStorageMatrix_Tests, reset_counts_sets_all_counters_to_zero) {
	TStorageYMatrix Y(1000, {1, 5});
	Y.insert_one(0);
	Y.insert_one(2);
	Y.insert_one(4);
	Y.add_to_counter(0);
	Y.add_to_counter(0);
	Y.add_to_counter(0);
	EXPECT_EQ(Y[0].get_counter(), 3);
	EXPECT_EQ(Y[2].get_counter(), 3);
	EXPECT_EQ(Y[4].get_counter(), 3);
	Y.reset_counts();
	EXPECT_EQ(Y[0].get_counter(), 0);
	EXPECT_EQ(Y[2].get_counter(), 0);
	EXPECT_EQ(Y[4].get_counter(), 0);
}

TEST(YStorageMatrix_Tests, reset_counts_does_not_affect_states) {
	TStorageYMatrix Y(1000, {1, 5});
	Y.insert_one(1);
	Y.insert_one(3);
	Y.add_to_counter(0);
	Y.reset_counts();
	EXPECT_TRUE(Y.is_one(1));
	EXPECT_TRUE(Y.is_one(3));
}

TEST(YStorageMatrix_Tests, set_state_flips_in_place) {
	TStorageYMatrix Y(1000, {1, 5});
	Y.insert_one(2);
	Y.add_to_counter(0); // counter of linear-2 cell -> 1
	Y.set_state(2, false);
	EXPECT_FALSE(Y.is_one(2));
	EXPECT_EQ(Y[2].get_counter(), 1); // flipping state preserves the counter
	Y.set_state(2, true);
	EXPECT_TRUE(Y.is_one(2));
	EXPECT_EQ(Y[2].get_counter(), 1);
}

TEST(YStorageMatrix_Tests, get_fraction_of_ones) {
	TStorageYMatrix Y(1000, {1, 5});
	Y.insert_one(3);
	for (size_t it = 0; it < 4; ++it) { Y.add_to_counter(0); } // counter -> 4
	// total_counts == n_iterations / thinning_factor == 1000
	EXPECT_DOUBLE_EQ(Y.get_fraction_of_ones(3), 4.0 / static_cast<double>(Y.get_total_counts()));
	// a cell that was never stored has fraction 0
	EXPECT_DOUBLE_EQ(Y.get_fraction_of_ones(0), 0.0);
}

TEST(YStorageMatrix_Tests, add_data) {
	TStorageYMatrix Y(1000, {2, 3});
	EXPECT_EQ(Y.total_size_of_container_space(), 6u);
	EXPECT_ANY_THROW(Y.insert_one(7));

	Y.insert_one(5);
	EXPECT_EQ(Y.get_full_Y_binary_vector(), (std::vector<uint8_t>{0, 0, 0, 0, 0, 1}));

	Y.set_to_zero(5); // linear index 5 -> (row 1, col 2)
	EXPECT_EQ(Y.get_full_Y_binary_vector(), (std::vector<uint8_t>{0, 0, 0, 0, 0, 0}));

	EXPECT_EQ(Y.get_multi_dimensional_index(5), (IndexArray{1, 2}));
	EXPECT_FALSE(Y.is_one(5));
	EXPECT_EQ(Y.get_thinning_factor(), 1u);

	// linear 5 is still stored (as a zero) until remove_zeros
	const auto entries = Y.get_stored_entries();
	ASSERT_EQ(entries.size(), 1u);
	EXPECT_EQ(entries[0].first, 5u);
}
