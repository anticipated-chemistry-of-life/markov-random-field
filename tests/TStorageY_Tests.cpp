#include "storages/y_storage/TStorageY.h"
#include "storages/y_storage/TStorageYMatrix.h"
#include "gtest/gtest.h"
#include <stdexcept>
#include <vector>

TEST(YStorage_Tests, flip_state) {
	TStorageY y;
	y.set_state(true);
	y.set_linear_index_in_Y_space(25);
	EXPECT_EQ(y.is_one(), true);
	EXPECT_EQ(y.get_linear_index_in_container_space(), 25);
	y.switch_state();

	EXPECT_EQ(y.is_one(), false);
	EXPECT_EQ(y.get_linear_index_in_container_space(), 25);
	y.switch_state();
	EXPECT_EQ(y.is_one(), true);
	EXPECT_EQ(y.get_linear_index_in_container_space(), 25);
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
	y.set_counter(std::numeric_limits<uint16_t>::max());
	EXPECT_EQ(y.get_counter(), std::numeric_limits<uint16_t>::max());
}

TEST(YStorage_Tests, set_state_true) {
	TStorageY y;
	y.set_state(true);
	EXPECT_EQ(y.is_one(), true);
}

TEST(YStorage_Tests, set_state_false) {
	TStorageY y;
	y.set_state(false);
	EXPECT_EQ(y.is_one(), false);
}

TEST(YStorage_Tests, set_coordinate) {
	TStorageY y;
	y.set_linear_index_in_Y_space(18987);
	EXPECT_EQ(y.get_linear_index_in_Y_space(), 18987);
}

TEST(YStorage_Tests, set_coordinate_zero) {
	TStorageY y;
	y.set_linear_index_in_Y_space(0);
	EXPECT_EQ(y.get_linear_index_in_Y_space(), 0);
}

TEST(YStorage_Tests, set_coordinate_max) {
	TStorageY y;
	y.set_linear_index_in_Y_space(MAX_LINEAR_INDEX_IN_Y_SPACE);
	EXPECT_EQ(y.get_linear_index_in_Y_space(), MAX_LINEAR_INDEX_IN_Y_SPACE);
}

TEST(YStorage_Tests, constructor_at_max_index_no_throw) {
	EXPECT_NO_THROW(TStorageY{MAX_LINEAR_INDEX_IN_Y_SPACE});
}

TEST(YStorage_Tests, constructor_exceeds_max_index_throws) {
	EXPECT_ANY_THROW(TStorageY{MAX_LINEAR_INDEX_IN_Y_SPACE + 1});
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
	EXPECT_EQ(y.is_one(), false);
	EXPECT_EQ(y.get_linear_index_in_Y_space(), 0);
}

TEST(YStorage_Tests, reset_counter) {
	TStorageY y;
	y.set_counter(1458);
	y.reset_counter();
	EXPECT_EQ(y.get_counter(), 0);
}

TEST(YStorageVector_Tests, remove_zeros_removes_zero_state_elements) {
	TStorageYMatrix Y(1000, {5});
	Y.insert_one(1);
	Y.insert_one(3);
	Y.insert_zero(2);
	EXPECT_EQ(Y.size(), 3);
	Y.remove_zeros();
	EXPECT_EQ(Y.size(), 2);
	EXPECT_EQ(Y[0].get_linear_index_in_container_space(), 1);
	EXPECT_EQ(Y[1].get_linear_index_in_container_space(), 3);
}

TEST(YStorageVector_Tests, remove_zeros_all_zeros_empties_vector) {
	TStorageYMatrix Y(1000, {5});
	Y.insert_zero(0);
	Y.insert_zero(2);
	Y.insert_zero(4);
	EXPECT_EQ(Y.size(), 3);
	Y.remove_zeros();
	EXPECT_EQ(Y.size(), 0);
}

TEST(YStorageVector_Tests, remove_zeros_all_ones_unchanged) {
	TStorageYMatrix Y(1000, {5});
	Y.insert_one(0);
	Y.insert_one(2);
	Y.insert_one(4);
	EXPECT_EQ(Y.size(), 3);
	Y.remove_zeros();
	EXPECT_EQ(Y.size(), 3);
}

TEST(YStorageVector_Tests, remove_zeros_set_to_zero_then_remove) {
	TStorageYMatrix Y(1000, {5});
	Y.insert_one(1);
	Y.insert_one(2);
	Y.insert_one(3);
	Y.set_to_zero(1); // sets element at vector position 1 (linear index 2) to zero
	Y.remove_zeros();
	EXPECT_EQ(Y.size(), 2);
	EXPECT_EQ(Y[0].get_linear_index_in_container_space(), 1);
	EXPECT_EQ(Y[1].get_linear_index_in_container_space(), 3);
}

TEST(YStorageVector_Tests, reset_counts_sets_all_counters_to_zero) {
	TStorageYMatrix Y(1000, {5});
	Y.insert_one(0);
	Y.insert_one(2);
	Y.insert_one(4);
	Y.add_to_counter(0);
	Y.add_to_counter(0);
	Y.add_to_counter(0);
	EXPECT_EQ(Y[0].get_counter(), 3);
	EXPECT_EQ(Y[1].get_counter(), 3);
	EXPECT_EQ(Y[2].get_counter(), 3);
	Y.reset_counts();
	EXPECT_EQ(Y[0].get_counter(), 0);
	EXPECT_EQ(Y[1].get_counter(), 0);
	EXPECT_EQ(Y[2].get_counter(), 0);
}

TEST(YStorageVector_Tests, reset_counts_does_not_affect_states) {
	TStorageYMatrix Y(1000, {5});
	Y.insert_one(1);
	Y.insert_one(3);
	Y.add_to_counter(0);
	Y.reset_counts();
	EXPECT_EQ(Y[0].is_one(), true);
	EXPECT_EQ(Y[1].is_one(), true);
}

TEST(YStorageVector_Tests, add_data) {
	TStorageYMatrix Y(1000, {2, 3});
	const auto &dimensions = Y.total_size_of_container_space();
	EXPECT_EQ(dimensions, 6);
	EXPECT_ANY_THROW(Y.insert_one(7));
	Y.insert_one(5);
	auto binary_vec          = Y.get_full_Y_binary_vector();
	std::vector<int> compare = {0, 0, 0, 0, 0, 1};
	EXPECT_EQ(binary_vec, compare);

	Y.set_to_zero(0);
	binary_vec = Y.get_full_Y_binary_vector();
	compare    = {0, 0, 0, 0, 0, 0};
	EXPECT_EQ(binary_vec, compare);

	auto multi_dim_index                     = Y.get_multi_dimensional_index(5);
	std::vector<size_t> multi_dim_index_true = {1, 2};
	EXPECT_EQ(multi_dim_index, multi_dim_index_true);
	EXPECT_EQ(Y.size(), 1);
	EXPECT_EQ(Y[0].get_linear_index_in_container_space(), 5);
	EXPECT_EQ(Y.get_thinning_factor(), 1);
	TStorageY x;
	EXPECT_THROW(x = Y.at(7), std::out_of_range);
}
