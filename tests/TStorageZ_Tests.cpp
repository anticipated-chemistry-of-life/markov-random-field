//
// Created by VISANI Marco on 17.10.2024.
//

#include "constants.h"
#include "storages/z_storage/TStorageZ.h"
#include "storages/z_storage/TStorageZMatrix.h"
#include "gtest/gtest.h"
#include <cstddef>
#include <cstdint>
#include <vector>

//-----------------------------------
// TStorageZ (single state byte; since Z migrated to a TSparseMatrix the cell no
// longer stores its own linear index, only the state bit)
//-----------------------------------

TEST(ZStorage_Tests, default_constructor_is_empty) {
	const TStorageZ z;
	EXPECT_FALSE(z.is_one());
	EXPECT_TRUE(z.is_empty());
}

TEST(ZStorage_Tests, constructor_from_state) {
	const TStorageZ one{true};
	EXPECT_TRUE(one.is_one());
	EXPECT_FALSE(one.is_empty());
	const TStorageZ zero{false};
	EXPECT_FALSE(zero.is_one());
	EXPECT_TRUE(zero.is_empty());
}

TEST(ZStorage_Tests, set_state) {
	TStorageZ z;
	z.set_state(true);
	EXPECT_TRUE(z.is_one());
	z.set_state(false);
	EXPECT_FALSE(z.is_one());
	// idempotent: setting the same state twice keeps it
	z.set_state(false);
	EXPECT_FALSE(z.is_one());
	z.set_state(true);
	z.set_state(true);
	EXPECT_TRUE(z.is_one());
}

TEST(ZStorage_Tests, switch_state) {
	TStorageZ z{true};
	EXPECT_TRUE(z.is_one());
	z.switch_state();
	EXPECT_FALSE(z.is_one());
	z.switch_state();
	EXPECT_TRUE(z.is_one());
}

TEST(ZStorage_Tests, is_empty_tracks_state) {
	TStorageZ z;
	EXPECT_TRUE(z.is_empty()); // default: state false
	z.set_state(true);
	EXPECT_FALSE(z.is_empty());
	z.set_state(false);
	EXPECT_TRUE(z.is_empty()); // state false == empty == TStorageZ{}
}

TEST(ZStorage_Tests, equality) {
	const TStorageZ a{true};
	const TStorageZ b{true};
	const TStorageZ c{false};
	EXPECT_TRUE(a == b);
	EXPECT_FALSE(a == c);
	EXPECT_TRUE(a != c);
}

TEST(ZStorage_Tests, size_is_one_byte) { EXPECT_EQ(sizeof(TStorageZ), 1u); }

//-----------------------------------
// TStorageZMatrix (sparse 2-D matrix; indices are linear indices in Z space)
// A single-row layout {1, N} makes the linear index equal the column.
//-----------------------------------

TEST(ZStorageMatrix_Tests, is_one_missing_cell_reads_false) {
	TStorageZMatrix Z({1, 5});
	EXPECT_FALSE(Z.is_one(0));
	EXPECT_FALSE(Z.is_one(4));
	EXPECT_TRUE(Z.empty());
	EXPECT_EQ(Z.size(), 0u);
}

TEST(ZStorageMatrix_Tests, insert_one_and_lookup) {
	TStorageZMatrix Z({1, 5});
	Z.insert_one(2);
	EXPECT_TRUE(Z.is_one(2));
	EXPECT_FALSE(Z.is_one(0));
	EXPECT_FALSE(Z.empty());
	EXPECT_EQ(Z.size(), 1u);
}

TEST(ZStorageMatrix_Tests, insert_one_by_multidim_index) {
	TStorageZMatrix Z({2, 3});
	const IndexArray idx{1, 2}; // row 1, col 2 -> linear 5
	EXPECT_EQ(Z.get_linear_index_in_Z_space(idx), 5u);
	Z.insert_one(idx);
	EXPECT_TRUE(Z.is_one(5));
	EXPECT_TRUE(Z.is_one(Z.get_linear_index_in_Z_space(idx)));
}

TEST(ZStorageMatrix_Tests, set_state_flips_in_place) {
	TStorageZMatrix Z({1, 5});
	Z.insert_one(2);
	Z.set_state(2, false);
	EXPECT_FALSE(Z.is_one(2));
	Z.set_state(2, true);
	EXPECT_TRUE(Z.is_one(2));
}

TEST(ZStorageMatrix_Tests, set_state_inserts_missing_cell) {
	TStorageZMatrix Z({1, 5});
	Z.set_state(3, true); // cell did not exist yet
	EXPECT_TRUE(Z.is_one(3));
	EXPECT_EQ(Z.size(), 1u);
}

TEST(ZStorageMatrix_Tests, remove_zeros_removes_zero_state_elements) {
	TStorageZMatrix Z({1, 5});
	Z.insert_one(1);
	Z.insert_one(3);
	Z.insert_zero(2);
	EXPECT_EQ(Z.get_stored_entries().size(), 3u);
	Z.remove_zeros();
	const auto entries = Z.get_stored_entries();
	ASSERT_EQ(entries.size(), 2u);
	EXPECT_EQ(entries[0].first, 1u);
	EXPECT_EQ(entries[1].first, 3u);
	EXPECT_TRUE(Z.is_one(1));
	EXPECT_TRUE(Z.is_one(3));
}

TEST(ZStorageMatrix_Tests, remove_zeros_all_zeros_empties_matrix) {
	TStorageZMatrix Z({1, 5});
	Z.insert_zero(0);
	Z.insert_zero(2);
	Z.insert_zero(4);
	EXPECT_EQ(Z.get_stored_entries().size(), 3u);
	Z.remove_zeros();
	EXPECT_EQ(Z.get_stored_entries().size(), 0u);
	EXPECT_TRUE(Z.empty());
}

TEST(ZStorageMatrix_Tests, remove_zeros_all_ones_unchanged) {
	TStorageZMatrix Z({1, 5});
	Z.insert_one(0);
	Z.insert_one(2);
	Z.insert_one(4);
	Z.remove_zeros();
	EXPECT_EQ(Z.size(), 3u);
	EXPECT_TRUE(Z.is_one(0));
	EXPECT_TRUE(Z.is_one(2));
	EXPECT_TRUE(Z.is_one(4));
}

TEST(ZStorageMatrix_Tests, get_stored_entries_ascending_linear_index) {
	TStorageZMatrix Z({2, 3}); // 2 rows, 3 cols (row-major linear index)
	Z.insert_one(5);           // (1, 2)
	Z.insert_one(0);           // (0, 0)
	Z.insert_one(4);           // (1, 1)
	Z.insert_one(1);           // (0, 1)
	const auto entries = Z.get_stored_entries();
	ASSERT_EQ(entries.size(), 4u);
	EXPECT_EQ(entries[0].first, 0u);
	EXPECT_EQ(entries[1].first, 1u);
	EXPECT_EQ(entries[2].first, 4u);
	EXPECT_EQ(entries[3].first, 5u);
}

TEST(ZStorageMatrix_Tests, add_data) {
	TStorageZMatrix Z({2, 3});
	EXPECT_EQ(Z.total_size_of_container_space(), 6u);
	EXPECT_ANY_THROW(Z.insert_one(7)); // linear index past the container size

	Z.insert_one(5);
	EXPECT_EQ(Z.get_full_Z_binary_vector(), (std::vector<size_t>{0, 0, 0, 0, 0, 1}));

	Z.set_state(5, false); // linear index 5 -> (row 1, col 2)
	EXPECT_EQ(Z.get_full_Z_binary_vector(), (std::vector<size_t>{0, 0, 0, 0, 0, 0}));

	EXPECT_EQ(Z.get_multi_dimensional_index(5), (IndexArray{1, 2}));
	EXPECT_FALSE(Z.is_one(5));

	// linear 5 is still stored (as a zero) until remove_zeros
	const auto entries = Z.get_stored_entries();
	ASSERT_EQ(entries.size(), 1u);
	EXPECT_EQ(entries[0].first, 5u);
	EXPECT_FALSE(entries[0].second.is_one());
}

TEST(ZStorageMatrix_Tests, insert_one_exceeds_size_throws) {
	TStorageZMatrix Z({2, 3}); // total size 6
	EXPECT_ANY_THROW(Z.insert_one(6));
	EXPECT_NO_THROW(Z.insert_one(5));
}

TEST(ZStorageMatrix_Tests, insert_zero_exceeds_size_throws) {
	TStorageZMatrix Z({2, 3}); // total size 6
	EXPECT_ANY_THROW(Z.insert_zero(6));
	EXPECT_NO_THROW(Z.insert_zero(5));
}

// A run of cells is a start, a count and a stride, and the caller does that arithmetic. A point
// lookup in a sorted-vector matrix costs a search, and answering one per cell is what this backend
// pays for a run; what it owes is the same answers a dense array gives.

// stride == 1: the varying dimension is the last one, so the cells are consecutive.
TEST(ZStorageMatrix_Tests, a_run_along_the_last_dimension_reads_consecutive_cells) {
	TStorageZMatrix Z({1, 6}); // 1 row, 6 cols
	Z.insert_one(1);
	Z.insert_one(3);
	Z.insert_zero(2); // stored but state zero

	// a stored zero and an absent cell both read as zero
	const std::vector<uint8_t> expected{0, 1, 0, 1, 0, 0};
	for (size_t k = 0; k < expected.size(); ++k) {
		EXPECT_EQ(Z.get_multi_dimensional_index(k), (IndexArray{0, k})) << "cell " << k;
		EXPECT_EQ(Z.is_one(k), expected[k] != 0) << "cell " << k;
	}
}

// stride > 1: the varying dimension is the first one, so the cells are one row apart.
TEST(ZStorageMatrix_Tests, a_run_along_the_first_dimension_steps_by_a_row) {
	TStorageZMatrix Z({3, 2}); // 3 rows, 2 cols -> a row is 2 cells wide
	Z.insert_one(2);           // (row 1, col 0) linear 2
	Z.insert_one(4);           // (row 2, col 0) linear 4
	Z.insert_one(3);           // (row 1, col 1) linear 3 -> a cell of the other column

	const std::vector<uint8_t> expected{0, 1, 1};
	for (size_t k = 0; k < expected.size(); ++k) {
		const size_t linear = k * 2;
		EXPECT_EQ(Z.get_multi_dimensional_index(linear), (IndexArray{k, 0})) << "cell " << k;
		EXPECT_EQ(Z.is_one(linear), expected[k] != 0) << "cell " << k;
	}
}

TEST(ZStorageMatrix_Tests, a_run_that_starts_partway_along_a_row_reads_from_there) {
	TStorageZMatrix Z({1, 6});
	Z.insert_one(2);
	Z.insert_one(4);
	Z.insert_one(5); // past the end of the run

	const std::vector<uint8_t> expected{1, 0, 1};
	for (size_t k = 0; k < expected.size(); ++k) {
		EXPECT_EQ(Z.is_one(2 + k), expected[k] != 0) << "cell " << k;
	}
}

TEST(ZStorageMatrix_Tests, insert_in_Z_bulk_merges_and_sorts) {
	TStorageZMatrix Z({1, 6});
	const std::vector<std::vector<size_t>> to_insert = {{1, 4}, {2}};
	Z.insert_in_Z(to_insert);
	EXPECT_TRUE(Z.is_one(1));
	EXPECT_TRUE(Z.is_one(2));
	EXPECT_TRUE(Z.is_one(4));
	EXPECT_FALSE(Z.is_one(0));
	const auto entries = Z.get_stored_entries();
	ASSERT_EQ(entries.size(), 3u);
	EXPECT_EQ(entries[0].first, 1u);
	EXPECT_EQ(entries[1].first, 2u);
	EXPECT_EQ(entries[2].first, 4u);
}

TEST(ZStorageMatrix_Tests, get_full_Z_binary_vector_multi_row) {
	TStorageZMatrix Z({2, 3});
	Z.insert_one(0); // (0, 0)
	Z.insert_one(4); // (1, 1)
	EXPECT_EQ(Z.get_full_Z_binary_vector(), (std::vector<size_t>{1, 0, 0, 0, 1, 0}));
}
