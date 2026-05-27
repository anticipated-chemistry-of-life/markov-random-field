//
// Created by VISANI Marco on 17.10.2024.
//

#include "storages/z_storage/TStorageZ.h"
#include "storages/z_storage/TStorageZVector.h"
#include "gtest/gtest.h"
#include <cstdint>
#include <limits>

TEST(TStorageZ_Tests, test_constructor) {
	const TStorageZ z(1000);
	EXPECT_EQ(z.is_one(), true);
	EXPECT_EQ(z.get_linear_index_in_Z_space(), 1000);
}

TEST(TStorageZ_Tests, flip_state) {
	TStorageZ z(1000);
	z.set_state(true);

	EXPECT_EQ(z.is_one(), true);

	z.switch_state();
	EXPECT_EQ(z.is_one(), false);
	z.switch_state();
	EXPECT_EQ(z.is_one(), true);
}

TEST(TStorageZ_Tests, test_set_coordinate) {
	TStorageZ z(1000);
	z.set_linear_index_in_Z_space(2000);
	EXPECT_EQ(z.get_linear_index_in_Z_space(), 2000);
}

TEST(TStorageZ_Tests, constructor_at_max_index_no_throw) {
	EXPECT_NO_THROW(TStorageZ{std::numeric_limits<int32_t>::max()});
}

TEST(TStorageZVector_Tests, insert_one_at_max_index_no_throw) {
	TStorageZVector Z({1});
	EXPECT_NO_THROW(Z.insert_one(static_cast<uint32_t>(std::numeric_limits<int32_t>::max())));
}

TEST(TStorageZVector_Tests, insert_one_exceeds_max_index_throws) {
	TStorageZVector Z({1});
	const uint32_t over_max = static_cast<uint32_t>(std::numeric_limits<int32_t>::max()) + 1;
	EXPECT_ANY_THROW(Z.insert_one(over_max));
}

TEST(TStorageZVector_Tests, insert_zero_exceeds_max_index_throws) {
	TStorageZVector Z({1});
	const uint32_t over_max = static_cast<uint32_t>(std::numeric_limits<int32_t>::max()) + 1;
	EXPECT_ANY_THROW(Z.insert_zero(over_max));
}

TEST(TStorageZ_Tests, test_set_state) {
	TStorageZ z(1000);
	z.set_state(false);
	EXPECT_EQ(z.is_one(), false);
	z.set_state(true);
	EXPECT_EQ(z.is_one(), true);

	z.set_state(false);
	EXPECT_EQ(z.is_one(), false);
	z.set_state(false);
	EXPECT_EQ(z.is_one(), false);
	z.set_state(true);
	EXPECT_EQ(z.is_one(), true);
	z.set_state(true);
	EXPECT_EQ(z.is_one(), true);
}
