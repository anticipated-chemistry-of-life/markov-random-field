//
// Created by VISANI Marco on 17.10.2024.
//

#include "TStorageZ.h"
#include "gtest/gtest.h"

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
