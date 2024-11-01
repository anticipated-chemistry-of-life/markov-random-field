//
// Created by VISANI Marco on 18.10.2024.
//

#include "TStorageZVector.h"
#include "coretools/devtools.h"
#include "gtest/gtest.h"

TEST(ZStorageVector_Tests, test_constructor) {
	const TStorageZVector z;
	EXPECT_EQ(z.is_one(1458), false);
}

TEST(ZStorageVector_Tests, set_to_one) {
	TStorageZVector z;
	z.set_to_one(1458);
	EXPECT_EQ(z.is_one(1458), true);
}

TEST(ZStorageVector_Tests, set_to_zero) {
	TStorageZVector z;
	z.set_to_zero(-1458);
	EXPECT_EQ(z.is_one(1458), false);
}

TEST(ZStorageVector_Tests, set_to_zero_2) {
	TStorageZVector z;
	z.set_to_one(1458);
	z.set_to_zero(1458);
	EXPECT_EQ(z.is_one(1458), false);
}

TEST(ZStorageVector_Tests, add_multiple_entries) {
	TStorageZVector z;
	z.set_to_one(1458);
	z.set_to_one(1459);
	z.set_to_one(1460);
	EXPECT_EQ(z.is_one(1458), true);
	EXPECT_EQ(z.is_one(1459), true);
	EXPECT_EQ(z.is_one(1460), true);
	z.set_to_zero(0);
	EXPECT_EQ(z.is_one(0), false);
}

TEST(ZStorageVector_Tests, min_int32_value) {
	TStorageZVector z;

	z.set_to_zero(std::numeric_limits<int32_t>::min());
	EXPECT_EQ(z.is_one(std::numeric_limits<int32_t>::min()), false);
}

TEST(ZStorageVector_Tests, max_int32_value) {
	TStorageZVector z;
	z.set_to_one(std::numeric_limits<int32_t>::max());
	EXPECT_EQ(z.is_one(std::numeric_limits<int32_t>::max()), true);
}

TEST(ZStorageVector_Tests, remove_zeros) {
	TStorageZVector z;
	z.set_to_one(1458);
	z.set_to_one(1459);
	z.set_to_one(1460);
	z.set_to_zero(1459);
	z.remove_zeros();

	EXPECT_EQ(z.size(), 2);
}
