#include "TStorageY.h"
#include "TStorageYSet.h"
#include "gtest/gtest.h"
#include <cmath>

TEST(YStorageYSet_Tests, insertions) {
	TStorageYSet y;
	y.initialize(100, {1000, 2000});
	y.insert_one(100);
	y.insert_one(105);
	y.insert_zero(203);
	EXPECT_EQ(y.is_one(100), true);
	EXPECT_EQ(y.is_one(106), false);
	EXPECT_EQ(y.is_one(203), false);
}
