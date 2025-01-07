#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <cstdint>
#include <vector>

#include "TLotus.h"
#include "TTree.h"
using namespace testing;

TEST(Tinput, test_reading_links) {
	TTree tree_1 = TTree(0);
	TTree tree_2 = TTree(1);

	tree_1.load_from_file("../tests/test_data/loading_tree.tsv", "species");
	tree_2.load_from_file("../tests/test_data/molecules.tsv", "molecules");

	std::vector<TTree> trees = {tree_1, tree_2};

	TLotus links(trees);

	links.load_from_file("../tests/test_data/links.tsv");

	EXPECT_EQ(links.get_TStorageYVector().size(), 6);

	std::vector<uint64_t> vector;
	for (size_t i = 0; i < links.get_TStorageYVector().size(); ++i) {
		vector.push_back(links.get_TStorageYVector()[i].get_linear_index_in_container_space());
	}

	EXPECT_THAT(vector, ElementsAre(1, 4, 6, 10, 11, 14));
}
