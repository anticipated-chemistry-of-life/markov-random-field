#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <cstdint>
#include <vector>

#include "TLotus.h"
#include "TTree.h"
#include "stattools/ParametersObservations/TParameter.h"
using namespace testing;

TEST(Tinput, test_reading_links) {
	TTree tree_1 = TTree(0, "../tests/test_data/loading_tree.tsv", "species");
	TTree tree_2 = TTree(1, "../tests/test_data/molecules.tsv", "molecules");

	std::vector<TTree> trees = {tree_1, tree_2};

	TLotus links(trees);

	links.load_from_file("../tests/test_data/links.tsv");

	EXPECT_EQ(links.get_Lotus().size(), 6);

	std::vector<uint64_t> vector;
	for (size_t i = 0; i < links.get_Lotus().size(); ++i) {
		vector.push_back(links.get_Lotus()[i].get_linear_index_in_container_space());
	}

	EXPECT_THAT(vector, ElementsAre(1, 4, 6, 10, 11, 14));
}
