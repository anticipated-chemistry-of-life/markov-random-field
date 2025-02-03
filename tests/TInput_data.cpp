#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <cstdint>
#include <vector>

#include "TLotus.h"
#include "TTree.h"
#include "Types.h"
#include "stattools/ParametersObservations/TParameter.h"
using namespace testing;

struct ModelDummy {
	PriorOnGamma prior_on_gamma{};
	stattools::TParameter<SpecGamma, TLotus<false>> gamma{&prior_on_gamma};
};

TEST(Tinput, test_reading_links) {

	std::vector<TTree> trees;
	trees.emplace_back(1, "../tests/test_data/molecules.tsv", "molecules");
	trees.emplace_back(0, "../tests/test_data/loading_tree.tsv", "species");

	TStorageYVector Y;
	ModelDummy model;
	TLotus links(trees, &model.gamma, Y, "simulations");

	links.load_from_file("../tests/test_data/links.tsv");

	EXPECT_EQ(links.get_Lotus().size(), 6);

	std::vector<uint64_t> vector;
	for (size_t i = 0; i < links.get_Lotus().size(); ++i) {
		vector.push_back(links.get_Lotus()[i].get_linear_index_in_container_space());
	}

	EXPECT_THAT(vector, ElementsAre(2, 6, 7, 9, 15, 16));
	stattools::instances::dagBuilder().clear();
}
