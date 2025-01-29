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
	stattools::TParameter<SpecGamma, TLotus> gamma{&prior_on_gamma};
};

TEST(Tinput, test_reading_links_should_fail) {

	std::vector<TTree> trees;
	trees.emplace_back(0, "../tests/test_data/loading_tree.tsv", "species");
	trees.emplace_back(1, "../tests/test_data/molecules.tsv", "molecules");

	TStorageYVector Y;
	ModelDummy model;
	TLotus links(trees, &model.gamma, Y);

	EXPECT_ANY_THROW(links.load_from_file("../tests/test_data/links.tsv"));
}
