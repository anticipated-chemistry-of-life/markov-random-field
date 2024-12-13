#include "TClique.h"
#include "TTree.h"
#include "Types.h"
#include "coretools/Main/TParameters.h"
#include "coretools/devtools.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
using namespace testing;

TEST(Tree_test, load_tree) {
	TTree tree(0);
	tree.load_from_file("../tests/test_data/loading_tree.tsv"); // The path looks weird because we actually run the
	                                                            // tests from the build directory
	uint roots = tree.number_of_roots();
	EXPECT_EQ(roots, 3);
	std::string mammal                 = std::string("mammal");
	const std::vector<size_t> children = tree.get_node(mammal).children_indices_in_tree();
	EXPECT_EQ(children.size(), 2);
	EXPECT_EQ(tree.get_node(mammal).get_branch_length_bin(), 24);
}

TEST(Tree_test, test_discretisation_of_branch_length) {
	coretools::instances::parameters().add("K", 10);
	TTree tree(0);
	tree.load_from_file("../tests/test_data/binning.tsv");
	std::vector<TypeBinBranches> vector = tree.get_all_binned_branch_lengths();
	EXPECT_THAT(vector, ElementsAre(3, 0, 1, 0, 2, 0));
}
