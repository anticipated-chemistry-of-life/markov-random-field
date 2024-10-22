#include "TTree.h"
#include "gtest/gtest.h"
#include <cstddef>
#include <string>
#include <vector>

TEST(Tree_test, load_tree) {
	TTree tree;
	tree.load_from_file("../tests/test_data/loading_tree.csv"); // The path looks weird because we actually run the
	                                                            // tests from the build directory
	uint roots = tree.count_roots();
	EXPECT_EQ(roots, 2);
	std::string mammal                 = std::string("mammal");
	const std::vector<size_t> children = tree.get_node(mammal).children();
	EXPECT_EQ(children.size(), 2);
	EXPECT_EQ(tree.get_node(mammal).branchLengthToParent(), 17.9);
}
