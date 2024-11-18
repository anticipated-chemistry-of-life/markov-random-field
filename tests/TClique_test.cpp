#include "TClique.h"
#include "TStorageYVector.h"
#include "TStorageZVector.h"
#include "TTree.h"
#include "gtest/gtest.h"
#include <vector>

// TEST(TCliqueTest, get_z) {
// 	// Set up the TTree
// 	TTree tree(0);
// 	TTree tree_2(1);
// 	tree.load_from_file("../tests/test_data/loading_tree.tsv");
// 	tree_2.load_from_file("../tests/test_data/loading_tree.tsv");

// 	// Set up the TStorageYVector and TStorageZVector
// 	std::vector<size_t> dimensions = {11, 11};
// 	TStorageYVector Y(1000, {6, 6});
// 	TStorageZVector Z({5, 5});

// 	// Set some values in Y and Z
// 	Y.set_to_one(0);
// 	Y.set_to_one(1);
// 	Y.set_to_one(8);
// 	Y.set_to_one(10);
// 	Z.set_to_one(2);
// 	Z.set_to_one(3);
// 	Z.set_to_one(6);
// 	Z.set_to_one(8);
// 	Z.set_to_one(9);

// 	// Call the update_Z_test function
// 	// std::vector<bool> result = clique.update_Z_test(Y, Z, tree);
// 	tree_2.initialize_cliques({tree, tree_2});
// 	auto &cliques = tree_2.get_cliques();

// 	// Verify the output
// 	std::vector<bool> expected = {false, false, true, true, true, true, false, false, false, false, false};
// 	EXPECT_EQ(cliques[0].fill_current_state(Y, Z, tree), expected);

// 	// Verify the output
// 	expected = {false, true, false, false, false, true, true, false, true, true, false};
// 	EXPECT_EQ(cliques[1].fill_current_state(Y, Z, tree), expected);
// }
