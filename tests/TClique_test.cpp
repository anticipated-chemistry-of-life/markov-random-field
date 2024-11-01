#include "TClique.h"
#include "TStorageYVector.h"
#include "TStorageZVector.h"
#include "TTree.h"
#include "gtest/gtest.h"
#include <vector>

TEST(TCliqueTest, UpdateZTest) {
	// Set up the TTree
	TTree tree(0);
	tree.load_from_file("../tests/test_data/loading_tree.tsv");

	// Set up the dimensions and number of nodes
	std::vector<size_t> start_index = {0};
	size_t variable_dimension       = 0;
	size_t n_nodes                  = tree.size();

	// Create a TClique instance
	TClique clique(start_index, variable_dimension, n_nodes);

	// Set up the TStorageYVector and TStorageZVector
	std::vector<size_t> dimensions = {11};
	TStorageYVector Y(1000, dimensions);
	TStorageZVector Z(dimensions);

	// Set some values in Y and Z
	Y.set_to_one(0);
	Y.set_to_one(1);
	Z.set_to_one(2);
	Z.set_to_one(3);

	// Call the update_Z_test function
	std::vector<bool> result = clique.update_Z_test(Y, Z, tree);

	// Verify the output
	std::vector<bool> expected = {false, false, true,  true,  true, true,
	                              false, false, false, false, false}; // Adjust this based on your expected output
	EXPECT_EQ(result, expected);
}
