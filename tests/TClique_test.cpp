#include "TClique.h"
#include "TStorageYVector.h"
#include "TStorageZVector.h"
#include "TTree.h"
#include "coretools/devtools.h"
#include "gtest/gtest.h"
#include <vector>

TEST(TCliqueTest, get_z) {
	// Set up the TTree
	TTree tree(0, "../tests/test_data/loading_tree.tsv", "species_1");
	TTree tree_2(1, "../tests/test_data/loading_tree.tsv", "species_2");

	// Set up the TStorageYVector and TStorageZVector
	TStorageYVector Y(1000, {6, 6});
	TStorageZVector Z({5, 5});

	// Set some values in Y and Z
	Y.insert_one(0);
	Y.insert_one(1);
	Z.insert_one(0);
	Z.insert_one(2);
	Z.insert_one(3);
	Z.insert_one(5);
	Z.insert_one(7);

	tree_2.initialize_cliques_and_Z({tree, tree_2});
	auto &cliques = tree_2.get_cliques();
	cliques[0].set_mus(0.3, 0.2);
	cliques[0].set_lambda();
	cliques[0].initialize(tree.get_a(), tree.get_delta(), tree.get_number_of_bins());
	cliques[0].update_Z(Y, Z, tree_2);
}
