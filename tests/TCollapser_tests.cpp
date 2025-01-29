#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <cstddef>
#include <cstdint>
#include <vector>

#include "TClique.h"
#include "TCollapser.h"
#include "TTree.h"
#include "Types.h"
#include "coretools/devtools.h"
#include "stattools/ParametersObservations/TParameter.h"
using namespace testing;

TEST(Collapser, keep_1_dim) {
	std::vector<TTree> trees;

	trees.emplace_back(0, "../tests/test_data/molecules.tsv", "molecules");
	trees.emplace_back(1, "../tests/test_data/loading_tree.tsv", "tissues");
	trees.emplace_back(2, "../tests/test_data/loading_tree.tsv", "species");

	TCollapser collapser(trees);
	std::vector<std::string> dimension_names_to_keep = {"species"};

	collapser.initialize(dimension_names_to_keep, "lotus");
}

TEST(Collapser, keep_2_dim) {
	std::vector<TTree> trees;

	trees.emplace_back(0, "../tests/test_data/molecules.tsv", "molecules");
	trees.emplace_back(1, "../tests/test_data/loading_tree.tsv", "tissues");
	trees.emplace_back(2, "../tests/test_data/loading_tree.tsv", "species");

	TCollapser collapser(trees);
	std::vector<std::string> dimension_names_to_keep = {"molecules", "species"};

	collapser.initialize(dimension_names_to_keep, "lotus");
}

TEST(Collapser, counter) {
	std::vector<TTree> trees;

	trees.emplace_back(0, "../tests/test_data/molecules.tsv", "molecules");
	trees.emplace_back(1, "../tests/test_data/loading_tree.tsv", "tissues");
	trees.emplace_back(2, "../tests/test_data/loading_tree.tsv", "species");

	for (auto &tree : trees) { tree.initialize_cliques_and_Z(trees); }

	trees[1].get_clique({0, 0, 0}).update_counter_leaves_state_1(2);

	TClique *clique = &trees[1].get_clique({0, 0, 0});
	ASSERT_EQ(clique->get_counter_leaves_state_1(), 2);

	TCollapser collapser(trees);
	std::vector<std::string> dimension_names_to_keep = {"molecules", "species"};

	collapser.initialize(dimension_names_to_keep, "lotus");

	EXPECT_EQ(collapser.x_is_one({0, 0, 0}), true);
	EXPECT_EQ(collapser.x_is_one({0, 1, 0}), true);
	EXPECT_EQ(collapser.x_is_one({0, 0, 1}), false);
}

TEST(Collapser, check_if_other_dims_are_not_collapsed) {
	/// this test checks that if the collapser collapses a specific dimension, then
	/// checking if x_is_one in other dimensions should return false
	std::vector<TTree> trees;

	trees.emplace_back(0, "../tests/test_data/molecules.tsv", "molecules");
	trees.emplace_back(1, "../tests/test_data/loading_tree.tsv", "tissues");
	trees.emplace_back(2, "../tests/test_data/loading_tree.tsv", "species");

	for (auto &tree : trees) { tree.initialize_cliques_and_Z(trees); }

	trees[0].get_clique({0, 0, 0}).update_counter_leaves_state_1(2);

	TCollapser collapser(trees);
	std::vector<std::string> dimension_names_to_keep = {"molecules", "species"};

	collapser.initialize(dimension_names_to_keep, "lotus");

	for (size_t i = 0; i < trees[0].get_number_of_leaves(); ++i) {
		for (size_t j = 0; j < trees[1].get_number_of_leaves(); ++j) {
			for (size_t k = 0; k < trees[2].get_number_of_leaves(); ++k) {
				EXPECT_EQ(collapser.x_is_one({i, j, k}), false);
			}
		}
	}
}

TEST(Collapser, old_state) {
	std::vector<TTree> trees;

	trees.emplace_back(0, "../tests/test_data/molecules.tsv", "molecules");
	trees.emplace_back(1, "../tests/test_data/loading_tree.tsv", "tissues");
	trees.emplace_back(2, "../tests/test_data/loading_tree.tsv", "species");

	for (auto &tree : trees) { tree.initialize_cliques_and_Z(trees); }
	trees[1].get_clique({0, 0, 0}).update_counter_leaves_state_1(1);

	TCollapser collapser(trees);
	std::vector<std::string> dimension_names_to_keep = {"molecules", "species"};
	collapser.initialize(dimension_names_to_keep, "lotus");

	EXPECT_EQ(collapser.x_is_one({0, 1, 0}, true), false);
	EXPECT_EQ(collapser.x_is_one({0, 1, 0}, false), true);
}

TEST(Collapser, old_state_is_0) {
	std::vector<TTree> trees;

	trees.emplace_back(0, "../tests/test_data/molecules.tsv", "molecules");
	trees.emplace_back(1, "../tests/test_data/loading_tree.tsv", "tissues");
	trees.emplace_back(2, "../tests/test_data/loading_tree.tsv", "species");

	for (auto &tree : trees) { tree.initialize_cliques_and_Z(trees); }
	trees[1].get_clique({0, 0, 0}).update_counter_leaves_state_1(0);

	TCollapser collapser(trees);
	std::vector<std::string> dimension_names_to_keep = {"molecules", "species"};
	collapser.initialize(dimension_names_to_keep, "lotus");

	EXPECT_EQ(collapser.x_is_one({0, 1, 0}, true), false);
	EXPECT_EQ(collapser.x_is_one({0, 1, 0}, false), false);
}

TEST(Collapser, old_state_is_bigger_than_1) {
	std::vector<TTree> trees;

	trees.emplace_back(0, "../tests/test_data/molecules.tsv", "molecules");
	trees.emplace_back(1, "../tests/test_data/loading_tree.tsv", "tissues");
	trees.emplace_back(2, "../tests/test_data/loading_tree.tsv", "species");

	for (auto &tree : trees) { tree.initialize_cliques_and_Z(trees); }
	trees[1].get_clique({0, 0, 0}).update_counter_leaves_state_1(10);

	TCollapser collapser(trees);
	std::vector<std::string> dimension_names_to_keep = {"molecules", "species"};
	collapser.initialize(dimension_names_to_keep, "lotus");

	EXPECT_EQ(collapser.x_is_one({0, 1, 0}, true), true);
	EXPECT_EQ(collapser.x_is_one({0, 1, 0}, false), true);
}

TEST(Collapser, collpase_2_dims) {
	std::vector<TTree> trees;

	trees.emplace_back(0, "../tests/test_data/molecules.tsv", "molecules");
	trees.emplace_back(1, "../tests/test_data/loading_tree.tsv", "tissues");
	trees.emplace_back(2, "../tests/test_data/loading_tree.tsv", "species");

	for (auto &tree : trees) { tree.initialize_cliques_and_Z(trees); }
	trees[1].get_clique({0, 0, 0}).update_counter_leaves_state_1(10);

	TCollapser collapser(trees);
	std::vector<std::string> dimension_names_to_keep = {"species"};
	collapser.initialize(dimension_names_to_keep, "lotus");

	EXPECT_EQ(collapser.x_is_one({0, 1, 0}, true), true);
	EXPECT_EQ(collapser.x_is_one({0, 1, 0}, false), true);
}

TEST(Collapser, collpase_2_dims_and_check_old_state_is_1) {
	std::vector<TTree> trees;

	trees.emplace_back(0, "../tests/test_data/molecules.tsv", "molecules");
	trees.emplace_back(1, "../tests/test_data/loading_tree.tsv", "tissues");
	trees.emplace_back(2, "../tests/test_data/loading_tree.tsv", "species");

	for (auto &tree : trees) { tree.initialize_cliques_and_Z(trees); }
	trees[1].get_clique({0, 0, 0}).update_counter_leaves_state_1(1);

	TCollapser collapser(trees);
	std::vector<std::string> dimension_names_to_keep = {"species"};
	collapser.initialize(dimension_names_to_keep, "lotus");

	EXPECT_EQ(collapser.x_is_one({0, 1, 0}, true), false);
	EXPECT_EQ(collapser.x_is_one({0, 1, 0}, false), true);
}
