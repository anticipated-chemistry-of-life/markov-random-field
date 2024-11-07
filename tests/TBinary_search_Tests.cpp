
#include "TStorageYVector.h"
#include "coretools/algorithms.h"
#include "coretools/devtools.h"
#include "update_current_state.h"
#include "gtest/gtest.h"
#include <random>
#include <vector>

TEST(Binary_search, easy_function) {
	// create the engine
	std::random_device rd{};
	std::mt19937 gen{rd()};

	double p = 0.1;
	std::bernoulli_distribution d(p);

	// We create a TStorageYVector object
	std::vector<size_t> dimensions{10000, 1000};
	TStorageYVector Y(1000, dimensions);

	// We fill randomly the TStorageYVector where the values are between 0 and 100*100
	for (size_t i = 0; i < coretools::containerProduct(dimensions); i++) {
		auto sample = d(gen);
		if (sample) { Y.set_to_one(i); }
	}

	for (size_t j = 0; j < dimensions[0]; ++j) {
		std::vector<size_t> multi_index_start{j, 0};
		std::vector<bool> res = fill_current_state_easy(Y, multi_index_start, dimensions[dimensions.size() - 1]);
		std::vector<bool> truth;
		truth.clear();
		auto linear_index = Y.get_linear_coordinate(multi_index_start);
		for (size_t i = 0; i < dimensions[dimensions.size() - 1]; ++i) { truth.push_back(Y.is_one(linear_index + i)); }
		EXPECT_EQ(res, truth);
	}
}

TEST(Binary_search, hard_function) {
	// create the engine
	std::random_device rd{};
	std::mt19937 gen{rd()};

	double p = 0.5;
	std::bernoulli_distribution d(p);

	// We create a TStorageYVector object
	std::vector<size_t> dimensions{4, 6};
	TStorageYVector Y(1000, dimensions);

	// We fill randomly the TStorageYVector where the values are between 0 and 100*100
	for (size_t i = 0; i < coretools::containerProduct(dimensions); i++) {
		auto sample = d(gen);
		if (sample) { Y.set_to_one(i); }
	}

	for (auto y : Y.get_vector()) { std::cout << y.get_coordinate() << " "; }
	for (size_t j = 0; j < dimensions[1]; ++j) {
		std::vector<size_t> multi_index_start{0, j};
		std::vector<bool> res = fill_current_state_hard(dimensions[0], Y, multi_index_start, dimensions[1],
		                                                coretools::containerProduct(dimensions));
		std::vector<bool> truth;
		truth.clear();
		auto linear_index = Y.get_linear_coordinate(multi_index_start);
		for (size_t i = 0; i < dimensions[0]; ++i) { truth.push_back(Y.is_one(linear_index + i * dimensions[1])); }
		EXPECT_EQ(res, truth);
	}
}
