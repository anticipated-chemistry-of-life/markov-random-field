
#include "TStorageYVector.h"
#include "coretools/algorithms.h"
#include "coretools/devtools.h"
#include "update_current_state.h"
#include "gtest/gtest.h"
#include <chrono>
#include <iostream>
#include <random>
#include <vector>

TEST(Binary_search, easy_function) {
	// create the engine
	std::random_device rd{};
	std::mt19937 gen{rd()};

	double p = 0.01;
	std::bernoulli_distribution d(p);

	// We create a TStorageYVector object
	std::vector<size_t> dimensions{1000, 200000};
	TStorageYVector Y(1000, dimensions);

	// We fill randomly the TStorageYVector where the values are between 0 and 100*100
	for (size_t i = 0; i < coretools::containerProduct(dimensions); i++) {
		auto sample = d(gen);
		if (sample) { Y.set_to_one(i); }
	}

	std::vector<long long> vec_duration_1;
	std::vector<long long> vec_duration_2;
	for (size_t j = 0; j < dimensions[0]; ++j) {
		std::vector<size_t> multi_index_start{j, 0};
		auto t1               = std::chrono::high_resolution_clock::now();
		std::vector<bool> res = fill_current_state_easy(Y, multi_index_start, dimensions[dimensions.size() - 1]);
		auto t2               = std::chrono::high_resolution_clock::now();
		std::vector<bool> truth;
		truth.clear();
		auto linear_index = Y.get_linear_coordinate(multi_index_start);
		auto t3           = std::chrono::high_resolution_clock::now();
		for (size_t i = 0; i < dimensions[dimensions.size() - 1]; ++i) { truth.push_back(Y.is_one(linear_index + i)); }
		auto t4        = std::chrono::high_resolution_clock::now();
		auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
		auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();
		vec_duration_1.push_back(duration1);
		vec_duration_2.push_back(duration2);

		EXPECT_EQ(res, truth);
	}
	std::cout << "Average duration for fill_current_state_easy: "
	          << std::accumulate(vec_duration_1.begin(), vec_duration_1.end(), 0) / vec_duration_1.size() << " ms"
	          << std::endl;
	std::cout << "Average duration for standard binary search: "
	          << std::accumulate(vec_duration_2.begin(), vec_duration_2.end(), 0) / vec_duration_2.size() << " ms"
	          << std::endl;
}

TEST(Binary_search, hard_function) {
	// create the engine
	std::random_device rd{};
	std::mt19937 gen{rd()};

	double p = 0.01;
	std::bernoulli_distribution d(p);

	// We create a TStorageYVector object
	std::vector<size_t> dimensions{1000000, 1000};
	TStorageYVector Y(1000, dimensions);

	// We fill randomly the TStorageYVector where the values are between 0 and 100*100
	for (size_t i = 0; i < coretools::containerProduct(dimensions); i++) {
		auto sample = d(gen);
		if (sample) { Y.set_to_one(i); }
	}

	std::vector<long long> vec_duration_1;
	std::vector<long long> vec_duration_2;
	std::vector<int> vec_lower_chunk;
	std::vector<int> vec_middle_chunk;
	std::vector<int> vec_upper_chunk;
	for (size_t j = 0; j < dimensions[1]; ++j) {
		std::vector<size_t> multi_index_start{0, j};
		auto t1                                            = std::chrono::high_resolution_clock::now();
		auto [res, lower_chunk, middle_chunk, upper_chunk] = fill_current_state_hard(
		    dimensions[0], Y, multi_index_start, dimensions[1], coretools::containerProduct(dimensions));
		auto t2 = std::chrono::high_resolution_clock::now();
		std::vector<bool> truth;
		truth.clear();
		auto linear_index = Y.get_linear_coordinate(multi_index_start);
		auto t3           = std::chrono::high_resolution_clock::now();
		for (size_t i = 0; i < dimensions[0]; ++i) { truth.push_back(Y.is_one(linear_index + i * dimensions[1])); }
		auto t4 = std::chrono::high_resolution_clock::now();
		EXPECT_EQ(res, truth);

		auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
		auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();
		vec_duration_1.push_back(duration1);
		vec_duration_2.push_back(duration2);
		vec_lower_chunk.push_back(lower_chunk);
		vec_middle_chunk.push_back(middle_chunk);
		vec_upper_chunk.push_back(upper_chunk);
	}
	std::cout << "Average duration for fill_current_state_hard: "
	          << std::accumulate(vec_duration_1.begin(), vec_duration_1.end(), 0) / vec_duration_1.size() << " ms"
	          << std::endl;
	std::cout << "Average duration for standard binary search: "
	          << std::accumulate(vec_duration_2.begin(), vec_duration_2.end(), 0) / vec_duration_2.size() << " ms"
	          << std::endl;
	std::cout << "Average lower chunk: "
	          << std::accumulate(vec_lower_chunk.begin(), vec_lower_chunk.end(), 0) / vec_lower_chunk.size()
	          << std::endl;
	std::cout << "Average middle chunk: "
	          << std::accumulate(vec_middle_chunk.begin(), vec_middle_chunk.end(), 0) / vec_middle_chunk.size()
	          << std::endl;
	std::cout << "Average upper chunk: "
	          << std::accumulate(vec_upper_chunk.begin(), vec_upper_chunk.end(), 0) / vec_upper_chunk.size()
	          << std::endl;
}
