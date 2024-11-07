
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
	std::vector<size_t> dimensions{100, 200};
	TStorageYVector Y(1000, dimensions);

	// We fill randomly the TStorageYVector where the values are between 0 and 100*100
	for (size_t i = 0; i < coretools::containerProduct(dimensions); i++) {
		auto sample = d(gen);
		if (sample) { Y.set_to_one(i); }
	}

	std::vector<size_t> multi_index{0, 0};
	auto res = fill_current_state_easy(Y, multi_index, 200);
}
