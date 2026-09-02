#include "constants.h"
#include "storages/y_storage/TStorageYMatrix.h"
#include "gtest/gtest.h"
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

// -------------------------------------------------------------------------
// Helpers
// -------------------------------------------------------------------------
//
// Opening a window is the whole of the sparse path's cost model. A point lookup in a sparse matrix
// costs a search, so the sparse window walks one line on open and answers every later read from
// what it found. Everything an update reads goes through that walk, and how it scales with density
// is what decides whether the sparse backend is usable at scale.
//
// Which line the window walks follows from the stride: a stride of one is a matrix row ("easy"),
// and a stride of one row width is a matrix column ("hard"). The two are timed apart because they
// are different walks over the same data.
//
// Correctness of the walk is a conformance question and is asserted over generated shapes in
// tests/TStorageConformance_Tests.cpp. The one check here guards the benchmark itself: it says the
// windows being timed, at these sizes and densities, hold what the field holds.

namespace {

// Build a TStorageYMatrix with Bernoulli(density) ones at each position.
// insert_in_Y takes batches of *linear indices* (the index is implicit in the
// matrix position; TStorageY does not store it).
TStorageYMatrix make_Y(const std::vector<size_t> &dims, double density, uint64_t seed = 42) {
	TStorageYMatrix Y;
	Y.initialize(/*n_iterations=*/1000, dims);
	const size_t total = Y.total_size_of_container_space();

	std::mt19937_64 rng(seed);
	std::bernoulli_distribution dist(density);

	std::vector<size_t> linear_indices;
	linear_indices.reserve(static_cast<size_t>(total * density * 1.2));
	for (size_t i = 0; i < total; ++i) {
		if (dist(rng)) { linear_indices.push_back(i); }
	}
	// already sorted because we iterate indices in ascending order
	std::vector<std::vector<size_t>> batch = {std::move(linear_indices)};
	Y.insert_in_Y(batch);
	return Y;
}

using Clock = std::chrono::high_resolution_clock;

struct BenchResult {
	double us_per_call;
	size_t reps;
};

// Run fn() reps times, preceded by a short warm-up, and return µs/call.
template<typename Fn> BenchResult timed(Fn &&fn, size_t reps = 300) {
	for (size_t i = 0; i < 5; ++i) { fn(); } // warm-up
	auto t0 = Clock::now();
	for (size_t i = 0; i < reps; ++i) { fn(); }
	auto t1   = Clock::now();
	double us = std::chrono::duration<double, std::micro>(t1 - t0).count();
	return {us / static_cast<double>(reps), reps};
}

void report(const std::string &label, const BenchResult &r) {
	std::cout << "    " << std::left << std::setw(45) << label << std::fixed << std::setprecision(2)
	          << r.us_per_call << " µs/call  (" << r.reps << " reps)\n";
}

} // namespace

// -------------------------------------------------------------------------
// Correctness: a window must hold what the point lookups answer
// -------------------------------------------------------------------------

TEST(OpenWindow_Matrix, matches_brute_force) {
	constexpr size_t dim0 = 200;
	constexpr size_t dim1 = 200;

	for (double density : {0.001, 0.05, 0.5}) {
		auto Y = make_Y({dim0, dim1}, density);

		// easy path: a whole matrix row (stride 1) -> linear = row * dim1 + k
		for (size_t row : {size_t{0}, size_t{37}, dim0 - 1}) {
			auto window = Y.open_window(IndexArray{row, 0}, dim1, /*stride=*/1);
			for (size_t k = 0; k < dim1; ++k) {
				const size_t linear = row * dim1 + k;
				EXPECT_EQ(window.linear_index(k), linear) << "easy row=" << row << " k=" << k;
				EXPECT_EQ(window.is_one(k), Y.is_one(linear)) << "easy row=" << row << " k=" << k;
			}
		}

		// hard path: a whole matrix column (stride dim1) -> linear = k * dim1 + col
		for (size_t col : {size_t{0}, size_t{37}, dim1 - 1}) {
			auto window = Y.open_window(IndexArray{0, col}, dim0, /*stride=*/dim1);
			for (size_t k = 0; k < dim0; ++k) {
				const size_t linear = k * dim1 + col;
				EXPECT_EQ(window.linear_index(k), linear) << "hard col=" << col << " k=" << k;
				EXPECT_EQ(window.is_one(k), Y.is_one(linear)) << "hard col=" << col << " k=" << k;
			}
		}
	}
}

// -------------------------------------------------------------------------
// Benchmarks
// -------------------------------------------------------------------------

// Easy path (stride 1): open a window over one full matrix row.
TEST(Benchmark_OpenWindow, easy_path) {
	constexpr size_t dim0  = 1000;
	constexpr size_t dim1  = 1000;
	constexpr size_t total = dim0 * dim1;

	const IndexArray start = {0, 0}; // row 0

	std::cout << "\n=== open_window — easy path (stride=1, along last dim) ===\n";
	std::cout << "    container: " << dim0 << " × " << dim1 << " = " << total << " total"
	          << "  n_cells=" << dim1 << "\n\n";

	for (double density : {0.001, 0.01, 0.05, 0.10, 0.30, 0.50}) {
		auto Y      = make_Y({dim0, dim1}, density);
		size_t sink = 0;
		auto r      = timed([&] {
			auto window = Y.open_window(start, dim1, /*stride=*/1);
			sink += window.is_one(0); // prevent dead-code elimination
		});
		report("density=" + std::to_string(density) +
		           "  stored=" + std::to_string(Y.number_of_ones()),
		       r);
		(void)sink;
	}
}

// Hard path (stride = dim1): open a window over one full matrix column.
TEST(Benchmark_OpenWindow, hard_path) {
	constexpr size_t dim0   = 1000;
	constexpr size_t dim1   = 1000;
	constexpr size_t total  = dim0 * dim1;
	constexpr size_t stride = dim1; // one row width

	const IndexArray start = {0, 0}; // column 0

	std::cout << "\n=== open_window — hard path (stride=" << stride << ", non-last dim) ===\n";
	std::cout << "    container: " << dim0 << " × " << dim1 << " = " << total << " total"
	          << "  n_cells=" << dim0 << "\n\n";

	for (double density : {0.001, 0.01, 0.05, 0.10, 0.30, 0.50}) {
		auto Y      = make_Y({dim0, dim1}, density);
		size_t sink = 0;
		auto r      = timed([&] {
			auto window = Y.open_window(start, dim0, stride);
			sink += window.is_one(0);
		});
		report("density=" + std::to_string(density) +
		           "  stored=" + std::to_string(Y.number_of_ones()),
		       r);
		(void)sink;
	}
}

// Easy vs. hard side-by-side at several densities.
TEST(Benchmark_OpenWindow, easy_vs_hard_comparison) {
	constexpr size_t dim0   = 1000;
	constexpr size_t dim1   = 1000;
	constexpr size_t stride = dim1;

	const IndexArray start = {0, 0};

	std::cout << "\n=== easy vs. hard comparison  (" << dim0 << "×" << dim1 << ") ===\n\n";
	std::cout << "    " << std::left << std::setw(10) << "density" << std::setw(12) << "stored"
	          << std::setw(18) << "easy (µs/call)" << std::setw(18) << "hard (µs/call)"
	          << "ratio (hard/easy)\n";
	std::cout << "    " << std::string(70, '-') << "\n";

	for (double density : {0.001, 0.01, 0.05, 0.10, 0.30, 0.50}) {
		auto Y = make_Y({dim0, dim1}, density);

		size_t sink = 0;
		auto easy   = timed([&] {
			auto window = Y.open_window(start, dim1, /*stride=*/1);
			sink += window.is_one(0);
		});
		auto hard   = timed([&] {
			auto window = Y.open_window(start, dim0, stride);
			sink += window.is_one(0);
		});
		(void)sink;

		const double ratio = hard.us_per_call / easy.us_per_call;
		std::cout << "    " << std::left << std::fixed << std::setprecision(2) << std::setw(10)
		          << density << std::setw(12) << Y.number_of_ones() << std::setw(18)
		          << easy.us_per_call << std::setw(18) << hard.us_per_call << ratio << "×\n";
	}
}
