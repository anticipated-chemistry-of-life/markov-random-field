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
// Since the migration to TSparseMatrix, fill_current_state is a member of
// TStorageYMatrix that walks a single matrix row (increment == 1, "easy") or a
// single matrix column (increment > 1, "hard"). The old sorted-vector variants
// (binary-search window / linear window) no longer exist, so there is nothing
// left to compare them against — these benchmarks now time the matrix fill
// directly and verify it against a brute-force per-cell lookup.

namespace {

// Build a TStorageYMatrix with Bernoulli(density) ones at each position.
// insert_in_Y now takes batches of *linear indices* (the index is implicit in
// the matrix position; TStorageY no longer stores it).
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
// Correctness: the matrix fill must agree with a brute-force per-cell lookup
// -------------------------------------------------------------------------

TEST(FillCurrentState_Matrix, matches_brute_force) {
	constexpr size_t dim0 = 200;
	constexpr size_t dim1 = 200;

	for (double density : {0.001, 0.05, 0.5}) {
		auto Y = make_Y({dim0, dim1}, density);
		std::vector<uint8_t> cur;
		std::vector<uint8_t> exists;
		std::vector<size_t> lin;

		// easy path: scan a full row (increment 1) -> linear = row * dim1 + k
		for (size_t row : {size_t{0}, size_t{37}, dim0 - 1}) {
			const IndexArray start = {row, 0};
			Y.fill_current_state(start, dim1, /*increment=*/1, cur, exists, lin);
			for (size_t k = 0; k < dim1; ++k) {
				const size_t linear = row * dim1 + k;
				EXPECT_EQ(static_cast<bool>(cur[k]), Y.is_one(linear))
				    << "easy row=" << row << " k=" << k;
				EXPECT_EQ(lin[k], linear);
			}
		}

		// hard path: scan a full column (increment dim1) -> linear = k * dim1 + col
		for (size_t col : {size_t{0}, size_t{37}, dim1 - 1}) {
			const IndexArray start = {0, col};
			Y.fill_current_state(start, dim0, /*increment=*/dim1, cur, exists, lin);
			for (size_t k = 0; k < dim0; ++k) {
				const size_t linear = k * dim1 + col;
				EXPECT_EQ(static_cast<bool>(cur[k]), Y.is_one(linear))
				    << "hard col=" << col << " k=" << k;
				EXPECT_EQ(lin[k], linear);
			}
		}
	}
}

// -------------------------------------------------------------------------
// Benchmarks
// -------------------------------------------------------------------------

// Easy path (increment = 1): scan one full matrix row.
TEST(Benchmark_FillCurrentState, easy_path) {
	constexpr size_t dim0  = 1000;
	constexpr size_t dim1  = 1000;
	constexpr size_t total = dim0 * dim1;

	const IndexArray start = {0, 0}; // scan row 0

	std::cout << "\n=== fill_current_state — easy path (increment=1, along last dim) ===\n";
	std::cout << "    container: " << dim0 << " × " << dim1 << " = " << total << " total"
	          << "  n_nodes=" << dim1 << "\n\n";

	for (double density : {0.001, 0.01, 0.05, 0.10, 0.30, 0.50}) {
		auto Y = make_Y({dim0, dim1}, density);
		std::vector<uint8_t> cur;
		std::vector<uint8_t> exists;
		std::vector<size_t> lin;
		size_t sink = 0;
		auto r      = timed([&] {
			Y.fill_current_state(start, dim1, /*increment=*/1, cur, exists, lin);
			sink += cur[0]; // prevent dead-code elimination
		});
		report("density=" + std::to_string(density) +
		           "  stored=" + std::to_string(Y.number_of_ones()),
		       r);
		(void)sink;
	}
}

// Hard path (increment = dim1): scan one full matrix column. This is the case the
// migration made cheap — a direct column walk instead of binary-search striding.
TEST(Benchmark_FillCurrentState, hard_path) {
	constexpr size_t dim0      = 1000;
	constexpr size_t dim1      = 1000;
	constexpr size_t total     = dim0 * dim1;
	constexpr size_t increment = dim1; // stride between consecutive row elements

	const IndexArray start = {0, 0}; // scan column 0

	std::cout << "\n=== fill_current_state — hard path (increment=" << increment
	          << ", non-last dim) ===\n";
	std::cout << "    container: " << dim0 << " × " << dim1 << " = " << total << " total"
	          << "  n_nodes=" << dim0 << "\n\n";

	for (double density : {0.001, 0.01, 0.05, 0.10, 0.30, 0.50}) {
		auto Y = make_Y({dim0, dim1}, density);
		std::vector<uint8_t> cur;
		std::vector<uint8_t> exists;
		std::vector<size_t> lin;
		size_t sink = 0;
		auto r      = timed([&] {
			Y.fill_current_state(start, dim0, increment, cur, exists, lin);
			sink += cur[0];
		});
		report("density=" + std::to_string(density) +
		           "  stored=" + std::to_string(Y.number_of_ones()),
		       r);
		(void)sink;
	}
}

// Easy vs. hard side-by-side at several densities.
TEST(Benchmark_FillCurrentState, easy_vs_hard_comparison) {
	constexpr size_t dim0      = 1000;
	constexpr size_t dim1      = 1000;
	constexpr size_t increment = dim1;

	const IndexArray start = {0, 0};

	std::cout << "\n=== easy vs. hard comparison  (" << dim0 << "×" << dim1 << ") ===\n\n";
	std::cout << "    " << std::left << std::setw(10) << "density" << std::setw(12) << "stored"
	          << std::setw(18) << "easy (µs/call)" << std::setw(18) << "hard (µs/call)"
	          << "ratio (hard/easy)\n";
	std::cout << "    " << std::string(70, '-') << "\n";

	for (double density : {0.001, 0.01, 0.05, 0.10, 0.30, 0.50}) {
		auto Y = make_Y({dim0, dim1}, density);
		std::vector<uint8_t> cur;
		std::vector<uint8_t> exists;
		std::vector<size_t> lin;

		size_t sink = 0;
		auto easy   = timed([&] {
			Y.fill_current_state(start, dim1, /*increment=*/1, cur, exists, lin);
			sink += cur[0];
		});
		auto hard   = timed([&] {
			Y.fill_current_state(start, dim0, increment, cur, exists, lin);
			sink += cur[0];
		});
		(void)sink;

		const double ratio = hard.us_per_call / easy.us_per_call;
		std::cout << "    " << std::left << std::fixed << std::setprecision(2) << std::setw(10)
		          << density << std::setw(12) << Y.number_of_ones() << std::setw(18)
		          << easy.us_per_call << std::setw(18) << hard.us_per_call << ratio << "×\n";
	}
}
