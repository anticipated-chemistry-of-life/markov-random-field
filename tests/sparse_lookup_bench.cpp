#include "constants.h"
#include "storages/cell_write.h"
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
// Reading a run of cells is the whole of the sparse path's cost model. A point lookup costs one
// hash of the linear index, and an update pays that once per cell -- where a sorted-vector matrix
// used to pay a search of one line, and a window before it walked the line once and answered every
// read from what it found.
//
// The cost no longer follows from which line a cell is in: a map keyed by the linear index knows
// nothing of rows and columns. What is left is the memory the run touches, so a run along the last
// dimension (consecutive linear indices) and one along the first (a row width apart) are timed
// apart to say what that costs, and how the answer moves with the density.
//
// Correctness of a run is a conformance question and is asserted over generated shapes in
// tests/TStorageConformance_Tests.cpp. The one check here guards the benchmark itself: it says the
// runs being timed, at these sizes and densities, read what the field holds.

namespace {

// Build a TStorageYMatrix with Bernoulli(density) ones at each position.
// insert_in_Y takes batches of *linear indices*, which is what the map is keyed
// by; TStorageY carries a state and a counter and no index of its own.
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

/// Reads a run of cells, exactly as an update reads one: a start, a count and a stride, and a
/// point lookup per cell. Returns the number of ones, so nothing here is dead code.
size_t read_run(const TStorageYMatrix &Y, size_t start, size_t n_cells, size_t stride) {
	size_t n_ones = 0;
	for (size_t k = 0; k < n_cells; ++k) { n_ones += Y.is_one(start + k * stride); }
	return n_ones;
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
// Correctness: a run must read what the point lookups answer
// -------------------------------------------------------------------------

TEST(SparseLookup_Matrix, a_run_reads_the_cells_its_arithmetic_names) {
	constexpr size_t dim0 = 200;
	constexpr size_t dim1 = 200;

	for (double density : {0.001, 0.05, 0.5}) {
		auto Y = make_Y({dim0, dim1}, density);

		// along the last dimension: a whole row (stride 1) -> linear = row * dim1 + k
		for (size_t row : {size_t{0}, size_t{37}, dim0 - 1}) {
			for (size_t k = 0; k < dim1; ++k) {
				const size_t linear = row * dim1 + k;
				EXPECT_EQ(Y.get_multi_dimensional_index(linear), (IndexArray{row, k}))
				    << "row=" << row << " k=" << k;
			}
		}

		// along the first dimension: a whole column (stride dim1) -> linear = k * dim1 + col
		for (size_t col : {size_t{0}, size_t{37}, dim1 - 1}) {
			for (size_t k = 0; k < dim0; ++k) {
				const size_t linear = k * dim1 + col;
				EXPECT_EQ(Y.get_multi_dimensional_index(linear), (IndexArray{k, col}))
				    << "col=" << col << " k=" << k;
			}
		}
	}
}

// -------------------------------------------------------------------------
// Benchmarks
// -------------------------------------------------------------------------

// Stride 1: read one full row, one cell at a time.
TEST(Benchmark_SparseLookup, a_run_along_the_last_dimension) {
	constexpr size_t dim0  = 1000;
	constexpr size_t dim1  = 1000;
	constexpr size_t total = dim0 * dim1;

	std::cout << "\n=== a run of point lookups — stride=1, along the last dimension ===\n";
	std::cout << "    container: " << dim0 << " × " << dim1 << " = " << total << " total"
	          << "  n_cells=" << dim1 << "\n\n";

	for (double density : {0.001, 0.01, 0.05, 0.10, 0.30, 0.50}) {
		auto Y      = make_Y({dim0, dim1}, density);
		size_t sink = 0;
		auto r      = timed([&] { sink += read_run(Y, /*start=*/0, dim1, /*stride=*/1); });
		report("density=" + std::to_string(density) +
		           "  stored=" + std::to_string(Y.number_of_ones()),
		       r);
		(void)sink;
	}
}

// Stride = dim1: read one full column, one cell at a time.
TEST(Benchmark_SparseLookup, a_run_along_the_first_dimension) {
	constexpr size_t dim0   = 1000;
	constexpr size_t dim1   = 1000;
	constexpr size_t total  = dim0 * dim1;
	constexpr size_t stride = dim1; // one row width

	std::cout << "\n=== a run of point lookups — along the first dimension, stride=" << stride
	          << " ===\n";
	std::cout << "    container: " << dim0 << " × " << dim1 << " = " << total << " total"
	          << "  n_cells=" << dim0 << "\n\n";

	for (double density : {0.001, 0.01, 0.05, 0.10, 0.30, 0.50}) {
		auto Y      = make_Y({dim0, dim1}, density);
		size_t sink = 0;
		auto r      = timed([&] { sink += read_run(Y, /*start=*/0, dim0, stride); });
		report("density=" + std::to_string(density) +
		           "  stored=" + std::to_string(Y.number_of_ones()),
		       r);
		(void)sink;
	}
}

// The two strides side by side at several densities.
TEST(Benchmark_SparseLookup, the_two_strides_compared) {
	constexpr size_t dim0   = 1000;
	constexpr size_t dim1   = 1000;
	constexpr size_t stride = dim1;

	std::cout << "\n=== the two strides compared  (" << dim0 << "×" << dim1 << ") ===\n\n";
	std::cout << "    " << std::left << std::setw(10) << "density" << std::setw(12) << "stored"
	          << std::setw(18) << "stride 1 (µs)" << std::setw(18) << "stride dim1 (µs)"
	          << "ratio\n";
	std::cout << "    " << std::string(70, '-') << "\n";

	for (double density : {0.001, 0.01, 0.05, 0.10, 0.30, 0.50}) {
		auto Y = make_Y({dim0, dim1}, density);

		size_t sink = 0;
		auto by_row = timed([&] { sink += read_run(Y, /*start=*/0, dim1, /*stride=*/1); });
		auto by_col = timed([&] { sink += read_run(Y, /*start=*/0, dim0, stride); });
		(void)sink;

		const double ratio = by_col.us_per_call / by_row.us_per_call;
		std::cout << "    " << std::left << std::fixed << std::setprecision(2) << std::setw(10)
		          << density << std::setw(12) << Y.number_of_ones() << std::setw(18)
		          << by_row.us_per_call << std::setw(18) << by_col.us_per_call << ratio << "×\n";
	}
}
