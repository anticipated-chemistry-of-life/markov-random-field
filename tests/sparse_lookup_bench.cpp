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
// Reading a run of cells is the whole of the sparse path's cost model. A point lookup in a
// sorted-vector matrix costs a search of one line, and an update now pays that search once per
// cell -- where a window used to walk the line once and answer every read from what it found. The
// window is gone, so this is what the sparse backend costs today, and it is the number a hash-map
// backing has to beat.
//
// Which line a lookup searches follows from the shape and not from the run: `get` searches
// whichever of the cell's row and column holds fewer entries. So a run along the last dimension
// ("easy") and one along the first ("hard") search the same kind of line, and the two are timed
// apart to say whether that is true in practice as well as on paper.
//
// Correctness of a run is a conformance question and is asserted over generated shapes in
// tests/TStorageConformance_Tests.cpp. The one check here guards the benchmark itself: it says the
// runs being timed, at these sizes and densities, read what the field holds.

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

		// easy path: a whole matrix row (stride 1) -> linear = row * dim1 + k
		for (size_t row : {size_t{0}, size_t{37}, dim0 - 1}) {
			for (size_t k = 0; k < dim1; ++k) {
				const size_t linear = row * dim1 + k;
				EXPECT_EQ(Y.get_multi_dimensional_index(linear), (IndexArray{row, k}))
				    << "easy row=" << row << " k=" << k;
			}
		}

		// hard path: a whole matrix column (stride dim1) -> linear = k * dim1 + col
		for (size_t col : {size_t{0}, size_t{37}, dim1 - 1}) {
			for (size_t k = 0; k < dim0; ++k) {
				const size_t linear = k * dim1 + col;
				EXPECT_EQ(Y.get_multi_dimensional_index(linear), (IndexArray{k, col}))
				    << "hard col=" << col << " k=" << k;
			}
		}
	}
}

// -------------------------------------------------------------------------
// Benchmarks
// -------------------------------------------------------------------------

// Easy path (stride 1): read one full matrix row, one cell at a time.
TEST(Benchmark_SparseLookup, easy_path) {
	constexpr size_t dim0  = 1000;
	constexpr size_t dim1  = 1000;
	constexpr size_t total = dim0 * dim1;

	std::cout << "\n=== a run of point lookups — easy path (stride=1, along last dim) ===\n";
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

// Hard path (stride = dim1): read one full matrix column, one cell at a time.
TEST(Benchmark_SparseLookup, hard_path) {
	constexpr size_t dim0   = 1000;
	constexpr size_t dim1   = 1000;
	constexpr size_t total  = dim0 * dim1;
	constexpr size_t stride = dim1; // one row width

	std::cout << "\n=== a run of point lookups — hard path (stride=" << stride
	          << ", non-last dim) ===\n";
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

// Easy vs. hard side-by-side at several densities.
TEST(Benchmark_SparseLookup, easy_vs_hard_comparison) {
	constexpr size_t dim0   = 1000;
	constexpr size_t dim1   = 1000;
	constexpr size_t stride = dim1;

	std::cout << "\n=== easy vs. hard comparison  (" << dim0 << "×" << dim1 << ") ===\n\n";
	std::cout << "    " << std::left << std::setw(10) << "density" << std::setw(12) << "stored"
	          << std::setw(18) << "easy (µs/call)" << std::setw(18) << "hard (µs/call)"
	          << "ratio (hard/easy)\n";
	std::cout << "    " << std::string(70, '-') << "\n";

	for (double density : {0.001, 0.01, 0.05, 0.10, 0.30, 0.50}) {
		auto Y = make_Y({dim0, dim1}, density);

		size_t sink = 0;
		auto easy   = timed([&] { sink += read_run(Y, /*start=*/0, dim1, /*stride=*/1); });
		auto hard   = timed([&] { sink += read_run(Y, /*start=*/0, dim0, stride); });
		(void)sink;

		const double ratio = hard.us_per_call / easy.us_per_call;
		std::cout << "    " << std::left << std::fixed << std::setprecision(2) << std::setw(10)
		          << density << std::setw(12) << Y.number_of_ones() << std::setw(18)
		          << easy.us_per_call << std::setw(18) << hard.us_per_call << ratio << "×\n";
	}
}
