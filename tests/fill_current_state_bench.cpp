#include "smart_binary_search.h"
#include "storages/y_storage/TStorageY.h"
#include "storages/y_storage/TStorageYVector.h"
#include "gtest/gtest.h"
#include <chrono>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

// -------------------------------------------------------------------------
// Helpers
// -------------------------------------------------------------------------

namespace {

// Build a TStorageYVector with Bernoulli(density) elements at each position.
// Uses insert_in_Y (single sort+merge) instead of repeated insert_one to avoid
// the O(N²) overhead of repeated sorted insertion.
TStorageYVector make_Y(const std::vector<size_t> &dims, double density, uint64_t seed = 42) {
	TStorageYVector Y(/*n_iterations=*/1000, dims);
	const size_t total = Y.total_size_of_container_space();

	std::mt19937_64 rng(seed);
	std::bernoulli_distribution dist(density);

	std::vector<TStorageY> elems;
	elems.reserve(static_cast<size_t>(total * density * 1.2));
	for (size_t i = 0; i < total; ++i) {
		if (dist(rng)) { elems.emplace_back(i); } // constructor sets state=1
	}
	// elems is already sorted because we iterate indices in ascending order
	std::vector<std::vector<TStorageY>> batch = {std::move(elems)};
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
// Benchmarks
// -------------------------------------------------------------------------

// Benchmark the easy path (AlongLastDim = true, increment = 1).
// Simulates scanning all K leaves along the last tree dimension for a fixed
// position in every other dimension.
TEST(Benchmark_FillCurrentState, easy_path) {
	// 2-D container: dim0 rows × dim1 columns. We scan one full row (n_nodes=dim1, increment=1).
	constexpr size_t dim0  = 1000;
	constexpr size_t dim1  = 1000;
	constexpr size_t total = dim0 * dim1;

	const std::vector<size_t> start = {0, 0}; // scan row 0

	std::cout << "\n=== fill_current_state — easy path (increment=1, along last dim) ===\n";
	std::cout << "    container: " << dim0 << " × " << dim1 << " = " << total << " total"
	          << "  n_nodes=" << dim1 << "\n\n";

	for (double density : {0.001, 0.01, 0.05, 0.10, 0.30, 0.50}) {
		auto Y      = make_Y({dim0, dim1}, density);
		size_t sink = 0;
		auto r      = timed([&] {
			auto res = fill_current_state<true>(Y, dim1, start, /*increment=*/1, total);
			sink += static_cast<size_t>(res.current_state[0]); // prevent dead-code elimination
		});
		report("density=" + std::to_string(density) + "  stored=" + std::to_string(Y.size()), r);
		(void)sink;
	}
}

// Benchmark the hard path (AlongLastDim = false, increment > 1).
// Simulates scanning all leaves along a non-last tree dimension, which means
// striding through the linear container space by `increment` at each step.
TEST(Benchmark_FillCurrentState, hard_path) {
	// 2-D container: dim0 rows × dim1 columns. We scan one full column
	// (n_nodes=dim0, increment=dim1 — stride between consecutive rows).
	constexpr size_t dim0      = 1000;
	constexpr size_t dim1      = 1000;
	constexpr size_t total     = dim0 * dim1;
	constexpr size_t increment = dim1; // stride between consecutive row elements

	const std::vector<size_t> start = {0, 0}; // scan column 0

	std::cout << "\n=== fill_current_state — hard path (increment=" << increment
	          << ", non-last dim) ===\n";
	std::cout << "    container: " << dim0 << " × " << dim1 << " = " << total << " total"
	          << "  n_nodes=" << dim0 << "\n\n";

	for (double density : {0.001, 0.01, 0.05, 0.10, 0.30, 0.50}) {
		auto Y      = make_Y({dim0, dim1}, density);
		size_t sink = 0;
		auto r      = timed([&] {
			auto res = fill_current_state<false>(Y, dim0, start, increment, total);
			sink += static_cast<size_t>(res.current_state[0]);
		});
		report("density=" + std::to_string(density) + "  stored=" + std::to_string(Y.size()), r);
		(void)sink;
	}
}

// -------------------------------------------------------------------------
// Alternative implementations of fill_current_state_hard for comparison
// -------------------------------------------------------------------------
//
// All variants below must produce the same current_state and exists_in_container
// as fill_current_state_hard (verified before timing).  They differ only in
// how they locate each query in the sorted container.
//
// Naming:
//   baseline  – current fill_current_state_hard (density probe + binary search in window)
//   opt1      – no density estimation; binary_search from current idx to end  (O(log N) / query)
//   linear    – no density estimation; linear scan from current idx           (O(increment·p) /
//   query) opt3      – same density probe as baseline; linear scan *inside* the window (replaces
//   inner bs)

namespace {

// Option 1 (proposed): drop density estimation entirely; for each query do a
// binary search from the current container position to the end.
template<typename Container>
CurrentStateResult hard_opt1_lower_bound(const Container &container, size_t n_nodes,
                                         const std::vector<size_t> &start_index, size_t increment,
                                         size_t /*total_size*/) {
	std::vector<int> state(n_nodes, false);
	std::vector<int> exists(n_nodes, false);
	std::vector<size_t> idx_vec(n_nodes, container.size());
	const size_t N = container.size();

	const auto lin_start         = container.get_linear_index_in_container_space(start_index);
	auto [found, idx, lin, last] = binary_search(container, lin_start, size_t{0}, lin_start + 1);
	idx_vec[0]                   = idx;
	if (found) {
		state[0]  = container[idx].is_one();
		exists[0] = true;
	}

	for (size_t i = 1; i < n_nodes; ++i) {
		if (last) return {state, exists, idx_vec};
		const auto target = lin_start + i * increment;
		if (target < lin) {
			idx_vec[i] = idx;
			continue;
		}
		std::tie(found, idx, lin, last) = binary_search(container, target, idx, N);
		idx_vec[i]                      = idx;
		if (found) {
			state[i]  = container[idx].is_one();
			exists[i] = true;
		}
	}
	return {state, exists, idx_vec};
}

// Option 3 (proposed): keep the density probe to narrow the search range, but
// replace the binary_search *inside the window* with a linear scan.
// Wins when the window [lower_bound, upper_bound] is small (sparse + large increment).
template<typename Container>
CurrentStateResult hard_opt3_linear_window(const Container &container, size_t n_nodes,
                                           const std::vector<size_t> &start_index, size_t increment,
                                           size_t total_size) {
	std::vector<int> state(n_nodes, false);
	std::vector<int> exists(n_nodes, false);
	std::vector<size_t> idx_vec(n_nodes, container.size());
	const size_t N = container.size();

	const auto lin_start         = container.get_linear_index_in_container_space(start_index);
	auto [found, idx, lin, last] = binary_search(container, lin_start, size_t{0}, lin_start + 1);
	idx_vec[0]                   = idx;
	if (found) {
		state[0]  = container[idx].is_one();
		exists[0] = true;
	}

	const double p   = static_cast<double>(N) / static_cast<double>(total_size);
	const double ip  = static_cast<double>(increment) * p;
	const double sig = 2.0 * std::sqrt(ip * (1.0 - p));
	const auto jr    = static_cast<size_t>(std::ceil(ip + sig));
	const size_t jl  = static_cast<size_t>(std::max(0.0, std::floor(ip - sig)));

	for (size_t i = 1; i < n_nodes; ++i) {
		const auto target = lin_start + i * increment;
		if (target < lin) {
			idx_vec[i] = idx;
			continue;
		}
		if (last) return {state, exists, idx_vec};

		auto ub = idx + jr;
		if (ub >= N) ub = N - 1;
		auto lb = idx + jl;
		if (lb >= N) lb = N - 1;
		const auto ub_lin = container[ub].get_linear_index_in_container_space();
		const auto lb_lin = container[lb].get_linear_index_in_container_space();

		// exact hits at the two probe positions (same as baseline)
		if (target == ub_lin) {
			idx        = ub;
			lin        = ub_lin;
			last       = (idx == N - 1);
			state[i]   = container[idx].is_one();
			exists[i]  = true;
			idx_vec[i] = idx;
			continue;
		}
		if (target == lb_lin) {
			idx        = lb;
			lin        = lb_lin;
			last       = (idx == N - 1);
			state[i]   = container[idx].is_one();
			exists[i]  = true;
			idx_vec[i] = idx;
			continue;
		}

		if (target > ub_lin) {
			// above window: binary search from ub to end (same as baseline)
			std::tie(found, idx, lin, last) = binary_search(container, target, ub, N);
		} else if (target > lb_lin /* && target < ub_lin */) {
			// inside window: linear scan [lb, ub) instead of binary search
			size_t pos = lb;
			while (pos < ub && container[pos].get_linear_index_in_container_space() < target) {
				++pos;
			}
			const auto pos_lin = container[pos].get_linear_index_in_container_space();
			found              = (pos_lin == target);
			idx                = pos;
			lin                = pos_lin;
			last               = (pos == N - 1);
		} else {
			// below window: binary search [idx, lb) (same as baseline)
			std::tie(found, idx, lin, last) = binary_search(container, target, idx, lb);
		}
		idx_vec[i] = idx;
		if (found) {
			state[i]  = container[idx].is_one();
			exists[i] = true;
		}
	}
	return {state, exists, idx_vec};
}

// Assert that two results agree on the semantically meaningful fields.
void assert_equal(const CurrentStateResult &ref, const CurrentStateResult &got,
                  const std::string &name) {
	ASSERT_EQ(ref.current_state, got.current_state) << name << ": current_state mismatch";
	ASSERT_EQ(ref.exists_in_container, got.exists_in_container)
	    << name << ": exists_in_container mismatch";
	ASSERT_EQ(ref.index_in_TStorageVector, got.index_in_TStorageVector)
	    << name << ": index_in_TStorageVector mismatch";
}

} // namespace

// -------------------------------------------------------------------------
// Variant comparison benchmark
// -------------------------------------------------------------------------

// Times four implementations of fill_current_state_hard against each other.
// First verifies correctness (all variants must match the baseline), then
// measures µs/call at several densities so trade-offs are visible.
TEST(Benchmark_FillCurrentState, hard_path_variants) {
	constexpr size_t dim0           = 1000;
	constexpr size_t dim1           = 1000;
	constexpr size_t total          = dim0 * dim1;
	constexpr size_t increment      = dim1;
	constexpr size_t n_nodes        = dim0;
	const std::vector<size_t> start = {0, 0};

	std::cout << "\n=== fill_current_state_hard: variant comparison ===\n";
	std::cout << "    container: " << dim0 << "×" << dim1 << " = " << total
	          << "  n_nodes=" << n_nodes << "  increment=" << increment << "\n\n";
	std::cout << "    " << std::left << std::setw(8) << "density" << std::setw(10) << "stored"
	          << std::setw(14) << "baseline" << std::setw(14) << "opt1_lb" << std::setw(14)
	          << "linear" << std::setw(14) << "opt3_win"
	          << "\n";
	std::cout << "    " << std::string(74, '-') << "\n";

	for (double density : {0.001, 0.01, 0.05, 0.10, 0.30, 0.50}) {
		auto Y = make_Y({dim0, dim1}, density);

		// --- correctness check (run once, asserts fire on mismatch) ---
		const auto ref  = fill_current_state<false>(Y, n_nodes, start, increment, total);
		const auto v1   = hard_opt1_lower_bound(Y, n_nodes, start, increment, total);
		const auto vlin = hard_linear_scan(Y, n_nodes, start, increment, total);
		const auto v3   = hard_opt3_linear_window(Y, n_nodes, start, increment, total);
		assert_equal(ref, v1, "opt1_lower_bound");
		assert_equal(ref, vlin, "linear_scan");
		assert_equal(ref, v3, "opt3_linear_window");

		// --- timing ---
		size_t sink = 0;
		auto t_base = timed([&] {
			auto r = fill_current_state<false>(Y, n_nodes, start, increment, total);
			sink += static_cast<size_t>(r.current_state[0]);
		});
		auto t_opt1 = timed([&] {
			auto r = hard_opt1_lower_bound(Y, n_nodes, start, increment, total);
			sink += static_cast<size_t>(r.current_state[0]);
		});
		auto t_lin  = timed([&] {
			auto r = hard_linear_scan(Y, n_nodes, start, increment, total);
			sink += static_cast<size_t>(r.current_state[0]);
		});
		auto t_opt3 = timed([&] {
			auto r = hard_opt3_linear_window(Y, n_nodes, start, increment, total);
			sink += static_cast<size_t>(r.current_state[0]);
		});
		(void)sink;

		auto fmt = [](double v) {
			std::ostringstream s;
			s << std::fixed << std::setprecision(2) << v << " µs";
			return s.str();
		};

		std::cout << "    " << std::left << std::fixed << std::setprecision(4) << std::setw(8)
		          << density << std::setw(10) << Y.size() << std::setw(14)
		          << fmt(t_base.us_per_call) << std::setw(14) << fmt(t_opt1.us_per_call)
		          << std::setw(14) << fmt(t_lin.us_per_call) << std::setw(14)
		          << fmt(t_opt3.us_per_call) << "\n";
	}
}

// Benchmark both paths at a realistic size with multiple densities so that
// easy vs. hard can be directly compared side-by-side.
TEST(Benchmark_FillCurrentState, easy_vs_hard_comparison) {
	constexpr size_t dim0      = 1000;
	constexpr size_t dim1      = 1000;
	constexpr size_t total     = dim0 * dim1;
	constexpr size_t increment = dim1;

	const std::vector<size_t> start = {0, 0};

	std::cout << "\n=== easy vs. hard comparison  (" << dim0 << "×" << dim1 << ") ===\n\n";
	std::cout << "    " << std::left << std::setw(10) << "density" << std::setw(12) << "stored"
	          << std::setw(18) << "easy (µs/call)" << std::setw(18) << "hard (µs/call)"
	          << "ratio (hard/easy)\n";
	std::cout << "    " << std::string(70, '-') << "\n";

	for (double density : {0.001, 0.01, 0.05, 0.10, 0.30, 0.50}) {
		auto Y = make_Y({dim0, dim1}, density);

		size_t sink = 0;
		auto easy   = timed([&] {
			auto res = fill_current_state<true>(Y, dim1, start, 1, total);
			sink += static_cast<size_t>(res.current_state[0]);
		});
		auto hard   = timed([&] {
			auto res = fill_current_state<false>(Y, dim0, start, increment, total);
			sink += static_cast<size_t>(res.current_state[0]);
		});
		(void)sink;

		const double ratio = hard.us_per_call / easy.us_per_call;
		std::cout << "    " << std::left << std::fixed << std::setprecision(2) << std::setw(10)
		          << density << std::setw(12) << Y.size() << std::setw(18) << easy.us_per_call
		          << std::setw(18) << hard.us_per_call << ratio << "×\n";
	}
}
