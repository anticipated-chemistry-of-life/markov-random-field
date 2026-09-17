//
// One cache line per thread, for the accumulators a parallel region writes once per cell.
//

#pragma once

/// `T`, alone on a cache line.
///
/// A `std::vector<T>` of per-thread accumulators puts two or three threads' slots inside one
/// 64-byte cache line, so a write by one thread takes the line away from every other thread
/// holding a slot in it. At one write per leaf pair that is the block update's largest single
/// cost: the two accumulator writes in `block_update::update_cell` carried 45% of its profile
/// samples, and a microbenchmark of those two writes alone cost 113-125 ns per leaf pair on shared
/// lines against 5.9 ns on private ones (16 threads, one per physical core).
///
/// The line width is written out rather than taken from
/// `std::hardware_destructive_interference_size`, whose use GCC warns about because the value is
/// an ABI commitment: a header that changes it breaks every translation unit compiled against the
/// old one. Every target this builds for -- x86-64 and AArch64 -- has a 64-byte line.
///
/// This is a stopgap. A per-thread accumulator written once per cell is the wrong shape whatever
/// its alignment, because the sum it holds depends on which cells the schedule handed that thread;
/// `TTree::log_node_state_density` keeps one slot per clique instead, for exactly that reason. The
/// accumulators this wraps move to one slot per row, written once per row, and then nothing here
/// has a user left.
template<typename T> struct alignas(64) TCacheLineSlot {
	T value{};
};
