//
// The uniform a cell draws, derived from where the cell is.
//

#ifndef ACOL_TCELLUNIFORMS_H
#define ACOL_TCELLUNIFORMS_H

#include "coretools/Main/TRandomGenerator.h"
#include <cstddef>
#include <cstdint>

/// Which draw of which container a stream of cell uniforms serves.
///
/// Two containers hold a cell at the same linear index, and two draws visit the same container. The
/// label keeps their uniforms apart; ADR-0007 says why each label is there.
enum class TCellStream : uint64_t {
	field,               ///< the field update
	node_state,          ///< one tree's node-state update
	field_at_start,      ///< the field a simulated chain starts from, drawn from the two trees
	node_state_at_start, ///< one tree's node state a simulated chain starts from
};

/// The uniforms one container's cells draw at one iteration.
///
/// A cell's uniform is a hash of the seed, the stream, the tree the container belongs to, the
/// iteration and the cell's linear index. The first four mix down to one key when the object is
/// built. The cell index is the counter. The mixing function is the SplitMix64 finaliser, and the
/// counter steps by SplitMix64's own increment, so the cells of one stream walk a SplitMix64
/// sequence.
///
/// Nothing in that derivation counts the draws that came before, so a cell's number does not move
/// with the thread count or with the order an update visits cells. ADR-0007 says why the sampler
/// needs that, and why two positions share a uniform only by chance.
class TCellUniforms {
private:
	/// The odd increment SplitMix64 walks its counter by: 2^64 divided by the golden ratio.
	static constexpr uint64_t _gamma = 0x9E3779B97F4A7C15ULL;

	/// The seed, the stream, the tree and the iteration, mixed down to one word. The cell index is
	/// all that is left to vary.
	uint64_t _key = 0;

	/// SplitMix64's finaliser. It avalanches: one input bit changes about half the output bits.
	static uint64_t _mix(uint64_t z) {
		z ^= z >> 30U;
		z *= 0xBF58476D1CE4E5B9ULL;
		z ^= z >> 27U;
		z *= 0x94D049BB133111EBULL;
		z ^= z >> 31U;
		return z;
	}

public:
	/// @param seed       the seed the run was started with
	/// @param stream     which draw of which container this serves
	/// @param iteration  the iteration the draw belongs to
	/// @param dimension  the tree whose node state this is. The field belongs to no tree, so the
	///                   field's streams leave it out.
	TCellUniforms(uint64_t seed, TCellStream stream, size_t iteration, size_t dimension = 0) {
		_key = _mix(seed + _gamma * static_cast<uint64_t>(stream));
		_key = _mix(_key + _gamma * static_cast<uint64_t>(dimension));
		_key = _mix(_key + _gamma * static_cast<uint64_t>(iteration));
	}

	/// The uniform on [0, 1) the cell at `linear_index_in_container_space` draws.
	[[nodiscard]] double at(size_t linear_index_in_container_space) const {
		const uint64_t bits =
		    _mix(_key + _gamma * static_cast<uint64_t>(linear_index_in_container_space));
		// The top 53 bits fill a double's mantissa, which covers [0, 1) in even steps and never
		// reaches 1. A uniform of exactly 1 would make a probability of 1 reject.
		return static_cast<double>(bits >> 11U) * 0x1.0p-53;
	}
};

/// The seed the run was started with.
///
/// Read it outside a parallel region. The generator is thread-local, and a worker thread holds a
/// copy of it whose seed no option fixed. Every stream is therefore built before the region that
/// uses it.
inline uint64_t run_seed() {
	return static_cast<uint64_t>(coretools::instances::randomGenerator().getSeed());
}

#endif // ACOL_TCELLUNIFORMS_H
