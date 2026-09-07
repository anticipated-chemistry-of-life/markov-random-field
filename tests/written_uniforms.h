//
// The uniform each cell draws, written out rather than hashed.
//
// `TCellUniforms` has a suite of its own. What every other suite needs is to name the number a
// given cell gets, so that a draw can be asserted against a threshold instead of against a
// frequency. This satisfies the `CellUniforms` concept and nothing else.
//

#pragma once

#include <cstddef>
#include <random>
#include <vector>

namespace uniforms {

/// A stream of uniforms a test writes, addressed by linear index.
class TWrittenUniforms {
private:
	std::vector<double> _values;

public:
	explicit TWrittenUniforms(size_t size, double value = 0.0) : _values(size, value) {}

	[[nodiscard]] double at(size_t linear_index) const { return _values.at(linear_index); }
	void set(size_t linear_index, double value) { _values.at(linear_index) = value; }

	/// Fills the whole stream, for a suite that asserts a frequency rather than one draw.
	void fill_from(std::mt19937_64 &rng) {
		std::uniform_real_distribution<double> uniform(0.0, 1.0);
		for (double &value : _values) { value = uniform(rng); }
	}
};

} // namespace uniforms
