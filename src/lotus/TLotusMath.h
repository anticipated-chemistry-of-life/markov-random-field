//
// The likelihood of one LOTUS cell.
//
// A LOTUS record is an occurrence someone reported in the literature, so its probability is not a
// flat error rate: a present pair is reported only if anyone looked, and how hard anyone looked is
// the research effort spent on that pair (CONTEXT.md, "Research effort", "LOTUS record").
//
// Research effort is a product over the trees of `1 - exp(-gamma_i * log(count_i + 1))`,
// where `count_i` is the paper count of the leaf that cell occupies in dimension `i` and `gamma_i`
// is that dimension's detection rate. Note there is one gamma *per tree*, not one
// overall. Given the latent state `x` of the cell, the emission is then
//
//     P(L = 1 | x = 1) = research effort          P(L = 1 | x = 0) = error rate
//     P(L = 0 | x = 1) = 1 - research effort      P(L = 0 | x = 0) = 1 - error rate
//
// so an absent pair is reported at a flat rate while a present one is reported at a rate that
// depends on where it sits. That asymmetry is the whole point of the source: a missing record for a
// well-studied pair is evidence of absence, and for an unstudied pair it is evidence of nothing.
//
// The transform `log(count + 1)` lives here rather than in the reader, because CONTEXT.md makes
// the raw publication count the input to research effort. Applying it while reading the file gave
// a reader that returned something other than paper counts, which the validation harness records
// as a trap: simulating LOTUS data from raw counts is indistinguishable from an inference bug.
//
// A grid is immutable. It is a value derived from (gammas, error rate, paper counts), so a
// Metropolis proposal on gamma builds a second one and the caller keeps whichever it accepts; there
// is no refresh hook and nothing to revert.
//
// Everything here is free-standing: it is deliberately NOT guarded by USE_LOTUS, so it compiles
// (and is unit-tested) in every build configuration.
//

#pragma once

#include "constants.h"
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace lotus_math {

/// The paper count as research effort consumes it: `log(count + 1)`.
///
/// A leaf with no papers gives 0, which drives its factor to `1 - exp(0) = 0`, so a pair nobody
/// has studied is never reported. That is the intended behaviour, not a degenerate case.
inline double occurrence_count(size_t paper_count) {
	return std::log(static_cast<double>(paper_count) + 1.0);
}

/// One clique of the observation model: the probability of a LOTUS record, per cell, given the
/// latent state of that cell.
class TReportingModel {
private:
	/// `_factor[dimension][leaf]` caches `1 - exp(-gamma_i * occurrence_count(count))`, so the
	/// per-cell path in the field update is a product of lookups with no exp() and no parameter access.
	std::vector<std::vector<double>> _factor;
	double _error_rate = 0.0;

public:
	/// @param gammas       One detection rate per tree.
	/// @param error_rate   The rate at which an absent pair is reported anyway.
	/// @param paper_counts Raw publication counts, `paper_counts[dimension][leaf]`.
	TReportingModel(const std::vector<double> &gammas, double error_rate,
	                const std::vector<std::vector<size_t>> &paper_counts)
	    : _error_rate(error_rate) {
		if (gammas.size() != paper_counts.size()) {
			throw std::invalid_argument(
			    "TReportingModel needs one gamma per dimension of paper counts.");
		}

		_factor.resize(gammas.size());
		for (size_t dimension = 0; dimension < gammas.size(); ++dimension) {
			const auto &counts = paper_counts[dimension];
			auto &factor       = _factor[dimension];
			factor.resize(counts.size());
			for (size_t leaf = 0; leaf < counts.size(); ++leaf) {
				factor[leaf] = 1.0 - std::exp(-gammas[dimension] * occurrence_count(counts[leaf]));
			}
		}
	}

	[[nodiscard]] size_t n_dimensions() const { return _factor.size(); }
	[[nodiscard]] double error_rate() const { return _error_rate; }

	/// The probability that a *present* pair at `cell` is reported. `cell` is indexed in leaf
	/// space: one leaf index per tree.
	[[nodiscard]] double research_effort(const IndexArray &cell) const {
		double product = 1.0;
		for (size_t dimension = 0; dimension < _factor.size(); ++dimension) {
			product *= _factor[dimension][cell[dimension]];
		}
		return product;
	}

	/// `P(L | x)` for the cell at `cell`.
	[[nodiscard]] double probability(bool x, bool L, const IndexArray &cell) const {
		if (!x) { return probability_absent(L); }
		const double effort = research_effort(cell);
		return L ? effort : 1.0 - effort;
	}

	/// `P(L | x = 0)`, which does not depend on where the cell is.
	///
	/// Separate from `probability` so a caller walking a sparse container can answer the absent
	/// case without building a cell index for it. Most cells are absent, so that conversion is the
	/// bulk of the work it would otherwise do.
	[[nodiscard]] double probability_absent(bool L) const {
		return L ? _error_rate : 1.0 - _error_rate;
	}
};

} // namespace lotus_math
