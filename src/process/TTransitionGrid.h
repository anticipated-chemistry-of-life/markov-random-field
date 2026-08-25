//
// The two-state process of one clique, discretised onto a bin grid.
//
// One clique carries one alpha and one nu (CONTEXT.md, "Alpha", "Nu"), which together define a
// two-state continuous-time Markov chain on the tree's branches:
//
//     Lambda = [[  -alpha * nu   ,     alpha * nu   ],
//               [ (1-alpha) * nu , (alpha - 1) * nu ]]
//
// Branches are only ever bins, so what the likelihood actually needs is one transition matrix per
// bin. This module holds exactly those, plus the stationary distribution the roots are drawn from.
//
// The grid is built by *recursion*, not by evaluating each bin's own branch length: the closed form
// is used exactly twice, for `P_0 = expm(Lambda * Delta / 2)` and the step `expm(Lambda * Delta)`,
// and then `P_k = P_{k-1} * step`. That divergence from the Python reference, which does evaluate
// per bin, is deliberate and load-bearing; see docs/adr/0003 before simplifying it away.
//
// A grid is immutable: it is a value derived from (alpha, nu, bins) and nothing else. A Metropolis
// proposal builds a second one and the caller keeps whichever it accepts, which is why there is no
// "try" state, no accept hook, and no `UseTry` template parameter anywhere in the interface.
//

#ifndef ACOL_TTRANSITIONGRID_H
#define ACOL_TTRANSITIONGRID_H

#include "TBinGrid.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

/// One clique's transition matrices, one per bin, plus its stationary distribution.
class TTransitionGrid {
private:
	/// Above this rate the grid short-circuits to exact (1-alpha, alpha) rows. Threshold from
	/// https://github.com/anticipated-chemistry-of-life/RMF-manuscript/blob/34ac59ccf4b6e70135dd22d13a8604d690965725/R-scripts/when_to_use_stationary_for_P.R
	///
	/// It was derived at a branch length of 1 -- the *mean* grid branch length -- where exp(-nu *
	/// t) is indeed below 1e-10. It is applied at every bin regardless, and bin 0 stands for Delta
	/// / 2 ~ 0.09, so at nu just past the threshold the shortest branch is off by ~6e-2 rather than
	/// 1e-10. Reaching 1e-10 at bin 0 would need nu > 248.
	///
	/// This does not weaken neutrality (ADR-0001): a pinned dimension gets exactly (0.5, 0.5) here
	/// by construction, so it is uninformative regardless of what the true process would have said.
	/// It does mean a species clique whose nu wanders past 25 sees its short branches jump.
	/// `TTransitionGrid_Tests` pins the asymmetry.
	static constexpr double STATIONARY_NU_THRESHOLD = 25.0;

	/// `_matrices[bin][from * 2 + to]` is P(child = to | parent = from) across that bin's branch.
	std::vector<std::array<double, 4>> _matrices;
	double _alpha = 0.0;

	/// `expm(Lambda * t)` in closed form. The generator has eigenvalues 0 and -nu, so its
	/// exponential collapses to a single decay term and needs no numerical routine.
	static std::array<double, 4> transition(double alpha, double nu, double t) {
		const double decay = std::exp(-nu * t);
		return {1.0 - alpha + alpha * decay, alpha * (1.0 - decay), (1.0 - alpha) * (1.0 - decay),
		        alpha + (1.0 - alpha) * decay};
	}

	/// Rows sum to 1, so only the first column is multiplied out and the second is the remainder.
	static std::array<double, 4> product(const std::array<double, 4> &first,
	                                     const std::array<double, 4> &second) {
		std::array<double, 4> out{};
		out[0] = first[0] * second[0] + first[1] * second[2];
		out[1] = 1.0 - out[0];
		out[2] = first[2] * second[0] + first[3] * second[2];
		out[3] = 1.0 - out[2];
		return out;
	}

public:
	TTransitionGrid(double alpha, double nu, const TBinGrid &bins) : _alpha(alpha) {
		_matrices.resize(bins.n_bins());

		if (nu > STATIONARY_NU_THRESHOLD) {
			for (auto &matrix : _matrices) { matrix = {1.0 - alpha, alpha, 1.0 - alpha, alpha}; }
			return;
		}

		const double delta = bins.delta();
		// The first bin is half a step in, since bin k stands for the length at its centre.
		_matrices[0]       = transition(alpha, nu, delta / 2.0);
		// Step matrix is also called matrix alpha in the supplementary notes.
		const auto step    = transition(alpha, nu, delta);
		// Only the two base matrices are evaluated directly; the rest of the grid is walked by
		// repeated multiplication. See docs/adr/0003 for why this is not a loop over `transition`.
		for (size_t k = 1; k < _matrices.size(); ++k) {
			_matrices[k] = product(_matrices[k - 1], step);
		}
	}

	/// `P(child = to | parent = from)` across a branch sitting in `bin`.
	[[nodiscard]] double probability(size_t bin, bool from, bool to) const {
		return _matrices[bin][static_cast<size_t>(from) * 2 + static_cast<size_t>(to)];
	}

	/// The stationary probability of a state, which is what a root is drawn from.
	[[nodiscard]] double stationary(bool state) const { return state ? _alpha : 1.0 - _alpha; }

	[[nodiscard]] double alpha() const { return _alpha; }
	[[nodiscard]] size_t n_bins() const { return _matrices.size(); }
};

#endif // ACOL_TTRANSITIONGRID_H
