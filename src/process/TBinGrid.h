//
// The bin <-> branch-length correspondence for one tree.
//
// Branch lengths are never continuous in this model: they are only ever bins (CONTEXT.md, "Bin").
// This module owns that correspondence and nothing else. It is a pure function of `n_bins`, with no
// knowledge of the tree's topology, of stattools parameters, or of the random generator: callers
// draw, this module consumes. That is what makes it testable without building an MCMC.
//
// The branch-length budget is the invariant tying the pieces together: a tree's bins always sum to
// `n_branches * n_bins / 2`, equivalently the mean grid branch length is exactly 1 (CONTEXT.md,
// "Branch-length budget"). `bins_from_lengths` establishes it and `step_direction` preserves it, so
// both live here rather than one here and one in the caller.
//
// Deliberately *not* here: `TypeBinnedBranchLengths::setMax`, which configures a stattools type and
// so cannot sit behind a pure interface. TTree keeps it as a one-line adapter, and this module
// carries its own `n_bins` instead of reading the type's bounds back out.
//

#ifndef ACOL_TBINGRID_H
#define ACOL_TBINGRID_H

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <vector>

/// The bin grid of one tree: bin width, the length each bin stands for, and the budget their sum
/// must hit.
class TBinGrid {
private:
	size_t _n_bins = 0;
	double _delta  = 0.0;
	std::vector<double> _grid_branch_lengths;

public:
	/// @param n_bins Number of bins on the grid. Must be at least 1.
	explicit TBinGrid(size_t n_bins) : _n_bins(n_bins) {
		if (n_bins == 0) { throw std::invalid_argument("TBinGrid needs at least one bin."); }

		_delta = 2.0 / ((double)n_bins + 1.0);

		// Bin k stands for the branch length at the *center* of the bin: the transition matrices
		// are built as exp(Lambda * Delta * (k + 0.5)), so Delta * (k + 0.5) is the length bin k
		// actually represents in the likelihood.
		_grid_branch_lengths.resize(n_bins);
		for (size_t k = 0; k < n_bins; ++k) {
			_grid_branch_lengths[k] = _delta * ((double)k + 0.5);
		}
	}

	[[nodiscard]] size_t n_bins() const { return _n_bins; }
	[[nodiscard]] double delta() const { return _delta; }
	/// Static because the lowest bin is 0 whatever the grid; `max_bin` is the one that moves.
	[[nodiscard]] static constexpr size_t min_bin() { return 0; }
	[[nodiscard]] size_t max_bin() const { return _n_bins - 1; }

	/// The continuous length each bin stands for, indexed by bin.
	[[nodiscard]] const std::vector<double> &grid_branch_lengths() const {
		return _grid_branch_lengths;
	}

	/// The continuous length bin `bin` stands for: its center, `Delta * (bin + 0.5)`.
	[[nodiscard]] double grid_branch_length(size_t bin) const {
		return _grid_branch_lengths.at(bin);
	}

	/// The bin a continuous length falls in. Bin k covers [k * Delta, (k + 1) * Delta), clamped to
	/// the grid; non-positive lengths land in bin 0.
	[[nodiscard]] size_t bin_from_length(double length) const {
		if (length <= 0.0) { return 0; }
		return std::min(static_cast<size_t>(length / _delta), _n_bins - 1);
	}

	/// The total the bins of a tree with `n_branches` branches must sum to, for life.
	[[nodiscard]] size_t budget(size_t n_branches) const { return (n_branches * _n_bins) / 2; }

	[[nodiscard]] bool is_on_budget(const std::vector<size_t> &bins) const {
		return std::accumulate(bins.begin(), bins.end(), size_t{0}) == budget(bins.size());
	}

	/// Walk `bins` onto the budget by +-1 steps.
	///
	/// `pick(n)` must return an index in [0, n). It is called once per iteration, including the
	/// iterations that pick a bin already at a boundary and make no progress, so a caller wiring a
	/// random generator in reproduces exactly the draw sequence this loop consumes.
	template<typename PickIndex>
	void repair_to_budget(std::vector<size_t> &bins, PickIndex pick) const {
		const size_t goal = budget(bins.size());
		size_t total      = std::accumulate(bins.begin(), bins.end(), size_t{0});

		while (total != goal) {
			const size_t ix = pick(bins.size());
			if (total < goal) {
				if (bins[ix] >= _n_bins - 1) { continue; } // already at the top, redraw
				bins[ix] += 1;
				total += 1;
			} else {
				if (bins[ix] == 0) { continue; } // already at the bottom, redraw
				bins[ix] -= 1;
				total -= 1;
			}
		}
	}

	/// Normalise `lengths` to mean 1, bin them, and repair onto the budget.
	///
	/// Every length must be strictly positive. In tree files that holds by construction: the reader
	/// rejects non-positive branch lengths, and the only zeros it writes are the roots' -- which
	/// have no branch and are filtered out by the caller, since knowing which node is a root is
	/// topology and does not belong here.
	template<typename PickIndex>
	[[nodiscard]] std::vector<size_t> bins_from_lengths(const std::vector<double> &lengths,
	                                                    PickIndex pick) const {
		if (lengths.empty()) {
			throw std::invalid_argument("TBinGrid was given no branch lengths.");
		}

		double sum = 0.0;
		for (const double length : lengths) {
			if (length <= 0.0) {
				throw std::invalid_argument(
				    "TBinGrid was given a non-positive branch length. Roots carry no branch and "
				    "must be filtered out before binning.");
			}
			sum += length;
		}
		const double average = sum / (double)lengths.size();

		std::vector<size_t> bins;
		bins.reserve(lengths.size());
		for (const double length : lengths) { bins.push_back(bin_from_length(length / average)); }

		repair_to_budget(bins, pick);
		return bins;
	}

	/// Whether a +-1 step on this pair has a free choice of direction.
	///
	/// Exposed so a caller can draw its coin only when the coin matters: at a boundary the
	/// direction is forced, and drawing anyway would consume a random number the pre-extraction
	/// code did not.
	[[nodiscard]] bool step_is_free(size_t first, size_t second) const {
		return first != min_bin() && second != min_bin() && first != max_bin() &&
		       second != max_bin();
	}

	/// The direction of a budget-conserving +-1 step on a pair of bins: `first` moves by the
	/// returned value and `second` by its negation, so the total is unchanged.
	///
	/// Returns 0 when both bins sit at the same boundary and no legal move exists. `decrease`
	/// picks the direction in the free case and is ignored otherwise; see `step_is_free`.
	[[nodiscard]] int step_direction(size_t first, size_t second, bool decrease) const {
		const size_t lowest  = min_bin();
		const size_t highest = max_bin();

		// Order matters: a pair with one bin at the minimum and the other at the maximum is legal
		// in both directions, and the minimum is checked first.
		if ((first == lowest && second == lowest) || (first == highest && second == highest)) {
			return 0; // both at the same boundary, no update possible
		}
		if (first == lowest) { return 1; }   // can only go right
		if (second == lowest) { return -1; } // second can only go right, so first goes left
		if (first == highest) { return -1; } // can only go left
		if (second == highest) { return 1; } // second can only go left, so first goes right
		return decrease ? -1 : 1;
	}
};

#endif // ACOL_TBINGRID_H
