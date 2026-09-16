//
// A tree's per-clique transition grids, stored as one contiguous array.
//
// TTransitionGrid (TTransitionGrid.h) is a self-contained value: built from an alpha, a nu and a
// bin grid, and owning its own matrices. That is what keeps it constructible and testable on its
// own (ADR-0003), and what a Metropolis proposal builds its candidate as.
//
// A tree keeps one of those per clique, though, and there can be as many cliques as the other
// tree has leaves -- hundreds of thousands on a real dataset. Storing them as
// `std::vector<std::optional<TTransitionGrid>>` puts every clique's bins in their own heap
// allocation, scattered wherever the allocator happened to put it, so a lookup on the block
// update's hot path (TTree::transition_grid_of_cell, asked twice per leaf pair) chases two
// pointers -- one for the slot, one for its bins -- to reach four doubles. This file trades that
// for one allocation per tree: every clique's bins live in one contiguous buffer, addressed by
// `clique * n_bins + bin`, with one alpha per clique beside it.
//
// `TTransitionGridView` reads that buffer with the exact surface `TTransitionGrid` offers
// (`probability`, `stationary`), so a caller that only reads a grid -- the node-state walk, the
// density, the branch-length ratio -- takes either one through the `TransitionGridLike` concept
// below and does not know which it got. Only `TTree` ever sees `TTransitionGridTable` itself.
//

#ifndef ACOL_TTRANSITIONGRIDTABLE_H
#define ACOL_TTRANSITIONGRIDTABLE_H

#include "coretools/Main/TError.h"
#include "tree/branch/TTransitionGrid.h"

#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <vector>

/// A read-only look at one clique's transition grid, backed by a `TTransitionGridTable`'s
/// contiguous storage rather than owning its own. Cheap to return and to copy: a pointer, a bin
/// count and an alpha, and no allocation.
class TTransitionGridView {
private:
	const std::array<double, 4> *_matrices; ///< this clique's `_n_bins` matrices, contiguous
	size_t _n_bins;
	double _alpha;

public:
	TTransitionGridView(const std::array<double, 4> *matrices, size_t n_bins, double alpha)
	    : _matrices(matrices), _n_bins(n_bins), _alpha(alpha) {}

	/// `P(child = to | parent = from)` across a branch sitting in `bin`. Same contract as
	/// `TTransitionGrid::probability`.
	[[nodiscard]] double probability(size_t bin, bool from, bool to) const {
		return _matrices[bin][static_cast<size_t>(from) * 2 + static_cast<size_t>(to)];
	}

	/// The stationary probability of a state. Same contract as `TTransitionGrid::stationary`.
	[[nodiscard]] double stationary(bool state) const { return state ? _alpha : 1.0 - _alpha; }

	[[nodiscard]] double alpha() const { return _alpha; }
	[[nodiscard]] size_t n_bins() const { return _n_bins; }
};

/// What a caller that only reads a process asks of it: `TTransitionGrid` and `TTransitionGridView`
/// both answer the same two questions, so either satisfies this and the node-state walk, the
/// density and the branch-length ratio do not have to know which one they were handed.
template<typename T>
concept TransitionGridLike = requires(const T &process, size_t bin, bool from, bool to) {
	{ process.probability(bin, from, to) } -> std::same_as<double>;
	{ process.stationary(from) } -> std::same_as<double>;
};

/// One tree's transition grids, one per clique, in one contiguous allocation.
///
/// A slot starts unset -- `resize` fills every alpha with NaN, which `is_set` reads back rather
/// than keeping a second array just to say so. `set` is how a grid gets installed: once a clique
/// when the parameters first exist, and again whenever a proposal on that clique's alpha or nu is
/// accepted. `at` is the read path the block update takes twice per leaf pair, which is what makes
/// the contiguous layout worth it (see the file comment).
class TTransitionGridTable {
private:
	std::vector<std::array<double, 4>> _matrices; ///< clique `c`, bin `b` -> `_matrices[c*_n_bins+b]`
	std::vector<double> _alpha;                    ///< one per clique; NaN marks "not set yet"
	size_t _n_bins = 0;

public:
	/// Sizes the table for `n_cliques` cliques of `n_bins` bins each. Every clique starts unset.
	void resize(size_t n_cliques, size_t n_bins) {
		_n_bins = n_bins;
		_matrices.assign(n_cliques * n_bins, std::array<double, 4>{});
		_alpha.assign(n_cliques, std::numeric_limits<double>::quiet_NaN());
	}

	[[nodiscard]] bool is_set(size_t c) const { return !std::isnan(_alpha[c]); }

	/// Installs `grid` at clique `c`, copying its bins into this table's contiguous storage.
	void set(size_t c, const TTransitionGrid &grid) {
		DEBUG_ASSERT(grid.n_bins() == _n_bins);
		const auto src = grid.matrices();
		std::copy(src.begin(), src.end(), _matrices.begin() + static_cast<ptrdiff_t>(c * _n_bins));
		_alpha[c] = grid.alpha();
	}

	/// Clique `c`'s grid. Throws rather than reading zeros if nothing has installed one yet --
	/// `TTransitionGrid`'s own guarantee (TTransitionGrid.h), carried over to its table.
	[[nodiscard]] TTransitionGridView at(size_t c) const {
		if (!is_set(c)) {
			throw coretools::TDevError("TTransitionGridTable: clique ", c, " has no grid yet.");
		}
		return {&_matrices[c * _n_bins], _n_bins, _alpha[c]};
	}
};

#endif // ACOL_TTRANSITIONGRIDTABLE_H
