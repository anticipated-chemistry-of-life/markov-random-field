//
// The indexed set of one tree's cliques: a clique number and the cell subscript it stands for.
//
// A clique of a tree carries a leaf in every dimension but that tree's own, which carries a 0
// (CONTEXT.md, "Clique"). This module owns the correspondence between a clique's number and that
// subscript, in both directions, and nothing else. It holds no storage, no topology, no parameter
// and no generator, so a test builds one from two numbers.
//
// One convention serves both directions: row-major, through coretools' `getSubscriptsAsArray` and
// `getLinearIndex`, which are exact inverses. Every storage in the codebase linearises the same way
// (storages/TDense.h, storages/TSparse.h), so a clique subscript and a cell subscript are read
// alike. ADR-0011 records why this convention and not the other.
//
// `NumDim` is a test affordance and not a seam. Production instantiates it once, at
// `NUMBER_OF_TREES`. The dimension count is a template parameter so that a test can build a clique
// space over three dimensions, which is the case that tells the two conventions apart, and which
// the shipping configuration cannot produce. Nothing about the model becomes N-tree by it: the
// block update, the eight-state table and `IndexArray` all still name two trees.
//

#pragma once

#include "constants.h"
#include "coretools/algorithms.h"

#include <array>
#include <cstddef>
#include <stdexcept>

/// The cliques of one tree, numbered.
///
/// @tparam NumDim The number of tree dimensions. Production takes the default; see the file header.
template<size_t NumDim = NUMBER_OF_TREES> class TCliqueSpace {
private:
	/// The extent of clique space in each dimension: the leaf count of every other tree, and a 1
	/// in this tree's own. The 1 is what puts a 0 in the own dimension of every subscript, and what
	/// makes the product below count the cliques.
	std::array<size_t, NumDim> _dimension_cliques{};

	/// The dimension this tree owns.
	size_t _dimension = 0;

public:
	/// @param leaf_counts The number of leaves of every tree, this one included.
	/// @param dimension   The dimension this tree owns. Must name one of `NumDim`.
	TCliqueSpace(const std::array<size_t, NumDim> &leaf_counts, size_t dimension)
	    : _dimension_cliques(leaf_counts), _dimension(dimension) {
		if (dimension >= NumDim) {
			throw std::invalid_argument("A tree of a clique space over " + std::to_string(NumDim) +
			                            " dimensions cannot own dimension " +
			                            std::to_string(dimension) + ".");
		}
		// The tree's own dimension does not vary along a clique, so clique space has extent 1
		// there. Establishing it here is what lets both directions below assume it.
		_dimension_cliques[_dimension] = 1;
	}

	/// The dimension the owning tree occupies.
	[[nodiscard]] size_t dimension() const { return _dimension; }

	/// How many cliques the tree has: the product of every other tree's leaf count.
	[[nodiscard]] size_t n_cliques() const {
		return coretools::containerProduct(_dimension_cliques);
	}

	/// The cell subscript of clique `clique`: a leaf in every dimension but the owning tree's,
	/// which carries a 0.
	[[nodiscard]] std::array<size_t, NumDim> index_of(size_t clique) const {
		return coretools::getSubscriptsAsArray(clique, _dimension_cliques);
	}

	/// The clique the cell at `cell` belongs to.
	///
	/// The owning dimension is zeroed before the conversion. Clique space has extent 1 there, so a
	/// caller that passes a real node or leaf index in that slot would otherwise leave the
	/// coordinate out of range. Zeroing here is what lets every caller hand over the whole cell and
	/// name no slot.
	[[nodiscard]] size_t clique_of(const std::array<size_t, NumDim> &cell) const {
		std::array<size_t, NumDim> along_the_clique = cell;
		along_the_clique[_dimension]                = 0;
		return coretools::getLinearIndex(along_the_clique, _dimension_cliques);
	}
};
