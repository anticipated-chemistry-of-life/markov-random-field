//
// The cells of one clique of one tree's node state, addressed by node index.
//

#pragma once

#include "constants.h"
#include "coretools/Main/TError.h"
#include "storages/cell_handle.h"
#include "storages/storage_backend.h"
#include "storages/storage_concepts.h"
#include "tree/TPhylogeny.h"

#include <cstddef>
#include <utility>
#include <vector>

/// The cell node `node` of a clique occupies: the clique's own multidimensional index with the
/// tree's dimension set to that node.
///
/// This is the clique's whole node-to-cell arithmetic, and it is written here once. A clique's
/// index carries a leaf in every dimension but its tree's own, which it carries a 0 in; setting
/// that dimension to a node index is what walks the clique. The rule is the same whichever
/// dimension the tree owns, so there is one function and not one per dimension.
[[nodiscard]] inline IndexArray clique_cell(const IndexArray &clique_index, size_t dimension,
                                            size_t node) {
	IndexArray cell = clique_index;
	cell[dimension] = node;
	return cell;
}

/// The distance in linear index between two consecutive nodes of a clique.
///
/// A clique's cells are a strided run of the node state -- the tree's dimension moves and every
/// other one stays -- so one subtraction names the run, whichever dimension the tree owns and
/// whatever the container's shape. A clique of a single node has no second cell, and any stride
/// describes a run of one.
template<BinaryStorage Storage>
[[nodiscard]] size_t clique_stride(const Storage &Z, const IndexArray &clique_index,
                                   size_t dimension, size_t n_nodes) {
	if (n_nodes < 2) { return 1; }
	return Z.get_linear_index_in_container_space(clique_cell(clique_index, dimension, 1)) -
	       Z.get_linear_index_in_container_space(clique_cell(clique_index, dimension, 0));
}

/// How a clique's cells are reached, which is the one thing the two node-state backends still do
/// differently.
///
/// A storage that answers `locate` gives up one cell at a time. The handle says whether the
/// storage holds that cell. An in-place write and a deferred insert are then one branch apart
/// (cell_handle.h).
///
/// A storage that answers no handle keeps its window. The sorted-vector matrix is the one such
/// storage left. It holds every cell twice, once in its row and once in its column, so it has no
/// single cell to point at. ADR-0006 argues the window it keeps.
///
/// Both spellings answer the same four questions. TCliqueView below never says which one it
/// holds.
template<typename Storage, bool = LocatableStorage<Storage>> class TCliqueCells;

/// The handle path. A read is a point lookup and a write is one branch.
template<typename Storage> class TCliqueCells<Storage, true> {
private:
	Storage *_Z          = nullptr;
	size_t _start_linear = 0;
	size_t _n_cells      = 0;
	size_t _stride       = 1;

	/// The cells the storage does not hold and a write turned into ones. The caller commits them
	/// once the parallel region ends. See ADR-0006.
	std::vector<size_t> _deferred_inserts;

public:
	TCliqueCells(Storage &Z, const IndexArray &first_cell, size_t n_cells, size_t stride)
	    : _Z(&Z), _start_linear(Z.get_linear_index_in_container_space(first_cell)),
	      _n_cells(n_cells), _stride(stride) {}

	[[nodiscard]] size_t size() const { return _n_cells; }

	[[nodiscard]] size_t linear_index(size_t k) const {
		DEBUG_ASSERT(k < _n_cells);
		return _start_linear + k * _stride;
	}

	[[nodiscard]] bool is_one(size_t k) const { return _Z->is_one(linear_index(k)); }

	void set_state(size_t k, bool state) {
		const size_t deferred_before = _deferred_inserts.size();
		write_or_defer(_Z->locate(linear_index(k)), state, _deferred_inserts);
		// A deferred cell still reads as zero here, because this read goes to the storage. A
		// post-order walk reads the state it has just given a child, so a storage that defers has
		// to bring that read-back with it. Every storage that answers a handle today holds every
		// cell of its container space, so nothing is deferred on this path. See ADR-0006.
		DEBUG_ASSERT(_deferred_inserts.size() == deferred_before);
	}

	[[nodiscard]] std::vector<size_t> take_deferred_inserts() {
		return std::exchange(_deferred_inserts, {});
	}
};

/// The window path, for a storage that has no cell to point at yet. The window walks the clique's
/// line once on open and buffers the writes the container cannot take, which is the same two
/// answers the handle gives one cell at a time.
template<typename Storage> class TCliqueCells<Storage, false> {
private:
	typename Storage::TWindow _window;

public:
	TCliqueCells(Storage &Z, const IndexArray &first_cell, size_t n_cells, size_t stride)
	    : _window(Z.open_window(first_cell, n_cells, stride)) {}

	[[nodiscard]] size_t size() const { return _window.size(); }
	[[nodiscard]] size_t linear_index(size_t k) const { return _window.linear_index(k); }
	[[nodiscard]] bool is_one(size_t k) const { return _window.is_one(k); }
	void set_state(size_t k, bool state) { _window.set_state(k, state); }

	[[nodiscard]] std::vector<size_t> take_deferred_inserts() {
		return _window.take_buffered_inserts();
	}
};

/// The cells of one clique of one tree's node state, addressed by node index.
///
/// A tree's update walks nodes; a node state is indexed by cell. The view is the one place that
/// turns the first into the second. It holds the node state, the clique's multidimensional index
/// and the tree's own dimension, and maps a node to a cell by setting that dimension to the node.
///
/// It reads one tree alone. A tree field is the leaf block of this very run of cells (ADR-0005),
/// so the node-state walk, the alpha and nu moves after it, the branch-length likelihood and the
/// density pass all address this and name neither the field nor the other tree.
///
/// It satisfies the clique-column concept the node-state draw and the node-state density are
/// written against. Those two headers reach a real node state through this. They reach a vector of
/// states through a test's own type. Neither header knows the difference.
///
/// The view keeps no copy of the states it assigns. It reads them back from the cells it wrote
/// them through, which is what lets a post-order walk read the states it has just given. A cell
/// the node state could not take is read back by the window that holds it.
template<typename Storage> class TCliqueView {
private:
	const TPhylogeny *_topology = nullptr;
	IndexArray _clique_index{};
	size_t _dimension = 0;
	TCliqueCells<Storage> _cells;

public:
	/// Opens a view over the cells of `clique_index`, the clique of the tree that owns
	/// `dimension`. The clique's index carries a leaf in every other dimension.
	TCliqueView(Storage &Z, const TPhylogeny &topology, const IndexArray &clique_index,
	            size_t dimension)
	    : _topology(&topology), _clique_index(clique_index), _dimension(dimension),
	      _cells(Z, clique_cell(clique_index, dimension, 0), topology.n_nodes(),
		         clique_stride(Z, clique_index, dimension, topology.n_nodes())) {}

	// A view may own writes that are not in the node state yet, so it is neither copied nor moved.
	// That is true of one backend only, and both are held to it, so that code written against
	// either compiles against the other.
	TCliqueView(const TCliqueView &)            = delete;
	TCliqueView &operator=(const TCliqueView &) = delete;
	TCliqueView(TCliqueView &&)                 = delete;
	TCliqueView &operator=(TCliqueView &&)      = delete;
	~TCliqueView()                              = default;

	/// The number of cells, which is every node of the tree -- leaves included.
	[[nodiscard]] size_t size() const { return _cells.size(); }

	/// The cell node `node` occupies, as a multidimensional index.
	[[nodiscard]] IndexArray cell_of(size_t node) const {
		return clique_cell(_clique_index, _dimension, node);
	}

	/// The state of `node`, from this tree's own node state.
	[[nodiscard]] bool is_one(size_t node) const { return _cells.is_one(node); }

	/// The linear index, in the node state's container space, of the cell `node` occupies. This is
	/// the cell's name to the stream of uniforms (ADR-0007).
	[[nodiscard]] size_t linear_index(size_t node) const { return _cells.linear_index(node); }

	/// Gives `node` its new state.
	void set_state(size_t node, bool state) {
		// The walk assigns internal nodes. A leaf's state is drawn with the field and the other
		// tree's leaf, as one eight-state block, and not here (ADR-0005). This is the only place
		// left that can catch a write to the wrong block, and it inverts when the walk covers
		// leaves.
		DEBUG_ASSERT(!_topology->is_leaf(node));
		// A write of the state the cell already carries is dropped. A sparse node state would
		// otherwise defer an insert for a cell it does not hold, and then throw that insert away.
		if (_cells.is_one(node) != state) { _cells.set_state(node, state); }
	}

	/// Hands out the cells the node state could not take, as linear indices in its container
	/// space. The caller commits them after the parallel region. That is the only exit a view
	/// inside one may take. ADR-0006 gives the argument.
	[[nodiscard]] std::vector<size_t> take_deferred_inserts() {
		return _cells.take_deferred_inserts();
	}
};

/// The view a tree opens over one clique of its own node state.
using TNodeStateCliqueView = TCliqueView<TNodeStateStorage>;
