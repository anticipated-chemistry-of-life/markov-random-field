//
// The cells of one clique of one tree's node state, addressed by node index.
//

#pragma once

#include "constants.h"
#include "coretools/Main/TError.h"
#include "storages/cell_write.h"
#include "storages/storage_backend.h"
#include "storages/storage_concepts.h"
#include "tree/TPhylogeny.h"

#include <cstddef>
#include <cstdint>
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

/// The cells of one clique, as the run of cells a node state holds them in.
///
/// There is one way in. A read is a point lookup; a write is the one branch in
/// storages/cell_write.h -- the node state takes it in place, or the cell waits in the deferred
/// list and the caller commits that list once the parallel region ends. Which of the two a
/// backend does is its own business, and neither the walk above nor this class asks.
///
/// A write the node state could not take is kept here as well, so that a later read on this run
/// sees it. That is what a post-order walk needs: it reaches a parent after its children and
/// reads the states they were just given. Without it the sparse and the dense backend would
/// compute different chains inside one update, which is what the parity gate exists to prevent.
///
/// A run writes each of its cells at most once, so a kept write is never written over. That is
/// the same rule `write_or_defer` states for the list it appends to.
template<typename Storage> class TCliqueCells {
private:
	Storage *_Z          = nullptr;
	size_t _start_linear = 0;
	size_t _n_cells      = 0;
	size_t _stride       = 1;

	/// The cells the node state does not hold and a write turned into ones. The caller commits
	/// them once the parallel region ends. See ADR-0006.
	std::vector<size_t> _deferred_inserts;

	/// The same cells, by position in the run, so that a read finds them. Empty until the first
	/// deferred write: a storage that holds every cell of its container space -- which every dense
	/// one does -- defers nothing and allocates nothing here.
	///
	/// A write of zero needs no entry. It is deferred by nobody, because a cell the node state
	/// does not hold already reads as zero.
	std::vector<uint8_t> _kept_ones;

public:
	TCliqueCells(Storage &Z, const IndexArray &first_cell, size_t n_cells, size_t stride)
	    : _Z(&Z), _start_linear(Z.get_linear_index_in_container_space(first_cell)),
	      _n_cells(n_cells), _stride(stride) {
		// The run is described here and read one cell at a time below, so a wrong start or stride
		// would otherwise be caught by the storage, a layer down from where it was made.
		DEBUG_ASSERT(n_cells == 0 ||
		             _start_linear + (n_cells - 1) * stride < Z.total_size_of_container_space());
	}

	[[nodiscard]] size_t size() const { return _n_cells; }

	[[nodiscard]] size_t linear_index(size_t k) const {
		DEBUG_ASSERT(k < _n_cells);
		return _start_linear + k * _stride;
	}

	[[nodiscard]] bool is_one(size_t k) const {
		DEBUG_ASSERT(k < _n_cells);
		if (!_kept_ones.empty() && _kept_ones[k] != 0) { return true; }
		return _Z->is_one(linear_index(k));
	}

	void set_state(size_t k, bool state) {
		if (write_or_defer(*_Z, linear_index(k), state, _deferred_inserts)) { return; }
		// The node state could not take this write, so the run keeps it. Undoing a kept write
		// would have to undo the deferred entry beside it, which nothing here can do without a
		// search -- so the once-per-cell rule above is asserted rather than defended. The window
		// this replaces kept a state per cell and could afford to.
		DEBUG_ASSERT(_kept_ones.empty() || _kept_ones[k] == 0);
		if (!state) { return; }
		if (_kept_ones.empty()) { _kept_ones.assign(_n_cells, 0); }
		_kept_ones[k] = 1;
	}

	[[nodiscard]] std::vector<size_t> take_deferred_inserts() {
		_kept_ones.clear();
		return std::exchange(_deferred_inserts, {});
	}
};

/// Which nodes the walk that holds a view writes.
///
/// A tree's update writes every node, leaves included: a tree field is the leaf block of that
/// tree's node state, and each tree draws its own (ADR-0005). So does a simulation's forward draw.
/// The chain start is the one walk that stops above the leaves -- it starts each internal node at
/// the state its children make most likely, and the leaves it reads were written before it ran. The
/// view is the only place that can catch a write to the wrong block, so it is told which of the two
/// walks holds it.
enum class TCliqueWrites : uint8_t { every_node, internal_nodes };

/// The cells of one clique of one tree's node state, addressed by node index.
///
/// A tree's update walks nodes; a node state is indexed by cell. The view is the one place that
/// turns the first into the second. It holds the node state, the clique's multidimensional index
/// and the tree's own dimension, and maps a node to a cell by setting that dimension to the node.
///
/// It reads one tree alone. A tree field is the leaf block of this very run of cells (ADR-0005),
/// so the node-state walk, the alpha and nu moves after it, the branch-length likelihood and the
/// density pass all address this and name neither the field nor the other tree. The walk draws its
/// leaves from a link its caller bound, and reaches that link through a seam rather than through
/// this view.
///
/// It satisfies the clique-column concept the node-state draw and the node-state density are
/// written against. Those two headers reach a real node state through this. They reach a vector of
/// states through a test's own type. Neither header knows the difference.
///
/// The view keeps no copy of the states the node state took. It reads those back from the cells
/// it wrote them through, which is what lets a post-order walk read the states it has just given.
/// A cell the node state could not take is the one exception, and TCliqueCells holds it.
template<typename Storage, TCliqueWrites Writes = TCliqueWrites::every_node> class TCliqueView {
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
		// The only place left that can catch a write to the wrong block. Which block that is
		// depends on the walk, which is what `Writes` says.
		if constexpr (Writes == TCliqueWrites::internal_nodes) {
			DEBUG_ASSERT(!_topology->is_leaf(node));
		}
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

/// The view a tree's update opens over one clique of its own node state, and the same view a
/// simulation's forward draw opens over those cells. Both write every node: a tree field is the
/// leaf block of the node state, and each tree draws its own.
using TNodeStateCliqueView = TCliqueView<TNodeStateStorage>;

/// The view the chain start opens over the same cells. It differs in one thing: the start writes
/// internal nodes only, because it starts each of them at the state its children make most likely
/// and the leaves below it already carry the states it reads.
using TNodeStateStartView = TCliqueView<TNodeStateStorage, TCliqueWrites::internal_nodes>;
