//
// Reading a previously written node state back in.
//

#pragma once

#include "constants.h"
#include "storages/storage_concepts.h"
#include "tree/io/node_state_columns.h"

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

/// Walk a node-state file and hand every cell it reports as present to `insert_the_cell`, as the
/// coordinates of that cell.
///
/// This is the whole of the reader that does not touch a storage: the file shape, the resolution
/// of a row to a cell, and every error either of them can raise. It streams, so a file with more
/// rows than fit in memory is still one row at a time.
///
/// Cells are resolved through the leading node-name columns, not through the linear index the
/// writer also emits, so such a file is portable across any change to how nodes are numbered and a
/// stale index becomes an error rather than a silent wrong answer.
///
/// `dimension_number_of_tree` says which column belongs to the tree this node state is for. That
/// column accepts any node of its tree, leaves included, since the node state spans them all
/// (ADR-0005); every other column indexes a leaf and still has to name one.
///
/// Throws coretools::TUserError if the file does not have one column per tree plus the index and
/// the state, if a name is not in the tree its column belongs to, or if a foreign column names
/// something other than a leaf.
void read_Z_cells_from_file(const std::string &filename,
                            const std::vector<TNodeStateColumn> &columns,
                            size_t dimension_number_of_tree,
                            const std::function<void(const IndexArray &)> &insert_the_cell);

/// Read a node-state file written by write_Z_to_file into `Z`, setting every cell the file reports
/// as present. See read_Z_cells_from_file for how a row becomes a cell, and for what is rejected.
///
/// The leaf rows are read like every other row, and the chain start keeps them: a supplied node
/// state is the run's to decide, tree field included (ADR-0019). A file written under the old model
/// names only internal nodes in its own column. It still loads, because resolution is by name, and
/// it leaves every leaf row at zero. Under --fix_Z that zero is held; otherwise the first block
/// update draws the leaf layer.
template<BinaryStorage NodeState>
void read_Z_from_file(const std::string &filename, NodeState &Z,
                      const std::vector<TNodeStateColumn> &columns,
                      size_t dimension_number_of_tree) {
	read_Z_cells_from_file(filename, columns, dimension_number_of_tree,
	                       [&Z](const IndexArray &cell) {
		                       Z.insert_one(Z.get_linear_index_in_container_space(cell));
	                       });
}
