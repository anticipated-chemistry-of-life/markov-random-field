//
// Reading a previously written node state back in.
//

#pragma once

#include "storages/storage_backend.h"
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

class TTree;

/// Read a node-state file written by write_Z_to_file into `Z`, setting every cell the file reports
/// as present.
///
/// Cells are resolved through the leading node-name columns, not through the linear index the
/// writer also emits, so such a file is portable across any change to how nodes are numbered and a
/// stale index becomes an error rather than a silent wrong answer.
///
/// `dimension_number_of_tree` says which column belongs to the tree this node state is for. That
/// column accepts any node of its tree, leaves included, since the node state spans them all
/// (ADR-0005); every other column indexes a leaf and still has to name one.
///
/// A file written under the old model names only internal nodes in its own column. It still loads,
/// because resolution is by name -- but it leaves every leaf row at zero, which is the migration
/// ADR-0004's closing consequence describes and #40 has to decide about.
///
/// Throws coretools::TUserError if the file does not have one column per tree plus the index and
/// the state, if a name is not in the tree its column belongs to, or if a foreign column names
/// something other than a leaf.
void read_Z_from_file(const std::string &filename, TInternalStateStorage &Z,
                      const std::vector<std::unique_ptr<TTree>> &trees,
                      size_t dimension_number_of_tree);
