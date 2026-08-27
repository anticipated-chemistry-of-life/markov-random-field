//
// Reading a previously written internal state back in.
//

#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

class TStorageZMatrix;
class TTree;

/// Read an internal-state file written by write_Z_to_file into `Z`, setting every cell the file
/// reports as present.
///
/// Cells are resolved through the leading node-name columns, not through the linear index the
/// writer also emits: `dimension_number_of_tree` says which column names an internal node and
/// which name leaves, and the rest is name lookup. That makes such a file portable across any
/// change to how nodes are numbered, and turns what would otherwise be a silent wrong answer --
/// an index that still resolves, but to a different node -- into an error.
///
/// Throws coretools::TUserError if the file does not have one column per tree plus the index and
/// the state, if a name is not in the tree its column belongs to, or if a name is in the tree but
/// in the wrong space: a leaf where the file wants an internal node, or the other way round.
void read_Z_from_file(const std::string &filename, TStorageZMatrix &Z,
                      const std::vector<std::unique_ptr<TTree>> &trees,
                      size_t dimension_number_of_tree);
