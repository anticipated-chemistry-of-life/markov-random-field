//
// Writing a tree's internal state to file.
//
// This is output, not model: it belongs beside the tree's other file handling rather than on the
// class, where it was a template in a header that no test could reach.
//

#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

class TTree;

/// Write every internal-state cell of `tree` as one row: one node name per tree, then the linear
/// index, then the state.
///
/// `dimension_number_of_tree` says which of `trees` this file belongs to; that dimension is
/// indexed in internal-node space and every other one in leaf space, which is what decides how a
/// cell's coordinates turn back into node names.
///
/// `write_full_Z` picks the cells: the whole container space, a missing cell reading as state 0,
/// or only the stored ones. Simulation, the only caller today, asks for the whole space.
void write_Z_to_file(const std::string &filename, const TTree &tree,
                     const std::vector<std::unique_ptr<TTree>> &trees,
                     size_t dimension_number_of_tree, bool write_full_Z);

/// Write the branch-length grid of `tree`, one row per bin.
///
/// A written internal state only ever names the bin a branch fell into, so this is what makes such
/// a file interpretable. It used to be emitted as a side effect of writing internal state; it is
/// its own file, so it is its own function.
void write_branch_length_grid(const TTree &tree);
