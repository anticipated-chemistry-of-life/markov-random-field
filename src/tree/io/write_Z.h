//
// Writing a tree's node state to file.
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

/// Write every node-state cell of `tree` as one row: one node name per tree, then the linear
/// index, then the state.
///
/// Every column of a cell's coordinates is a node index -- this tree's own dimension spans all its
/// nodes, and in the others a leaf's index in leaf space is its node index (ADR-0004) -- so turning
/// coordinates back into node names needs no knowledge of which column belongs to which tree.
///
/// `write_full_Z` picks the cells: the whole container space, a missing cell reading as state 0,
/// or only the stored ones. Simulation, the only caller today, asks for the whole space.
void write_Z_to_file(const std::string &filename, const TTree &tree,
                     const std::vector<std::unique_ptr<TTree>> &trees, bool write_full_Z);

/// Write the branch-length grid of `tree`, one row per bin.
///
/// A written node state only ever names the bin a branch fell into, so this is what makes such a
/// file interpretable. It used to be emitted as a side effect of writing the node state; it is
/// its own file, so it is its own function.
void write_branch_length_grid(const TTree &tree);
