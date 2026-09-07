//
// The columns a node-state file names its cells through.
//

#pragma once

#include "constants.h"

#include <memory>
#include <string>
#include <vector>

class TPhylogeny;
class TTree;

/// One column of a node-state file: the tree whose node names that column holds, and that tree's
/// name.
///
/// The reader and the writer want a topology and a name, and nothing else a tree carries. Taking
/// the pair rather than the tree is what lets both of them run in a test. A phylogeny is a value.
/// A tree stands up three live stattools parameters and a command line.
struct TNodeStateColumn {
	const TPhylogeny *topology;
	std::string tree_name;
};

/// The columns of a node-state file over these trees, in file order.
///
/// One column per tree, in the order the trees are in. The reader and the writer both assume that
/// order, and a node state's own dimension is the column at its own tree's index.
[[nodiscard]] std::vector<TNodeStateColumn>
node_state_columns(const std::vector<std::unique_ptr<TTree>> &trees);

/// Throws coretools::TDevError unless there is exactly one column per tree.
///
/// A cell index is a fixed-size IndexArray, so more columns than it holds would run a fill loop
/// past its end. The reader and the writer both assume one column per dimension. This is a guard
/// on that assumption, and not a limit either of them imposes.
void throw_unless_one_column_per_tree(const std::vector<TNodeStateColumn> &columns);

/// The node names of one cell, one name per column.
///
/// Every coordinate of a cell is a node index -- a node state's own dimension spans every node of
/// its tree, and in the others a leaf's index in leaf space is its node index (ADR-0004) -- so
/// this needs no knowledge of which column belongs to which tree.
[[nodiscard]] std::vector<std::string> node_names_of(const IndexArray &multidim_index,
                                                     const std::vector<TNodeStateColumn> &columns);
