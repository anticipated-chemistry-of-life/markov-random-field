//
// Writing a tree's node state to file.
//
// This is output, not model: it belongs beside the tree's other file handling rather than on the
// tree class, where no test could reach it. It takes a node state and one column per tree, so a
// test builds what it writes.
//

#pragma once

#include "coretools/Files/TOutputFile.h"
#include "storages/storage_concepts.h"
#include "tree/io/node_state_columns.h"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

class TTree;

/// Write every node-state cell of `Z` as one row: one node name per column, then the linear index,
/// then the state.
///
/// The node state spans every node of its own tree, so a written file carries the leaf rows and
/// not only the internal ones (ADR-0005). Every cell of the container space gets a row, a missing
/// cell reading as state 0, so the file reads the same under either backend.
template<BinaryStorage NodeState>
void write_Z_to_file(const std::string &filename, const NodeState &Z,
                     const std::vector<TNodeStateColumn> &columns) {
	throw_unless_one_column_per_tree(columns);

	std::vector<std::string> header;
	header.reserve(columns.size() + 2);
	for (const auto &column : columns) { header.push_back(column.tree_name); }
	header.emplace_back("position");
	header.emplace_back("Z_state");

	coretools::TOutputFile file(filename, header, "\t");

	for (size_t i = 0; i < Z.total_size_of_container_space(); ++i) {
		const std::array<size_t, 2> line{i, Z.is_one(i)};
		file.writeln(node_names_of(Z.get_multi_dimensional_index(i), columns), line);
	}
}

/// Write the branch-length grid of `tree`, one row per bin.
///
/// A written node state only ever names the bin a branch fell into, so this is what makes such a
/// file interpretable. It used to be emitted as a side effect of writing the node state; it is
/// its own file, so it is its own function.
void write_branch_length_grid(const TTree &tree);
