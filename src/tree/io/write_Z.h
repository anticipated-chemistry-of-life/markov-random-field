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
/// not only the internal ones (ADR-0005).
///
/// `write_full_Z` picks the cells: the whole container space, a missing cell reading as state 0,
/// or only the stored ones. Simulation, the only caller today, asks for the whole space.
template<BinaryFieldStorage NodeState>
void write_Z_to_file(const std::string &filename, const NodeState &Z,
                     const std::vector<TNodeStateColumn> &columns, bool write_full_Z) {
	throw_unless_one_column_per_tree(columns);

	std::vector<std::string> header;
	header.reserve(columns.size() + 2);
	for (const auto &column : columns) { header.push_back(column.tree_name); }
	header.emplace_back("position");
	header.emplace_back("Z_state");

	coretools::TOutputFile file(filename, header, "\t");

	const auto write_cell = [&](size_t linear_index_in_Z_space, bool state) {
		const std::array<size_t, 2> line{linear_index_in_Z_space, state};
		file.writeln(node_names_of(Z.get_multi_dimensional_index(linear_index_in_Z_space), columns),
		             line);
	};

	// Which cells, and nothing else, is what the two write paths differ in: the whole container
	// space, where a missing cell reads as state 0 and so a point lookup covers both cases, or
	// only the stored entries, which come in ascending linear-index order.
	if (write_full_Z) {
		for (size_t i = 0; i < Z.total_size_of_container_space(); ++i) {
			write_cell(i, Z.is_one(i));
		}
	} else {
		for (const auto &[linear_index_in_Z_space, storage] : Z.get_stored_entries()) {
			write_cell(linear_index_in_Z_space, storage.is_one());
		}
	}
}

/// Write the branch-length grid of `tree`, one row per bin.
///
/// A written node state only ever names the bin a branch fell into, so this is what makes such a
/// file interpretable. It used to be emitted as a side effect of writing the node state; it is
/// its own file, so it is its own function.
void write_branch_length_grid(const TTree &tree);
