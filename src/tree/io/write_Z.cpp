#include "tree/io/write_Z.h"

#include "Types.h"
#include "coretools/Files/TOutputFile.h"
#include "storages/storage_backend.h"
#include "tree/TTree.h"

#include <array>

namespace {

/// The node names of one cell, one per tree. Writing this out once per set of cells is what made
/// the two write paths look like two different things.
std::vector<std::string> node_names_of(const IndexArray &multidim_index,
                                       const std::vector<std::unique_ptr<TTree>> &trees) {
	std::vector<std::string> names;
	names.reserve(multidim_index.size());
	for (size_t idx = 0; idx < multidim_index.size(); ++idx) {
		// Every column holds a node index: this tree's own dimension spans all its nodes, and in
		// the others a leaf's index in leaf space is its node index (ADR-0004).
		names.push_back(trees[idx]->get_node_id(multidim_index[idx]));
	}
	return names;
}

} // namespace

void write_Z_to_file(const std::string &filename, const TTree &tree,
                     const std::vector<std::unique_ptr<TTree>> &trees, bool write_full_Z) {
	std::vector<std::string> header;
	header.reserve(trees.size() + 2);
	for (const auto &t : trees) { header.push_back(t->get_tree_name()); }
	header.emplace_back("position");
	header.emplace_back("Z_state");

	const auto &Z = tree.get_Z();
	coretools::TOutputFile file(filename, header, "\t");

	const auto write_cell = [&](size_t linear_index_in_Z_space, bool state) {
		const std::array<size_t, 2> line{linear_index_in_Z_space, state};
		file.writeln(node_names_of(Z.get_multi_dimensional_index(linear_index_in_Z_space), trees),
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

void write_branch_length_grid(const TTree &tree) {
	const std::vector<std::string> header = {"grid_position", "branch_length"};
	coretools::TOutputFile file(
	    "acol_simulated_" + tree.get_tree_name() + "_branch_length_grid.txt", header, "\t");

	const auto &grid_branch_lengths = tree.grid_branch_lengths();
	for (size_t i = 0; i < grid_branch_lengths.size(); ++i) {
		file.writeln(i, grid_branch_lengths[i]);
	}
}
