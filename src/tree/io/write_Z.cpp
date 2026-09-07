#include "tree/io/write_Z.h"

#include "coretools/Files/TOutputFile.h"
#include "tree/TTree.h"

#include <cstddef>
#include <string>
#include <vector>

void write_branch_length_grid(const TTree &tree) {
	const std::vector<std::string> header = {"grid_position", "branch_length"};
	coretools::TOutputFile file(
	    "acol_simulated_" + tree.get_tree_name() + "_branch_length_grid.txt", header, "\t");

	const auto &grid_branch_lengths = tree.grid_branch_lengths();
	for (size_t i = 0; i < grid_branch_lengths.size(); ++i) {
		file.writeln(i, grid_branch_lengths[i]);
	}
}
