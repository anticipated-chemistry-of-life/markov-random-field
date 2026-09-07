#include "tree/io/node_state_columns.h"

#include "constants.h"
#include "coretools/Main/TError.h"
#include "tree/TPhylogeny.h"
#include "tree/TTree.h"

#include <cstddef>

std::vector<TNodeStateColumn> node_state_columns(const std::vector<std::unique_ptr<TTree>> &trees) {
	std::vector<TNodeStateColumn> columns;
	columns.reserve(trees.size());
	for (const auto &tree : trees) {
		columns.push_back(TNodeStateColumn{&tree->phylogeny(), tree->get_tree_name()});
	}
	return columns;
}

void throw_unless_one_column_per_tree(const std::vector<TNodeStateColumn> &columns) {
	if (columns.size() != NUMBER_OF_TREES) {
		throw coretools::TDevError("A node-state file was given ", columns.size(),
		                           " columns, but a cell index holds ", NUMBER_OF_TREES, ".");
	}
}

std::vector<std::string> node_names_of(const IndexArray &multidim_index,
                                       const std::vector<TNodeStateColumn> &columns) {
	std::vector<std::string> names;
	names.reserve(multidim_index.size());
	// One name per coordinate, so the loop is bounded by the cell and not by the column list.
	// throw_unless_one_column_per_tree is what says the two have the same length.
	for (size_t idx = 0; idx < multidim_index.size(); ++idx) {
		names.push_back(columns[idx].topology->id_of(multidim_index[idx]));
	}
	return names;
}
