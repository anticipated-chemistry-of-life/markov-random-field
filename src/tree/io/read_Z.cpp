#include "tree/io/read_Z.h"

#include "constants.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Main/TError.h"
#include "tree/TPhylogeny.h"
#include "tree/io/node_state_columns.h"

namespace {

/// One column's node name, as the index that column's dimension is addressed by.
size_t index_in_dimension(const TNodeStateColumn &column, const std::string &node_name,
                          bool own_column) {
	const TPhylogeny &topology = *column.topology;
	if (!topology.contains(node_name)) {
		throw coretools::TUserError("Node '", node_name, "' does not exist in tree '",
		                            column.tree_name, "'.");
	}

	const size_t node = topology.index_of(node_name);
	// The column belonging to this node state's own tree accepts any node of it, leaves included:
	// the node state now spans them. Every other column indexes a leaf, so it still has to be one.
	if (!own_column && !topology.is_leaf(node)) {
		throw coretools::TUserError("Node '", node_name, "' of tree '", column.tree_name,
		                            "' is an internal node, but this column holds leaves.");
	}
	// Either way the answer is the node index: for a leaf, its index in leaf space is its node
	// index (ADR-0004), and for this tree's own column the node state is indexed by node.
	return node;
}

} // namespace

void read_Z_cells_from_file(const std::string &filename,
                            const std::vector<TNodeStateColumn> &columns,
                            size_t dimension_number_of_tree,
                            const std::function<void(const IndexArray &)> &insert_the_cell) {
	coretools::TInputFile file(filename, coretools::FileType::Header);

	// One name column per tree, then the linear index and the state.
	const size_t expected_columns = columns.size() + 2;
	if (file.numCols() != expected_columns) {
		throw coretools::TUserError("The file '", filename, "' for setting Z must have ",
		                            expected_columns, " columns, but has ", file.numCols(), " !");
	}

	throw_unless_one_column_per_tree(columns);

	const size_t state_column = columns.size() + 1;
	for (; !file.empty(); file.popFront()) {
		// An absent cell reads as state 0, so a row saying 0 says nothing this loop has to act on.
		if (!file.get<bool>(state_column)) { continue; }

		IndexArray multidim_index{};
		for (size_t idx = 0; idx < columns.size(); ++idx) {
			multidim_index[idx] = index_in_dimension(columns[idx], std::string(file.get(idx)),
			                                         idx == dimension_number_of_tree);
		}
		insert_the_cell(multidim_index);
	}
}
