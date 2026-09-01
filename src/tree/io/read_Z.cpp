#include "tree/io/read_Z.h"

#include "constants.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Main/TError.h"
#include "storages/storage_backend.h"
#include "tree/TPhylogeny.h"
#include "tree/TTree.h"

namespace {

/// One column's node name, as an index into the space that column is written in.
size_t index_in_dimension(const TTree &tree, const std::string &node_name, bool internal_space) {
	const auto &topology = tree.phylogeny();
	if (!topology.contains(node_name)) {
		throw coretools::TUserError("Node '", node_name, "' does not exist in tree '",
		                            tree.get_tree_name(), "'.");
	}

	const size_t node = topology.index_of(node_name);
	// The column belonging to this node state's own tree accepts any node of it, leaves included:
	// the node state now spans them. Every other column indexes a leaf, so it still has to be one.
	if (!internal_space && !topology.is_leaf(node)) {
		throw coretools::TUserError("Node '", node_name, "' of tree '", tree.get_tree_name(),
		                            "' is an internal node, but this column holds leaves.");
	}
	// Either way the answer is the node index: for a leaf, its index in leaf space is its node
	// index (ADR-0004), and for this tree's own column the node state is indexed by node.
	return node;
}

} // namespace

void read_Z_from_file(const std::string &filename, TNodeStateStorage &Z,
                      const std::vector<std::unique_ptr<TTree>> &trees,
                      size_t dimension_number_of_tree) {
	coretools::TInputFile file(filename, coretools::FileType::Header);

	// One name column per tree, then the linear index and the state.
	const size_t expected_columns = trees.size() + 2;
	if (file.numCols() != expected_columns) {
		throw coretools::TUserError("The file '", filename, "' for setting Z must have ",
		                            expected_columns, " columns, but has ", file.numCols(), " !");
	}

	// A cell index is a fixed-size IndexArray, so more trees than it holds would run the fill loop
	// below past its end. The rest of the update makes the same assumption -- this is a guard on it,
	// not a limit this reader imposes.
	if (trees.size() != NUMBER_OF_TREES) {
		throw coretools::TDevError("read_Z_from_file was given ", trees.size(),
		                           " trees, but a cell index holds ", NUMBER_OF_TREES, ".");
	}

	const size_t state_column = trees.size() + 1;
	for (; !file.empty(); file.popFront()) {
		// An absent cell reads as state 0, so a row saying 0 says nothing this loop has to act on.
		if (!file.get<bool>(state_column)) { continue; }

		IndexArray multidim_index{};
		for (size_t idx = 0; idx < trees.size(); ++idx) {
			multidim_index[idx] = index_in_dimension(*trees[idx], std::string(file.get(idx)),
			                                         idx == dimension_number_of_tree);
		}
		Z.insert_one(multidim_index);
	}
}
