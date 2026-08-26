#include "tree/io/read_Z.h"

#include "constants.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Main/TError.h"
#include "storages/z_storage/TStorageZMatrix.h"
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

	const size_t node  = topology.index_of(node_name);
	const size_t index = internal_space ? topology.internal_index(node) : topology.leaf_index(node);
	if (index == TPhylogeny::NOT_IN_SPACE) {
		throw coretools::TUserError(
		    "Node '", node_name, "' of tree '", tree.get_tree_name(), "' is ",
		    internal_space ? "a leaf, but this column holds internal nodes."
		                   : "an internal node, but this column holds leaves.");
	}
	return index;
}

} // namespace

void read_Z_from_file(const std::string &filename, TStorageZMatrix &Z,
                      const std::vector<std::unique_ptr<TTree>> &trees,
                      size_t dimension_number_of_tree) {
	coretools::TInputFile file(filename, coretools::FileType::Header);

	// One name column per tree, then the linear index and the state.
	const size_t expected_columns = trees.size() + 2;
	if (file.numCols() != expected_columns) {
		throw coretools::TUserError("The file '", filename, "' for setting Z must have ",
		                            expected_columns, " columns, but has ", file.numCols(), " !");
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
