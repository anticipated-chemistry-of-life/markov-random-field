#include "TInput.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/devtools.h"
#include <cstddef>
#include <vector>

void TLinks::load_from_file(const std::string &filename) {
	coretools::instances::logfile().listFlush("Reading links from file '", filename, "' ...");
	coretools::TInputFile file(filename, coretools::FileType::Header);

	if (file.numCols() != this->_trees.size()) {
		UERROR("File '", filename, "' is expected to have", this->_trees.size(), "columns, but has ", file.numCols(),
		       " !");
	}

	for (; !file.empty(); file.popFront()) {
		std::vector<size_t> coordinate(this->_trees.size());
		for (size_t i = 0; i < this->_trees.size(); ++i) {
			std::string node = std::string(file.get(i));
			if (!this->in_any_tree(node)) { UERROR("Node '", node, "' doesn't exist in any of the provided trees !"); }
			if (!this->is_leaf_in_any_tree(node)) {
				UERROR("Node '", node, "' is not a leaf ! So far, we have defined our model to only accept leaves.");
			}

			// we get the index within the leaves of the tree
			coordinate[this->_tree_index_of_node(node)] =
			    this->_trees[this->_tree_index_of_node(node)].get_index_within_leaves(
			        this->_trees[this->_tree_index_of_node(node)].get_node_index(node));
		}

		// we currently have the multidimensional index in the Y space. We need to convert it to a linear index in Y
		// space and store it in the data_Y vector.
		OUT(coordinate);
		size_t linear_index_in_Y_space = this->_data_Y.get_linear_index_in_Y_space(coordinate);
		OUT(linear_index_in_Y_space);
		this->_data_Y.insert_one(linear_index_in_Y_space);
	}
};
