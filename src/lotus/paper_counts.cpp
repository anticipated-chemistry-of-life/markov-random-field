#include "lotus/paper_counts.h"

#include "coretools/Files/TInputFile.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TParameters.h"
#include "tree/TPhylogeny.h"

std::vector<size_t> read_paper_counts(const std::string &tree_name, const TPhylogeny &topology) {
	const std::string parameter_name = tree_name + "_paper_counts";
	if (!coretools::instances::parameters().exists(parameter_name)) {
		throw coretools::TUserError("Parameter '", parameter_name,
		                            "' not found. Please provide it.");
	}

	const auto filename = coretools::instances::parameters().get<std::string>(parameter_name);
	coretools::TInputFile file(filename, coretools::FileType::Header);

	if (file.numCols() != 2) {
		throw coretools::TUserError("File '", filename, "' is expected to have 2 columns, but has ",
		                            file.numCols(), " !");
	}

	// A leaf the file never mentions keeps its zero.
	std::vector<size_t> paper_counts(topology.n_leaves(), 0);

	for (; !file.empty(); file.popFront()) {
		const std::string leaf_name = std::string(file.get(0));
		const auto count            = file.get<size_t>(1);

		const size_t node_index = topology.index_of(leaf_name);
		if (!topology.is_leaf(node_index)) {
			throw coretools::TUserError("All nodes should be leaves.");
		}
		const size_t leaf_index = topology.leaf_index(node_index);

		if (leaf_index >= paper_counts.size()) {
			throw coretools::TUserError("Leaf index ", leaf_index,
			                            " is out of bounds for paper counts vector of size ",
			                            paper_counts.size(), ".");
		}

		paper_counts[leaf_index] = count;
	}

	return paper_counts;
}
