//
// Shared parsing helpers for the sparse TSV data files (LOTUS and the simple error model).
//
// Both sources use the same on-disk format: a header naming one tree per column, then one row per
// cell whose state is 1, giving that cell's leaf node id in each tree. Cells not listed are 0.
// Both name every tree, in tree order, so both files have exactly one column per tree. The header
// validation and the node-name -> leaf-index resolution are shared and live here.
//
// Not guarded by any USE_* macro: these helpers are needed by whichever data sources are compiled
// in, and by none of them exclusively.
//

#pragma once

#include "coretools/Files/TInputFile.h"
#include "coretools/Main/TError.h"
#include "tree/TTree.h"
#include <algorithm>
#include <cstddef>
#include <deque>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace sparse_data_file {

/// Tree names in tree order -- the header of a sparse data file that names every tree.
[[nodiscard]] inline std::vector<std::string>
header_from_trees(const std::vector<std::unique_ptr<TTree>> &trees) {
	std::vector<std::string> header;
	header.reserve(trees.size());
	for (const auto &tree : trees) { header.push_back(tree->get_tree_name()); }
	return header;
}

/// A file names every tree, exactly once, in tree order: column i is matched positionally against
/// the i-th tree. Throws TUserError otherwise.
///
/// Both sources are indexed on every tree, so the count is pinned here rather than by each caller.
/// That is also what turns the order check below into an exact one: with the counts equal, a header
/// that is a subsequence of the tree names is the tree names.
inline void validate_header_against_trees(const coretools::TInputFile &file,
                                          const std::vector<std::unique_ptr<TTree>> &trees,
                                          std::string_view filename) {
	if (file.header().size() != trees.size()) {
		throw coretools::TUserError("File '", filename, "' has ", file.header().size(),
		                            " columns but there are ", trees.size(),
		                            " trees. This data source is indexed on every tree, so the "
		                            "file must have exactly one column per tree.");
	}

	const auto tree_names = header_from_trees(trees);

	for (const auto &header_name : file.header()) {
		if (std::find(tree_names.begin(), tree_names.end(), header_name) == tree_names.end()) {
			throw coretools::TUserError("Header '", header_name, "' in file '", filename,
			                            "' does not match any tree name!");
		}
	}

	// Check that the headers are a subsequence of the tree names: walk the trees in order and pop
	// each header that matches. Anything left over was out of order.
	std::deque<std::string> remaining_headers(file.header().begin(), file.header().end());
	for (const auto &tree_name : tree_names) {
		if (remaining_headers.empty()) { break; }
		if (remaining_headers.front() == tree_name) { remaining_headers.pop_front(); }
	}
	if (!remaining_headers.empty()) {
		throw coretools::TUserError("Headers in file '", filename,
		                            "' should be ordered in the same way as the trees.");
	}
}

/// Leaf index of `node_name` within `tree`. Throws TUserError if the node is unknown (from
/// TTree::get_node_index) or if it is an internal node: the data files may only reference leaves,
/// since the observed matrices are indexed in leaf space.
[[nodiscard]] inline size_t leaf_index_or_throw(const TTree &tree, const std::string &node_name) {
	if (!tree.isLeaf(tree.get_node_index(node_name))) {
		throw coretools::TUserError("Node '", node_name, "' in tree '", tree.get_tree_name(),
		                            "' is not a leaf !");
	}
	return tree.get_index_within_leaves(node_name);
}

} // namespace sparse_data_file
