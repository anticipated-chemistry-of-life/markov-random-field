//
// Writing one tree field's posterior to file.
//
// The field writes its own posterior, and this writes the two that stand behind it. Three files
// rather than three columns of one, because they are three variables: what the tree field says and
// what the field says differ by the link, and a chain riding the ADR-0005 ridge separates them
// (ADR-0005, derivation 3).
//
// It takes a posterior, a node state and one column per tree, so a test builds what it writes --
// the same seam write_Z_to_file takes.
//

#pragma once

#include "coretools/Files/TOutputFile.h"
#include "field/tree_field_posterior.h"
#include "storages/storage_concepts.h"
#include "tree/io/node_state_columns.h"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

/// Write one tree field's posterior: one row per leaf pair, naming the two leaves, the state the
/// chain ended at, and the fraction of counted iterations the cell was a one.
///
/// The rows are the field's own leaf-pair space, so `position` is the field's linear index and not
/// the node state's. A tree field and the field are addressed at the same `(row, column)`
/// (ADR-0005), so `position` is the same name in this file and in the field's posterior. Join the
/// two on it, and never by row number: each file drops its own uninformative cells, and with a
/// non-zero error probability those two sets differ.
///
/// A cell that is not a one now and was never counted a one carries no posterior, and is left out
/// -- the rule the field's posterior file follows, so that both files say the same thing under
/// either backend.
template<FieldStorage Field, BinaryStorage NodeState>
void write_tree_field_posterior(const std::string &filename, const Field &Y, const NodeState &Z,
                                const TTreeFieldPosterior &posterior,
                                const std::vector<TNodeStateColumn> &columns) {
	throw_unless_one_column_per_tree(columns);

	std::vector<std::string> header;
	header.reserve(columns.size() + 3);
	header.emplace_back("position");
	header.emplace_back("Z_state");
	for (const auto &column : columns) { header.push_back(column.tree_name); }
	header.emplace_back("fraction_of_one");

	coretools::TOutputFile file(filename, header, "\t");

	for (size_t cell = 0; cell < posterior.size(); ++cell) {
		const IndexArray leaf_pair = Y.get_multi_dimensional_index(cell);
		const bool state           = Z.is_one(Z.get_linear_index_in_container_space(leaf_pair));
		if (!state && posterior.get_counter(cell) == 0) { continue; }
		const std::array<size_t, 2> line{cell, state};
		file.writeln(line, node_names_of(leaf_pair, columns), posterior.get_fraction_of_ones(cell));
	}
}
