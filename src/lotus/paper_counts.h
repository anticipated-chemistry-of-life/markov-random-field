//
// Per-leaf publication counts, the constant half of LOTUS' research effort.
//
// This sits beside the reporting model that consumes it rather than on the tree: its only
// connection to a tree is that it resolves a leaf name to a leaf index, and it was on the tree
// only because the tree class was where everything to do with tree files ended up.
//

#pragma once

#include <cstddef>
#include <string>
#include <vector>

class TPhylogeny;

/// Raw publication counts per leaf, in leaf order, exactly as the file states them.
///
/// The file is named by the `<tree_name>_paper_counts` command-line parameter. Leaves the file
/// does not mention keep a count of zero.
///
/// The `log(count + 1)` that research effort consumes is applied by `lotus_math::TReportingModel`,
/// not here: CONTEXT.md makes the raw count the input to research effort, and a reader that
/// quietly returned logs was a documented trap in the validation harness.
///
/// Throws coretools::TUserError if the parameter is missing, the file does not have two columns,
/// a named node is not a leaf, or the leaf index it resolves to is out of range.
[[nodiscard]] std::vector<size_t> read_paper_counts(const std::string &tree_name,
                                                    const TPhylogeny &topology);
