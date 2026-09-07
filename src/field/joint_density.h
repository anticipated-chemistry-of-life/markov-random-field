//
// The joint density, by factor.
//
// ADR-0005 factors the model as
//
//     p(Z_s | theta_s) * p(Z_m | theta_m) * p(Y | Z_s, Z_m, omega) * p(L, D | Y)
//
// and every factor is a proper conditional density, so their product is one too. That is what makes
// this number the no-drift instrument: it is a density that moves only with the configuration and
// the parameters, where the old sum of the two trees' likelihoods was not a density at all
// (ADR-0002).
//
// The factors are kept apart rather than added up on the spot, because the trace writes one column
// each. A single total says a chain drifts; the columns say which factor is dragging it.
//

#pragma once

#include "constants.h"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace joint_density {

/// One iteration's joint density, factor by factor. One entry of `node_states` per tree, in
/// tree order.
struct TJointDensity {
	std::array<double, NUMBER_OF_TREES> node_states{};
	/// `log p(Y | Z_s, Z_m, omega)`, from the six link counters.
	double link = 0.0;
	/// `log p(L, D | Y)`: the LOTUS records and the simple error model, whichever of them is
	/// compiled in. Zero for a simulated chain, which draws from the prior and has no data to
	/// score.
	double data = 0.0;

	[[nodiscard]] double total() const {
		double sum = link + data;
		for (const double tree : node_states) { sum += tree; }
		return sum;
	}
};

/// The columns of the trace file, in file order: one per tree, then the link, the data and the
/// total. `tree_names` is in tree order and names the `node_states` entries.
[[nodiscard]] inline std::vector<std::string>
trace_header(const std::vector<std::string> &tree_names) {
	std::vector<std::string> header;
	header.reserve(tree_names.size() + 3);
	for (const auto &name : tree_names) { header.push_back(name + "_node_state"); }
	header.emplace_back("link");
	header.emplace_back("data");
	header.emplace_back("joint_density");
	return header;
}

/// One row of the trace file, in the order `trace_header` names.
[[nodiscard]] inline std::vector<double> trace_row(const TJointDensity &density) {
	std::vector<double> row;
	row.reserve(density.node_states.size() + 3);
	for (const double tree : density.node_states) { row.push_back(tree); }
	row.push_back(density.link);
	row.push_back(density.data);
	row.push_back(density.total());
	return row;
}

} // namespace joint_density
