//
// The joint density, as the trace file carries it.
//
// The arithmetic is a sum, and the tests below say so. What is worth pinning is that the header
// and the row stay in the same order: a mismatch there mislabels a whole trace and nothing in a
// run would say so.
//

#include "field/joint_density.h"
#include "gtest/gtest.h"

#include <string>
#include <vector>

namespace {

using joint_density::TJointDensity;

TEST(JointDensity, adds_both_trees_the_link_and_the_data) {
	const TJointDensity density{.node_states = {-3.0, -5.0}, .link = -0.5, .data = -1.25};
	EXPECT_DOUBLE_EQ(density.total(), -9.75);
}

TEST(JointDensity, writes_one_column_per_tree_then_the_link_the_data_and_the_total) {
	const std::vector<std::string> tree_names = {"species", "molecules"};
	const auto header                         = joint_density::trace_header(tree_names);
	EXPECT_EQ(header, (std::vector<std::string>{"species_node_state", "molecules_node_state",
	                                            "link", "data", "joint_density"}));
}

TEST(JointDensity, puts_every_value_in_the_column_the_header_names) {
	const TJointDensity density{.node_states = {-3.0, -5.0}, .link = -0.5, .data = -1.25};
	const auto row = joint_density::trace_row(density);

	ASSERT_EQ(row.size(), joint_density::trace_header({"species", "molecules"}).size());
	EXPECT_EQ(row, (std::vector<double>{-3.0, -5.0, -0.5, -1.25, -9.75}));
}

} // namespace
