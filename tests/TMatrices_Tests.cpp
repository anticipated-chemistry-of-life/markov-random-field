#include "TBranchLengths.h"
#include "TTree.h"
#include "coretools/Main/TParameters.h"
#include "coretools/devtools.h"
#include "gtest/gtest.h"
#include <cstddef>
#include <string>
#include <vector>

TEST(TMatrices_tests, test_matrix) {
	coretools::instances::parameters().add("K", 10);
	TTree tree;
	tree.load_from_file("../tests/test_data/binning.tsv");
	TBranchLengths branchLengths;
	branchLengths.initialize(tree.get_all_branch_lengths());
	branchLengths.set_lambda(0.3, 0.2);
	branchLengths.compute_matrices();

	// TODO : finish test
}
