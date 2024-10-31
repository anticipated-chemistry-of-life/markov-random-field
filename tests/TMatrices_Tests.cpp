#include "TClique.h"
#include "TTree.h"
#include "armadillo"
#include "coretools/Main/TParameters.h"
#include "coretools/devtools.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <cstddef>
#include <string>
#include <vector>
using namespace testing;

TEST(TMatrices_tests, test_matrix) {
	coretools::instances::parameters().add("K", 10);
	TTree tree(0);
	tree.load_from_file("../tests/test_data/binning.tsv");
	TClique clique;
	clique.set_lambda();
	clique.initialize(tree.get_a(), tree.get_delta(), tree.get_number_of_bins());
	branchLengths.compute_matrices();
	auto p1           = clique.get_matrices().get_matrices()[1].get_matrix();
	arma::mat p1_true = arma::zeros(2, 2);
	p1_true[0]        = 0.970737654700428;
	p1_true[1]        = 0.019508230199714;
	p1_true[2]        = 0.029262345299572;
	p1_true[3]        = 0.980491769800286;
	EXPECT_TRUE(arma::approx_equal(p1, p1_true, "absdiff", 1e-14));

	auto p2           = clique.get_matrices().get_matrices()[2].get_matrix();
	arma::mat p2_true = arma::zeros(2, 2);
	p2_true[0]        = 0.942902450821576;
	p2_true[1]        = 0.038065032785616;
	p2_true[2]        = 0.057097549178424;
	p2_true[3]        = 0.961934967214384;
	EXPECT_TRUE(arma::approx_equal(p2, p2_true, "absdiff", 1e-14));
}
