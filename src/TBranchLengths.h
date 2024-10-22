//
// Created by madleina on 22.10.24.
//

#ifndef ACOL_TBRANCHLENGTHS_H
#define ACOL_TBRANCHLENGTHS_H

#include "Types.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Types/probability.h"
#include "coretools/algorithms.h"
#include <armadillo>

class TMatrix {
private:
	arma::mat _mat;

public:
	TMatrix() { _mat.zeros(2, 2); }

	void set_from_matrix_exponential(const arma::mat &Lambda) { _mat = arma::expmat(Lambda); }

	void set_from_product(const TMatrix &First, const TMatrix &Second) {
		// do matrix multiplication
		// _mat(0,0) = First(0,0) * Second(0,0) + First(0,1) * Second(1,0)
		_mat[0] = First[0] * Second[0] + First[2] * Second[1];
		_mat[2] = 1.0 - _mat[0];
		// TODO: complete
	}

	constexpr double operator()(size_t i, size_t j) const {
		// get element from matrix
		// 0  2
		// 1  3
		return _mat[i + 2 * j];
	}

	double operator[](size_t i) const { return _mat[i]; }
};

class TMatrices {
private:
	std::vector<TMatrix> _matrices;
	TMatrix _matrix_alpha;
	arma::mat _lambda_c;

public:
	TMatrices() = default;
	TMatrices(size_t NumBins) { resize(NumBins); }

	void resize(size_t NumBins) { _matrices.resize(NumBins); }

	void set(double a, double Delta) {
		// calculate matrix exponential for first bin
		// P0 = ... TODO complete

		// calculate matrix exponential of scaling matrix
		// TODO complete

		// do recursion
		for (size_t k = 1; k < _matrices.size(); ++k) {
			_matrices[k].set_from_product(_matrices[k - 1], _matrix_alpha);
		}
	}
};

class TBranchLengths {
private:
	coretools::Probability _a;
	coretools::Probability _b;
	double _delta;

	std::vector<TypeBinBranches> _binned_branch_lengths;
	TMatrices _matrices;

public:
	TBranchLengths() {
		// read a, b and K from command-line
		_a       = coretools::instances::parameters().get("a", coretools::Probability(0.0));
		_b       = coretools::instances::parameters().get("b", coretools::Probability(0.1));
		size_t K = coretools::instances::parameters().get("K", 100);
		if (K >= (size_t)std::numeric_limits<TypeBinBranches>::max()) {
			UERROR("More bins (", K, ") required than type allows (", std::numeric_limits<TypeBinBranches>::max(),
			       ")! Please decrease K or change type of bins.");
		}

		_matrices.resize(K);

		// calculate Delta
		_delta = ((double)_b - (double)_a) / (double)K;
	}

	void set_branch_lengths(std::vector<double> BranchLengths) {
		// normalize such that they sum to one
		coretools::normalize(BranchLengths);

		// TODO: create the grid (min + Delta, min + 2*Delta, ..., max)
		// TODO: assign each branch length to its bin based on the grid
	}

	void update_branch_lengths(size_t ix_1, size_t ix_2) {
		// ...
	}
};

#endif // ACOL_TBRANCHLENGTHS_H
