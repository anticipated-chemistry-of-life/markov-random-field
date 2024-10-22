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
#include <vector>

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
		_mat[1] = First[1] * Second[0] + First[3] * Second[1];
		_mat[3] = 1.0 - _mat[1];
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
	explicit TMatrices(size_t NumBins) { resize(NumBins); }

	void resize(size_t NumBins) { _matrices.resize(NumBins); }
	[[nodiscard]] size_t size() const { return _matrices.size(); }

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

	void discretize_branch_lengths(std::vector<double> BranchLengths) {
		if (_matrices.size() > BranchLengths.size()) {
			UERROR("Number of bins (", _matrices.size(), ") is larger than number of branch lengths (",
			       BranchLengths.size(), ")! Please decrease K or provide a bigger graph.");
		}
		// normalize such that they sum to one
		coretools::normalize(BranchLengths);

		std::vector<double> grid;
		double upper_bound = _a;
		while (upper_bound < _b) {
			grid.push_back(upper_bound);
			upper_bound += _delta;
		}
		grid.push_back(_b);
		grid.push_back(std::numeric_limits<double>::max());

		// assign each branch length to its bin
		_binned_branch_lengths.resize(BranchLengths.size());

		for (size_t i = 0; i < BranchLengths.size(); ++i) {
			// find bin
			auto it                   = std::lower_bound(grid.begin(), grid.end(), BranchLengths[i]);
			size_t bin                = std::distance(grid.begin(), it);
			_binned_branch_lengths[i] = bin;
		}
	}

	void update_branch_lengths(size_t ix_1, size_t ix_2) {
		// ...
	}
};

#endif // ACOL_TBRANCHLENGTHS_H
