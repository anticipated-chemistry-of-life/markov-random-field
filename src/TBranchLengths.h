//
// Created by madleina on 22.10.24.
//

#ifndef ACOL_TBRANCHLENGTHS_H
#define ACOL_TBRANCHLENGTHS_H

#include "Types.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Types/probability.h"
#include "coretools/algorithms.h"
#include "coretools/devtools.h"
#include <algorithm>
#include <armadillo>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>

class TMatrix {
private:
	arma::mat _mat;

public:
	TMatrix() { _mat.zeros(2, 2); }

	void set_from_matrix_exponential(const arma::mat &Lambda) { _mat = arma::expmat(Lambda); }
	arma::mat get_matrix() const { return _mat; }

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
	arma::mat _lambda_c = arma::zeros(2, 2);

public:
	TMatrices() = default;
	explicit TMatrices(size_t NumBins) { resize(NumBins); }

	void resize(size_t NumBins) { _matrices.resize(NumBins); }
	[[nodiscard]] size_t size() const { return _matrices.size(); }
	const std::vector<TMatrix> &get_matrices() const { return _matrices; }

	void set_lambda(double mu_c_1, double mu_c_0) {
		_lambda_c[0] = -mu_c_1;
		_lambda_c[1] = mu_c_0;
		_lambda_c[2] = mu_c_1;
		_lambda_c[3] = -mu_c_0;
	}

	void set(double a, double Delta) {
		// calculate matrix exponential for first bin
		TMatrix P_0;
		P_0.set_from_matrix_exponential(_lambda_c * a);

		// calculate matrix exponential of scaling matrix
		TMatrix _matrix_alpha;
		_matrix_alpha.set_from_matrix_exponential(_lambda_c * Delta);

		_matrices[0] = P_0;
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

	void _discretize_branch_lengths(std::vector<double> BranchLengths) {
		// normalize such that they sum to one
		coretools::normalize(BranchLengths);

		std::vector<double> grid(_matrices.size());
		for (size_t k = 0; k < _matrices.size(); ++k) { grid[k] = (_a + _delta * (k + 1)); }

		_binned_branch_lengths.resize(BranchLengths.size());

		for (size_t i = 0; i < BranchLengths.size(); ++i) {
			// find bin
			auto it = std::lower_bound(grid.begin(), grid.end(), BranchLengths[i]);
			if (it == grid.end()) {
				// last bin
				_binned_branch_lengths[i] = grid.size() - 1;
			} else {
				_binned_branch_lengths[i] = std::distance(grid.begin(), it);
			}
		}
	}

public:
	TBranchLengths() = default;
	[[nodiscard]] double get_a() const { return (double)_a; }
	[[nodiscard]] double get_b() const { return (double)_b; }
	[[nodiscard]] double get_delta() const { return _delta; }
	void compute_matrices() { _matrices.set((double)_a, _delta); }
	void set_lambda(double mu_c_1, double mu_c_0) { _matrices.set_lambda(mu_c_1, mu_c_0); }
	const TMatrices &get_matrices() const { return _matrices; }

	void initialize(const std::vector<double> &branch_length) {
		// read a, b and K from command-line
		_a               = coretools::instances::parameters().get("a", coretools::Probability(0.0));
		double default_b = std::min(1.0, 1.0 / (double)branch_length.size() * 10);
		_b               = coretools::instances::parameters().get("b", coretools::Probability(default_b));
		size_t K         = coretools::instances::parameters().get("K", 100);
		if (K >= (size_t)std::numeric_limits<TypeBinBranches>::max()) {
			UERROR("More bins (", K, ") required than type allows (", std::numeric_limits<TypeBinBranches>::max(),
			       ")! Please decrease K or change type of bins.");
		}

		_matrices.resize(K);

		// calculate Delta
		_delta = ((double)_b - (double)_a) / (double)K;

		_discretize_branch_lengths(branch_length);
	}

	[[nodiscard]] size_t size() const { return _binned_branch_lengths.size(); }

	const std::vector<TypeBinBranches> &get_binned_branch_lengths() const { return _binned_branch_lengths; }

	void update_branch_lengths(size_t ix_1, size_t ix_2) {
		// ...
	}
};

#endif // ACOL_TBRANCHLENGTHS_H
