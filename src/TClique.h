//
// Created by madleina on 22.10.24.
//

#ifndef ACOL_TBRANCHLENGTHS_H
#define ACOL_TBRANCHLENGTHS_H

#include "TStorageYVector.h"
#include "TStorageZVector.h"
#include <armadillo>
#include <cstddef>
#include <vector>
class TTree; // forward declaration
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

class TClique {
private:
	TMatrices _matrices;
	double _mu_c_1;
	double _mu_c_0;
	std::vector<size_t> _start_index;
	size_t _variable_dimension;
	size_t _n_nodes;

public:
	TClique() = default;
	TClique(const std::vector<size_t> &start_index, size_t variable_dimension, size_t n_nodes) {
		_start_index        = start_index;
		_variable_dimension = variable_dimension;
		_n_nodes            = n_nodes;
	};
	void initialize(double a, double delta, size_t n_bins) {
		_matrices.resize(n_bins);
		_matrices.set(a, delta);
	}
	void set_lambda() { _matrices.set_lambda(_mu_c_1, _mu_c_0); }
	const TMatrices &get_matrices() const { return _matrices; }

	void update_Z(const TStorageYVector &Y, const TStorageZVector &Z, const TTree &tree);
};

#endif // ACOL_TBRANCHLENGTHS_H
