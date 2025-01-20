//
// Created by madleina on 22.10.24.
//

#ifndef ACOL_TBRANCHLENGTHS_H
#define ACOL_TBRANCHLENGTHS_H

#include "TCurrentState.h"
#include "TStorageYVector.h"
#include "TStorageZ.h"
#include "TStorageZVector.h"
#include "TTree.h"
#include "coretools/Math/TAcceptOddsRation.h"
#include "coretools/Math/TSumLog.h"
#include "coretools/devtools.h"
#include <armadillo>
#include <cstddef>
#include <tuple>
#include <vector>

/**
 * @brief Class to store the matrix exponential of the scaling matrix and the matrix exponential of the rate matrix for
 * each bin.
 */
class TMatrix {
private:
	arma::mat _mat;

public:
	TMatrix() { _mat.zeros(2, 2); }

	/** @brief Set _mat to the matrix exponential of the provided matrix.
	 * @param Lambda Matrix to calculate the matrix exponential from.
	 */
	void set_from_matrix_exponential(const arma::mat &Lambda) { _mat = arma::expmat(Lambda); }

	/** @brief Get the matrix.
	 * @return Matrix.
	 */
	[[nodiscard]] arma::mat get_matrix() const { return _mat; }

	/** @brief Perform matrix multiplication of two matrices. Since the matrix rows have to sum to 1, the second value
	 * of the row can easily be calculted.
	 * @param First
	 * @param Second
	 */
	void set_from_product(const TMatrix &First, const TMatrix &Second) {
		// do matrix multiplication
		// _mat(0,0) = First(0,0) * Second(0,0) + First(0,1) * Second(1,0)
		_mat[0] = First[0] * Second[0] + First[2] * Second[1];
		_mat[2] = 1.0 - _mat[0];
		_mat[1] = First[1] * Second[0] + First[3] * Second[1];
		_mat[3] = 1.0 - _mat[1];
	}

	double operator()(size_t i, size_t j) const {
		// get element from matrix
		// 0  2
		// 1  3
		return _mat[i + 2 * j];
	}

	double operator[](size_t i) const { return _mat[i]; }

	void print() const {
		uint cols = _mat.n_cols;
		uint rows = _mat.n_rows;

		std::cout << "---------------" << std::endl;
		for (uint i = 0; i < rows; i++) {
			for (uint j = 0; j < cols; j++) { std::cout << _mat(i, j) << " "; }
			std::cout << std::endl;
		}
		std::cout << "---------------" << std::endl;
	}
};

/** @brief Class to store the matrices for each bin.
 */
class TMatrices {
private:
	std::vector<TMatrix> _matrices;
	TMatrix _matrix_alpha;
	arma::mat _lambda_c = arma::zeros(2, 2);

public:
	TMatrices() = default;
	explicit TMatrices(size_t NumBins) { resize(NumBins); }

	/** @brief Resize the vector of matrices to the given number of bins.
	 * @param NumBins
	 */
	void resize(size_t NumBins) { _matrices.resize(NumBins); }

	/** @brief Get the number of matrices.
	 * @return Number of matrices.
	 */
	[[nodiscard]] size_t size() const { return _matrices.size(); }

	/** @brief Get the vector of matrices.
	 * @return Vector of matrices.
	 */
	const std::vector<TMatrix> &get_matrices() const { return _matrices; }
	const TMatrix &operator[](size_t i) const { return _matrices[i]; }

	/** @brief Set the matrix lambda for the clique given the two rate parameters.
	 * @param mu_c_1
	 * @param mu_c_0
	 */
	void set_lambda(double mu_c_1, double mu_c_0) {
		_lambda_c[0] = -mu_c_1;
		_lambda_c[1] = mu_c_0;
		_lambda_c[2] = mu_c_1;
		_lambda_c[3] = -mu_c_0;
	}

	/** @brief Set the matrices for each bin. Instead of calculating the matrix exponential for each bin, the matrix
	 * exponential of the scaling matrix is calculated once and then multiplied with the previous matrix. This is
	 * mathematically equevalent to calculating the matrix exponential for each bin.
	 * @param a
	 * @param Delta
	 */
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
	static void print_mat(const TMatrix &my_matrix) {
		uint cols = my_matrix.get_matrix().n_cols;
		uint rows = my_matrix.get_matrix().n_rows;

		std::cout << "---------------" << std::endl;
		for (uint i = 0; i < rows; i++) {
			for (uint j = 0; j < cols; j++) { std::cout << my_matrix(i, j) << " "; }
			std::cout << std::endl;
		}
		std::cout << "---------------" << std::endl;
	}
};

/** Class representing a clique in our model. A clique is defined as having a set of nodes that are all leaves in
 * all dimensions except one. Each clique has a set of matrices, the change rate parameters, and the start index of the
 * nodes in the tree. The start index, the variable dimension, and the number of nodes are needed to get the correct
 * indices in our multidimensional space Y and Z.
 */
class TClique {
private:
	// transition matrix and parameters
	TMatrices _matrices;
	double _mu_c_1;
	double _mu_c_0;

	// info about size and dimensionality of clique
	std::vector<size_t> _start_index_in_leaves_space;
	size_t _variable_dimension;
	size_t _n_nodes;
	size_t _increment;

	// count the number of leaves with value 1 in a clique
	TypeCounter1 _counter_leaves_state_1 = 0;

	void _update_current_state(TStorageZVector &Z, const TCurrentState &current_state, size_t index_in_tree,
	                           bool new_state, std::vector<TStorageZ> &linear_indices_in_Z_space_to_insert,
	                           const TTree &tree) const;

	static void _calculate_log_prob_root(double stationary_0, std::array<coretools::TSumLogProbability, 2> &sum_log);
	void _calculate_log_prob_node_to_children(size_t index_in_tree, const TTree &tree,
	                                          const TCurrentState &current_state,
	                                          std::array<coretools::TSumLogProbability, 2> &sum_log) const;

	template<typename ContainerStates> // can either be TSheet or TCurrentStates
	bool _getState(const ContainerStates &states, size_t parent_index_in_tree,
	               size_t leaf_index_in_tree_of_last_dim) const {
		if constexpr (std::is_same_v<ContainerStates, TSheet>) { // is a sheet
			return states.get(parent_index_in_tree, leaf_index_in_tree_of_last_dim);
		} else { // TCurrentState
			return states.get(parent_index_in_tree);
		}
	}

public:
	TClique(const std::vector<size_t> &start_index, size_t variable_dimension, size_t n_nodes, size_t increment);
	~TClique() = default;

	/// @brief Initialize the matrices for the clique.
	/// @param a The lower bound of the bin.
	/// @param delta The bin width.
	/// @param n_bins The number of bins.
	void initialize(double a, double delta, size_t n_bins) {
		_matrices.resize(n_bins);
		_matrices.set(a, delta);
	}

	void set_mus(double mu_c_1, double mu_c_0) {
		_mu_c_1 = mu_c_1;
		_mu_c_0 = mu_c_0;
	}

	double get_stationary_probability(bool state) const {
		if (state) { return _mu_c_1 / (_mu_c_1 + _mu_c_0); }
		return _mu_c_0 / (_mu_c_1 + _mu_c_0);
	}

	/// @brief Set the rate parameters for the clique.
	void set_lambda() { _matrices.set_lambda(_mu_c_1, _mu_c_0); }

	/// @brief Returns the matrices for the clique.
	/// @return The class containing the matrices.
	[[nodiscard]] const TMatrices &get_matrices() const { return _matrices; }

	/// @brief Update the Z dimension for this clique.
	/// @param Y The current state of the Y dimension.
	/// @param Z The current state of the Z dimension.
	/// @param tree The tree.
	std::vector<TStorageZ> update_Z(const TStorageYVector &Y, TStorageZVector &Z, const TTree &tree) const;

	size_t get_number_of_nodes() const { return _n_nodes; }

	template<typename ContainerStates> // can either be TSheet or TCurrentStates
	void calculate_log_prob_parent_to_node(size_t index_in_tree, const TTree &tree,
	                                       size_t leaf_index_in_tree_of_last_dim, const ContainerStates &states,
	                                       std::array<coretools::TSumLogProbability, 2> &sum_log) const {
		// calculates log P(node = 0 | parent) and log P(node = 1 | parent)
		const auto &node                  = tree.get_node(index_in_tree);
		auto bin_length                   = node.get_branch_length_bin();
		const auto &matrix_for_bin        = _matrices[bin_length];
		const size_t parent_index_in_tree = node.parentIndex_in_tree();
		for (size_t i = 0; i < 2; ++i) { // loop over possible values (0 or 1) of the node
			const bool state_of_parent = _getState(states, parent_index_in_tree, leaf_index_in_tree_of_last_dim);
			sum_log[i].add(matrix_for_bin(state_of_parent, i));
		}
	}

	void update_counter_leaves_state_1(bool new_state, bool old_state);
	void update_counter_leaves_state_1(int difference);
};

bool sample(std::array<coretools::TSumLogProbability, 2> &sum_log);

#endif // ACOL_TBRANCHLENGTHS_H
