//
// Created by madleina on 22.10.24.
//

#ifndef ACOL_TBRANCHLENGTHS_H
#define ACOL_TBRANCHLENGTHS_H

#include "TCurrentState.h"
#include "TStorageYVector.h"
#include "TStorageZ.h"
#include "TStorageZVector.h"
#include "Types.h"
#include "coretools/Main/TError.h"
#include "coretools/Math/TSumLog.h"
#include "coretools/devtools.h"
#include <armadillo>
#include <cstddef>
#include <iomanip>
#include <stdexcept>
#include <unistd.h>
#include <vector>

class TTree;

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
	void set_from_matrix_exponential(const arma::mat &Lambda) {
		try {
			_mat = arma::expmat(Lambda);
		} catch (const std::runtime_error &e) {
			DEVERROR("can't perform matrix exponential from the following matrix : ");
			Lambda.print();
		}
	}

	/** @brief Get the matrix.
	 * @return Matrix.
	 */
	[[nodiscard]] arma::mat get_matrix() const { return _mat; }

	/** @brief Perform matrix multiplication of two matrices. Since the matrix rows have to sum to 1, the second
	 * value of the row can easily be calculted.
	 * @param First
	 * @param Second
	 */
	void set_from_product(const TMatrix &First, const TMatrix &Second) {
		// do matrix multiplication
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
		auto cols = _mat.n_cols;
		auto rows = _mat.n_rows;

		std::cout << "---------------" << std::endl;
		for (auto i = 0; i < rows; i++) {
			for (auto j = 0; j < cols; j++) { std::cout << std::setprecision(20) << _mat(i, j) << " "; }
			std::cout << std::endl;
		}
		std::cout << "---------------" << std::endl;
		OUT("Row sum: ", _mat(0, 0) + _mat(0, 1), " ", _mat(1, 0) + _mat(1, 1));
	}
};

/** @brief Class to store the matrices for each bin.
 */
class TMatrices {
private:
	std::vector<TMatrix> _matrices;
	arma::mat _lambda_c = arma::zeros(2, 2);

	static double _Delta;

	/** @brief Set the matrices for each bin. Instead of calculating the matrix exponential for each bin, the matrix
	 * exponential of the scaling matrix is calculated once and then multiplied with the previous matrix. This is
	 * mathematically equevalent to calculating the matrix exponential for each bin.
	 */
	void _fill_matrices() {
		// calculate matrix exponential for first bin
		TMatrix P_0;
		P_0.set_from_matrix_exponential(_lambda_c * _Delta / 2.0);

		// calculate matrix exponential of scaling matrix
		TMatrix matrix_alpha;
		matrix_alpha.set_from_matrix_exponential(_lambda_c * _Delta);

		_matrices[0] = P_0;
		// do recursion
		for (size_t k = 1; k < _matrices.size(); ++k) { _matrices[k].set_from_product(_matrices[k - 1], matrix_alpha); }
	}

public:
	TMatrices() = default;
	explicit TMatrices(size_t NumBins, double Delta) { resize(NumBins, Delta); }

	/** @brief Resize the vector of matrices to the given number of bins.
	 * @param NumBins
	 */
	void resize(size_t NumBins, double Delta) {
		_matrices.resize(NumBins);
		_Delta = Delta;
	}

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
	void set_lambda(double alpha, TypeNu nu) {
		_lambda_c[0] = (-alpha) * nu;
		_lambda_c[1] = (1 - alpha) * nu;
		_lambda_c[2] = alpha * nu;
		_lambda_c[3] = (alpha - 1) * nu;

		_fill_matrices();
	}

	static void print_mat(const TMatrix &my_matrix) {
		auto cols = my_matrix.get_matrix().n_cols;
		auto rows = my_matrix.get_matrix().n_rows;

		std::cout << "---------------" << std::endl;
		for (auto i = 0; i < rows; i++) {
			for (auto j = 0; j < cols; j++) { std::cout << my_matrix(i, j) << " "; }
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
	using TypeParamBinBranches = stattools::TParameter<SpecBinnedBranches, TTree>;

	// transition matrix and parameters
	TMatrices _cur_matrices;
	TMatrices _try_matrices;

	// info about size and dimensionality of clique
	std::vector<size_t> _start_index_in_leaves_space;
	size_t _variable_dimension;
	size_t _n_nodes;
	size_t _increment;

	// count the number of leaves with value 1 in a clique
	TypeCounter1 _counter_leaves_state_1 = 0;

	/// @brief Calculates the log probability of the root given the stationary distribution
	static void _calculate_log_prob_root(double stationary_0, std::array<coretools::TSumLogProbability, 2> &sum_log);

	template<typename ContainerStates> // can either be TSheet or TCurrentStates
	inline bool _getState(const ContainerStates &states, size_t parent_index_in_tree,
	                      size_t leaf_index_in_tree_of_last_dim) const {
		if constexpr (std::is_same_v<ContainerStates, TSheet>) { // is a sheet
			return states.get(parent_index_in_tree, leaf_index_in_tree_of_last_dim);
		} else { // TCurrentState
			return states.get(parent_index_in_tree);
		}
	}

	void _update_current_state(TStorageZVector &Z, TCurrentState &current_state, size_t index_in_tree, bool new_state,
	                           std::vector<TStorageZ> &linear_indices_in_Z_space_to_insert, const TTree *tree) const;

	/// @brief Calculates the log probability of a node to its children
	void _calculate_log_prob_node_to_children(
	    size_t index_in_tree, const TTree *tree, const TCurrentState &current_state,
	    std::array<coretools::TSumLogProbability, 2> &sum_log, const TypeParamBinBranches *binned_branch_lengths,
	    const std::vector<size_t> &leaves_and_internal_nodes_without_roots_indices) const;

	/// @brief Sets Z given the maximal likelihood given its children. This was created to avoid that Z is stuck in a
	/// state and cannot change.
	/// @param node_index The index of the internal node we want to set
	/// @param current_state The current state of the clique.
	/// @param Z the Z vector of that tree (i.e that clique)
	/// @param tree the tree of interest
	/// @param binned_branch_lengths the vector of branch length
	/// @param leaves_and_internal_nodes_without_roots_indices Same as the variable name
	/// @param linear_indices_in_Z_space_to_insert Same as the variable name
	void _set_Z_to_MLE(size_t node_index, TCurrentState &current_state, TStorageZVector &Z, const TTree *tree,
	                   const TypeParamBinBranches *binned_branch_lengths,
	                   const std::vector<size_t> &leaves_and_internal_nodes_without_roots_indices,
	                   std::vector<TStorageZ> &linear_indices_in_Z_space_to_insert) const;

	static size_t _get_parent_index(size_t index_in_tree, const TTree *tree);

public:
	TClique(const std::vector<size_t> &start_index, size_t variable_dimension, size_t n_nodes, size_t increment);
	~TClique() = default;

	/// @brief Initialize the matrices for the clique.
	/// @param delta The bin width.
	/// @param n_bins The number of bins.
	void initialize(double delta, size_t n_bins) {
		_cur_matrices.resize(n_bins, delta);
		_try_matrices.resize(n_bins, delta);
	}

	/// Gets the stationary probability for state 0 or 1.
	static double get_stationary_probability(bool state, double alpha) {
		if (state) { return alpha; }
		return 1.0 - alpha;
	}

	/// @brief Set the rate parameters for the clique.
	void set_lambda(TypeAlpha alpha, TypeNu nu) {
		_cur_matrices.set_lambda(alpha, nu);
		_try_matrices = _cur_matrices;
	}

	/// @brief Sets the "try Matrix" to the given values
	void update_lambda(double alpha, double nu) { _try_matrices.set_lambda(alpha, nu); }

	/// @brief If the "try matrix" is accepted, then we change our "current matrix" to the "try matrix"
	void accept_update_mu() { _cur_matrices = _try_matrices; }

	/// @brief Returns the matrices for the clique.
	/// @return The class containing the matrices.
	[[nodiscard]] const TMatrices &get_matrices() const { return _cur_matrices; }

	/// @brief Update the Z dimension for this clique.
	/// @param Y The current state of the Y dimension.
	/// @param Z The current state of the Z dimension.
	/// @param tree The tree.
	std::vector<TStorageZ> update_Z(std::vector<double> &joint_log_prob_density, TCurrentState &current_state,
	                                TStorageZVector &Z, const TTree *tree, TypeAlpha alpha,
	                                const TypeParamBinBranches *binned_branch_lengths,
	                                const std::vector<size_t> &leaves_and_internal_nodes_without_roots_indices) const;

	std::vector<TStorageZ>
	initialize_Z_from_children(TCurrentState &current_state, TStorageZVector &Z, const TTree *tree,
	                           const TypeParamBinBranches *binned_branch_lengths,
	                           const std::vector<size_t> &leaves_and_internal_nodes_without_roots_indices) const;

	TCurrentState create_current_state(const TStorageYVector &Y, TStorageZVector &Z, const TTree &tree);

	/// @brief Return the number of nodes in the clique
	/// @return Return the number of nodes in the clique
	size_t get_number_of_nodes() const { return _n_nodes; }

	/// @brief Gets the matrix at the corresponding bin
	/// @param bin_length: a size_t the is the index where to get the matrix
	/// @return A reference to the asked matrix
	template<bool UseTry> const TMatrix &get_matrix(size_t bin_length) const {
		if constexpr (UseTry) { return _try_matrices[bin_length]; }
		return _cur_matrices[bin_length];
	}

	/// @brief Calculates the log probability of a node to its parent.
	template<typename ContainerStates, bool UseTry = false> // ContainerStates can either be TSheet or TCurrentStates
	void calculate_log_prob_parent_to_node(size_t index_in_tree, TypeBinnedBranchLengths binned_branch_length,
	                                       const TTree *tree, size_t leaf_index_in_tree_of_last_dim,
	                                       const ContainerStates &states,
	                                       std::array<coretools::TSumLogProbability, 2> &sum_log) const {
		const size_t parent_index_in_tree = _get_parent_index(index_in_tree, tree);
		const auto &matrix_for_bin        = get_matrix<UseTry>(binned_branch_length);
		for (size_t i = 0; i < 2; ++i) { // loop over possible values (0 or 1) of the node
			const bool state_of_parent = _getState(states, parent_index_in_tree, leaf_index_in_tree_of_last_dim);
			sum_log[i].add(matrix_for_bin(state_of_parent, i));
		}
	}

	/// Updates the counter of the clique. This is used by the collapser.
	void update_counter_leaves_state_1(bool new_state, bool old_state);
	void update_counter_leaves_state_1(int difference);
	size_t get_counter_leaves_state_1() const { return _counter_leaves_state_1; }

	/// @return Returns the jump size of the clique
	size_t get_increment() const { return _increment; }

	/// @return Returns the start index in the leaf space of that specific clique
	const std::vector<size_t> &get_start_index_in_leaf_space() const { return _start_index_in_leaves_space; }

	template<bool UseTryMatrix>
	double calculate_prob_to_parent(size_t index_in_tree, const TTree *tree,
	                                TypeBinnedBranchLengths binned_branch_length,
	                                const TCurrentState &current_state) const {
		// always use cur matrix
		const auto &matrix = get_matrix<UseTryMatrix>(binned_branch_length);

		size_t parent_index = _get_parent_index(index_in_tree, tree);

		bool parent_state  = current_state.get(parent_index);
		bool child_state   = current_state.get(index_in_tree);
		double parent_prob = matrix(parent_state, child_state); // from parent_state to child_state
		return parent_prob;
	}
};

bool sample(std::array<coretools::TSumLogProbability, 2> &sum_log);
bool sample(double log_prob_0, double log_prob_1);
inline size_t getLinearIndexSkippingDim(const std::vector<size_t> &index, size_t skip_dim,
                                        const std::vector<size_t> &dims) {
	// check if size matches and if coordinates are within dimensions
	assert(index.size() == dims.size());
	assert(([dims = std::as_const(dims), &index = std::as_const(index)]() constexpr {
		for (size_t i = 0; i < dims.size(); i++) {
			if (index[i] >= dims[i]) { return false; }
		}
		return true;
	})());

	size_t linear_index;
	if (skip_dim == index.size() - 1) {
		linear_index = 0;
	} else {
		linear_index = index.back();
	}
	size_t prod = 1;

	for (size_t i = dims.size() - 1; i > 0; i--) {
		prod *= dims[i];
		if (i - 1 == skip_dim) { continue; };
		linear_index += prod * index[i - 1];
	}
	return linear_index;
}

#endif // ACOL_TBRANCHLENGTHS_H
