//
// The "simple error model" data source.
//
// D is a binary matrix with exactly the same dimensions as the latent field Y: each cell of D is a
// direct observation of the corresponding cell of Y, correct with probability 1 - epsilon and
// inverted with probability epsilon. Unlike LOTUS, D carries no notion of
// research effort -- it is deliberately the simplest possible data source, so that switching every
// other source off isolates whether a convergence problem lives in a data likelihood or in the
// Markov field itself.
//
// This class is not a stattools box. Its parameter (epsilon_simple_model) hangs off TDataModel,
// which owns this object and forwards the MCMC callbacks; see TDataModel.h.
//

#pragma once

#include "simple_error_model/TSimpleErrorModelMath.h"

#ifdef USE_SIMPLE_ERROR_MODEL

#include "Types.h"
#include "constants.h"
#include "stattools/ParametersObservations/TParameter.h"
#include "storages/storage_backend.h"
#include "tree/TTree.h"
#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

class TSimpleErrorModel {
public:
	using TypeParamEpsilon = stattools::TParameter<SpecEpsilonSimpleModel, TDataModel>;

private:
	/// trees are a const ref: we neither modify nor copy them
	const std::vector<std::unique_ptr<TTree>> &_trees;

	/// the observed data, same dimensions as Y
	TFieldStorage _D;

	/// error rate; owned by TModel, updated by stattools
	TypeParamEpsilon *_epsilon = nullptr;

	/// Number of cells where D and Y currently disagree, and the total number of cells. The
	/// likelihood depends on the data only through these two numbers (see TSimpleErrorModelMath.h),
	/// so keeping the count up to date after each field update makes every epsilon move O(1).
	size_t _n_disagree  = 0;
	size_t _total_cells = 0;

	[[nodiscard]] double _eps() const { return (double)_epsilon->value(); }

public:
	TSimpleErrorModel(const std::vector<std::unique_ptr<TTree>> &trees, TypeParamEpsilon *epsilon);
	~TSimpleErrorModel() = default;

	/// Sizes D from the trees. Called in both inference and simulation: during simulation D starts
	/// empty and is filled by simulate_D_from_Y.
	void initialize_storage();

	void load_from_file(const std::string &filename);

	/// Synchronises the disagreement count with the current Y. Must be called once before the first
	/// likelihood evaluation, and after anything that changes Y outside of an update.
	void guess_initial_values(const TFieldStorage &Y);

	// --- hooks used by the field update (see TMarkovField::_update_Y) ---

	/// D's cells for one row of the field, as a window: `n_cells` cells from `start_index`, one
	/// after the other. D has the field's dimensions, so the field's index is already D's. The
	/// update opens one window per species leaf and reads it. Nothing writes D.
	[[nodiscard]] TFieldStorage::TWindow open_row(const IndexArray &start_index, size_t n_cells) {
		return _D.open_window(start_index, n_cells, /*stride=*/1);
	}

	/// prob[0] = P(D_cell | Y = 0), prob[1] = P(D_cell | Y = 1). `observed_state` is the state D
	/// holds for the cell, which the caller reads from the window above.
	void probabilities_for_Y_update(bool observed_state, std::array<double, 2> &prob) const {
		simple_error_model::probabilities_for_both_Y_states(observed_state, _eps(), prob);
	}

	/// Installs the disagreement count accumulated over a full field update.
	void set_n_disagree(size_t n_disagree) { _n_disagree = n_disagree; }

	// --- likelihood ---

	[[nodiscard]] double log_likelihood() const {
		return simple_error_model::log_likelihood_from_counts(_total_cells, _n_disagree, _eps());
	}

	/// Log-likelihood ratio between two epsilon values. O(1): the disagreement count does not
	/// depend on epsilon, so nothing has to be re-scanned.
	[[nodiscard]] double log_likelihood_ratio(double old_eps, double new_eps) const {
		return simple_error_model::log_likelihood_ratio(_total_cells, _n_disagree, old_eps,
		                                                new_eps);
	}

	// --- simulation ---

	/// Draws every cell of D from the corresponding cell of Y at the current error rate.
	void simulate_D_from_Y(const TFieldStorage &Y);

	/// Writes the simulated D as <prefix>_simulated_simple_data.tsv: a header naming every tree,
	/// then one row of leaf node ids per cell whose state is 1 (same sparse format as the LOTUS
	/// file, so it can be read straight back by load_from_file).
	void write_simulated_D(const std::string &prefix) const;

	// --- accessors ---

	[[nodiscard]] const TFieldStorage &get_D() const { return _D; }
	[[nodiscard]] size_t n_disagree() const { return _n_disagree; }
};

#endif // USE_SIMPLE_ERROR_MODEL
