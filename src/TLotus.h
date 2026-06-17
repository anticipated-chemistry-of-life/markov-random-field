#pragma once

#include "TCollapser.h"
#include "TCurrentState.h"
#include "TMarkovField.h"
#include "Types.h"
#include "cli.h"
#include "constants.h"
#include "ntfy/TNtfyNotifier.h"
#include "stattools/ParametersObservations/TParameter.h"
#include "stattools/Priors/TPriorBase.h"
#include "storages/y_storage/TStorageYMatrix.h"
#include "tree/TTree.h"
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

class TLotus : public stattools::prior::TBaseLikelihoodPrior<stattools::TObservationBase, TypeLotus,
                                                             NumDimLotus> {
public:
	// some type aliases, for better readability
	using BoxType = TLotus;
	using Base =
	    stattools::prior::TBaseLikelihoodPrior<stattools::TObservationBase, TypeLotus, NumDimLotus>;
	using typename Base::Storage;
	using typename Base::UpdatedStorage;

	using TypeParamGamma     = stattools::TParameter<SpecGamma, BoxType>;
	using TypeParamErrorRate = stattools::TParameter<SpecErrorRate, BoxType>;

private:
	// trees should be a const ref because we don't want to change the trees and don't want to copy
	// them
	const std::vector<std::unique_ptr<TTree>> &_trees;

	// data
	TStorageYMatrix _L;

	/// The occurences is the log(count + 1) of the paper counts
	std::vector<std::vector<double>> _occurrence_counters;

	// Markov field
	TMarkovField _markov_field;

	// Markov field parameter (only needed for stattools purposes to build a valid DAG)
	const std::vector<std::unique_ptr<stattools::TParameter<SpecMarkovField, TLotus>>>
	    &_markov_field_stattools_param;

	// how to collapse
	TCollapser _collapser;

	// parameters gamma
	TypeParamGamma *_gamma = nullptr;

	// Error rate of lotus
	TypeParamErrorRate *_error_rate = nullptr;

	// Error rate of simple error model
	double _epsilon = 0.0001;

	// temporary values
	double _oldLL;
	double _curLL;
	TCurrentState _tmp_state_along_last_dim;

	// output file
	std::string _prefix;

	// simulate or infer?
	bool _simulate = false;

	// notifications
	TNtfyNotifier _notifier;

	// counts how many burnin rounds have finished (oneBurninHasFinished no longer
	// receives the round number from stattools)
	size_t _burnin_round = 0;

	// monotonic counter of markov field updates (the update function registered via
	// addFuncToUpdate is no longer passed the iteration number by stattools)
	size_t _mrf_update_iteration = 0;

	// private functions
	[[nodiscard]] double
	_calculate_research_effort(const std::vector<size_t> &index_in_collapsed_space) const;
	[[nodiscard]] double
	_calculate_probability_of_L_given_x(bool x, bool L,
	                                    const std::vector<size_t> &index_in_collapsed_space) const;
	[[nodiscard]] double
	_calculate_probability_of_L_given_x(bool x, bool L,
	                                    size_t linear_index_in_collapsed_space) const;
	[[nodiscard]] double _calculate_log_likelihood_of_L_no_collapsing() const;
	[[nodiscard]] double _calculate_log_likelihood_of_L_do_collapse() const;
	void _simulateUnderPrior(Storage *) override;
	[[nodiscard]] std::pair<bool, size_t>
	_get_state_of_Y(size_t i, size_t index_in_Y,
	                const std::vector<std::pair<size_t, TStorageY>> &y_entries) const;
	[[nodiscard]] std::pair<bool, size_t>
	_get_state_of_L(size_t i, size_t index_in_L,
	                const std::vector<std::pair<size_t, TStorageY>> &l_entries) const;

	double _return_error_rate(bool L) const;

public:
	TLotus(std::vector<std::unique_ptr<TTree>> &trees, TypeParamGamma *gamma,
	       TypeParamErrorRate *error_rate, size_t n_iterations,
	       const std::vector<std::unique_ptr<stattools::TParameter<SpecMarkovField, TLotus>>>
	           &markov_field_stattools_param,
	       std::string prefix, bool simulate);
	~TLotus() override = default;

	[[nodiscard]] std::string name() const override;
	void initialize() override;
	void load_from_file(const std::string &filename);
	void guessInitialValues() override;

	void burninHasFinished() override;
	void oneBurninHasFinished() override;
	void MCMCHasFinished() override;

	[[nodiscard]] double calculate_log_likelihood_of_L() const;
	[[nodiscard]] double getSumLogPriorDensity(const Storage &) const override;

	void fill_tmp_state_along_last_dim(const std::vector<size_t> &start_index_clique_along_last_dim,
	                                   size_t K);
	void calculate_LL_update_Y(const std::vector<size_t> &index_in_leaves_space,
	                           size_t index_for_tmp_state, bool old_state,
	                           std::array<double, 2> &prob) const;
	void update_cur_LL(double cur_LL);

	void update_markov_field();
	[[nodiscard]] double calculateLLRatio(TypeParamGamma *, size_t /*Index*/);
	double calculateLLRatio(TypeParamErrorRate *, size_t /*Index*/);
	void updateTempVals(TypeParamGamma *, size_t /*Index*/, bool Accepted);
	void updateTempVals(TypeParamErrorRate *, size_t /*Index*/, bool Accepted);

	[[nodiscard]] const TStorageYMatrix &get_Lotus() const;

	static std::string get_filename_lotus() { return ProgramOptions::LOTUS_FILENAME; }
	[[nodiscard]] const TCurrentState &get_tmp_state_along_last_dim() const {
		return _tmp_state_along_last_dim;
	}

	[[nodiscard]] const TMarkovField &get_markov_field() const { return _markov_field; }
};
