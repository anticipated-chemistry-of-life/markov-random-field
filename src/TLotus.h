#pragma once

#include "TCollapser.h"
#include "TMarkovField.h"
#include "TStorageYVector.h"
#include "TTree.h"
#include "Types.h"
#include "coretools/Files/TOutputFile.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Types/probability.h"
#include "stattools/ParametersObservations/TParameter.h"
#include "stattools/Priors/TPriorBase.h"
#include <cstddef>
#include <deque>
#include <string>
#include <utility>
#include <vector>

class TLotus : public stattools::prior::TBaseLikelihoodPrior<TypeLotus, NumDimLotus> {
public:
	// some type aliases, for better readability
	using BoxType = TLotus;
	using Base    = stattools::prior::TBaseLikelihoodPrior<TypeLotus, NumDimLotus>;
	using typename Base::Storage;
	using typename Base::UpdatedStorage;

	using TypeParamGamma = stattools::TParameter<SpecGamma, BoxType>;

private:
	// trees should be a const ref because we don't want to change the trees and don't want to copy them
	const std::vector<std::unique_ptr<TTree>> &_trees;

	// data
	TStorageYVector _L;
	std::vector<std::vector<size_t>> _occurrence_counters;

	// Markov field
	TMarkovField _markov_field;

	// Markov field parameter (only needed for stattools purposes to build a valid DAG)
	const std::vector<std::unique_ptr<stattools::TParameter<SpecMarkovField, TLotus>>> &_markov_field_stattools_param;

	// how to collapse
	TCollapser _collapser;

	// parameters gamma
	TypeParamGamma *_gamma = nullptr;
	double _epsilon        = 0.0001;

	// temporary values
	double _oldLL;
	double _curLL;
	TCurrentState _tmp_state_along_last_dim;

	// output file
	std::string _prefix;

	// private functions
	double _calculate_research_effort(const std::vector<size_t> &index_in_collapsed_space) const;
	double _calculate_probability_of_L_given_x(bool x, bool L,
	                                           const std::vector<size_t> &index_in_collapsed_space) const;
	double _calculate_probability_of_L_given_x(bool x, bool L, size_t linear_index_in_collapsed_space) const;
	double _calculate_log_likelihood_of_L_no_collapsing() const;
	double _calculate_log_likelihood_of_L_do_collapse() const;
	void _simulateUnderPrior(Storage *) override;

public:
	TLotus(std::vector<std::unique_ptr<TTree>> &trees, TypeParamGamma *gamma, size_t n_iterations,
	       const std::vector<std::unique_ptr<stattools::TParameter<SpecMarkovField, TLotus>>>
	           &markov_field_stattools_param,
	       std::string prefix);
	~TLotus() override = default;

	[[nodiscard]] std::string name() const override;
	void initialize() override;
	void load_from_file(const std::string &filename);
	void guessInitialValues() override;

	void burninHasFinished() override;
	void MCMCHasFinished() override;

	double calculate_log_likelihood_of_L() const;
	double getSumLogPriorDensity(const Storage &) const override;

	void fill_tmp_state_along_last_dim(const std::vector<size_t> &start_index_clique_along_last_dim, size_t K);
	void calculate_LL_update_Y(const std::vector<size_t> &index_in_leaves_space, bool old_state,
	                           std::array<double, 2> &prob) const;
	void update_cur_LL(double cur_LL);

	void update_markov_field(size_t iteration);
	[[nodiscard]] double calculateLLRatio(TypeParamGamma *, size_t /*Index*/, const Storage &);
	void updateTempVals(TypeParamGamma *, size_t /*Index*/, bool Accepted);

	const TStorageYVector &get_Lotus() const;

	static std::string get_filename_lotus() { return coretools::instances::parameters().get("lotus"); }
};
