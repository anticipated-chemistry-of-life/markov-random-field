#pragma once

#include "TCollapser.h"
#include "TStorageYVector.h"
#include "TTree.h"
#include "Types.h"
#include "stattools/ParametersObservations/TParameter.h"
#include "stattools/Priors/TPriorBase.h"
#include <cstddef>
#include <string>
#include <vector>

// For this class, we need to enforce that the first tree will always be the one of the species and
// the second tree will always be the one of the molecules. The rest of the trees, we won't care.
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
	const std::vector<TTree> &_trees;

	// data
	TStorageYVector _L;
	std::vector<std::vector<size_t>> _occurrence_counters;

	// how to collapse
	TCollapser _collapser;

	// parameters gamma
	TypeParamGamma *_gamma = nullptr;

	// temporary values
	double _oldLL;
	double _curLL;
	TCurrentState _tmp_state_along_last_dim;

	// private functions
	double _calculate_research_effort(const std::vector<size_t> &index_in_collapsed_space) const;
	double _calculate_probability_of_L_given_x(bool x, bool L,
	                                           const std::vector<size_t> &index_in_collapsed_space) const;

	void _simulateUnderPrior(Storage *) override;

public:
	TLotus(const std::vector<TTree> &trees, TypeParamGamma *gamma);
	TLotus(const std::vector<TTree> &trees);
	~TLotus() override = default;

	[[nodiscard]] std::string name() const override;
	void initialize() override;

	void load_from_file(const std::string &filename);

	double calculate_log_likelihood_of_L() const;

	double getSumLogPriorDensity(const Storage &) const override;

	void fill_tmp_state_along_last_dim(const std::vector<size_t> &start_index_in_leaves);
	void calculate_LL_update_Y(const std::vector<size_t> &index_in_leaves_space, bool new_state, bool old_state,
	                           std::array<coretools::TSumLogProbability, 2> &sum_log);

	[[nodiscard]] double calculateLLRatio(TypeParamGamma *, size_t Index, const Storage &);
	void updateTempVals(TypeParamGamma *, size_t Index, bool Accepted);

	void guessInitialValues() override;

	const TStorageYVector &get_Lotus() const { return _L; }
};
