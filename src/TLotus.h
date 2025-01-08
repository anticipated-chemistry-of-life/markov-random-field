#pragma once

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
	using TypeParamDelta = stattools::TParameter<SpecDelta, BoxType>;

private:
	// trees should be a const ref because we don't want to change the trees and don't want to copy them
	const std::vector<TTree> &_trees;
	TStorageYVector _L_sm;
	TStorageYVector _x_sm;
	std::vector<size_t> _species_counter;
	std::vector<size_t> _molecules_counter;

	// parameters gamma and delta
	TypeParamGamma *_gamma = nullptr;
	TypeParamDelta *_delta = nullptr;

	// temporary values
	double _oldLL;
	double _curLL;

	// private functions
	static std::vector<size_t> _get_dimensions_Lotus_space(const std::vector<TTree> &trees) {
		return {trees[0].get_number_of_leaves(), trees[1].get_number_of_leaves()};
	}

	/// To construct the _data_X_of_Lotus we will need to have as input a vector of TStorageYVector but this time with
	/// as many dimensions as there are trees.
	void _initialize_x_sm(const TStorageYVector &Y);

	double _calculate_probability_of_L_sm(bool x_sm, bool L_sm, size_t linear_index_in_L_space) const;

	void _simulateUnderPrior(Storage *) override;

public:
	TLotus(const std::vector<TTree> &trees, TypeParamGamma *gamma, TypeParamDelta *delta);
	~TLotus() = default;

	[[nodiscard]] std::string name() const override;
	void initialize() override;

	void load_from_file(const std::string &filename);

	double calculate_research_effort(size_t linear_index_in_L_space) const;
	double calculate_log_likelihood_of_L() const;

	double getSumLogPriorDensity(const Storage &) const override;

	[[nodiscard]] double calculateLLRatio(TypeParamGamma *, size_t /*Index*/, const Storage &);
	[[nodiscard]] double calculateLLRatio(TypeParamDelta *, size_t /*Index*/, const Storage &);

	void updateTempVals(TypeParamGamma *, size_t /*Index*/, bool Accepted);
	void updateTempVals(TypeParamDelta *, size_t /*Index*/, bool Accepted);

	void guessInitialValues() override;

	void set_x(bool state, size_t linear_index_in_L_space);
};
