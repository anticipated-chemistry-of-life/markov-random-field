#include "constants.h"
#include "mass_spec/feature_likelihood.h"
#include "mass_spec/msms_run.h"
#include "gtest/gtest.h"
#include <stdexcept>
#include <vector>

TEST(TFeatureLikelihood_Tests, set_linear_index_1) {
	TFeatureLikelihood storage(16777214, 125);
	EXPECT_EQ(storage.get_molecule_index(), 16777214);
	EXPECT_EQ(storage.get_binned_likelihood(), 125);
}

TEST(TFeatureLikelihood_Tests, set_linear_index_2) {
	TFeatureLikelihood storage(0, 200);
	EXPECT_EQ(storage.get_molecule_index(), 0);
	EXPECT_EQ(storage.get_binned_likelihood(), 200);
}

TEST(TFeatureLikelihood_Tests, set_linear_index_3) {
	TFeatureLikelihood storage(1, 200);
	EXPECT_EQ(storage.get_molecule_index(), 1);
	EXPECT_EQ(storage.get_binned_likelihood(), 200);
}

TEST(TFeatureLikelihood_Tests, set_unknown_molecule) {
	auto storage = TFeatureLikelihood::new_unknown_molecule(128);
	EXPECT_TRUE(storage.is_unknown_molecule());
	EXPECT_EQ(storage.get_binned_likelihood(), 128);
	EXPECT_EQ(storage.get_molecule_index(), TFeatureLikelihood::get_unknown_molecule_index());
	EXPECT_EQ(storage.get_molecule_index(), MAX_NUMBER_OF_MOLECULES);

	auto storage2 = TFeatureLikelihood::new_unknown_molecule(0);
	EXPECT_TRUE(storage2.is_unknown_molecule());
	EXPECT_EQ(storage2.get_binned_likelihood(), 0);
	EXPECT_EQ(storage2.get_molecule_index(), TFeatureLikelihood::get_unknown_molecule_index());
	EXPECT_EQ(storage2.get_molecule_index(), MAX_NUMBER_OF_MOLECULES);

	auto storage3 = TFeatureLikelihood(0, 255);
	EXPECT_FALSE(storage3.is_unknown_molecule());
	EXPECT_EQ(storage3.get_binned_likelihood(), 255);
	EXPECT_EQ(storage3.get_molecule_index(), 0);
}

TEST(TFeatureLikelihood_Tests, set_binned_likelihood_1) {
	TFeatureLikelihood storage(0, 255);
	EXPECT_EQ(storage.get_binned_likelihood(), 255);
	EXPECT_EQ(storage.get_molecule_index(), 0);
}

TEST(TFeatureLikelihood_Tests, set_binned_likelihood_2) {
	TFeatureLikelihood storage(0, 0);
	EXPECT_EQ(storage.get_binned_likelihood(), 0);
}

TEST(TFeatureLikelihood_Tests, set_binned_likelihood_3) {
	TFeatureLikelihood storage(0, 128);
	EXPECT_EQ(storage.get_binned_likelihood(), 128);
	EXPECT_EQ(storage.get_molecule_index(), 0);
}

TEST(TFeatureLikelihood_Tests, combined) {
	TFeatureLikelihood storage(16777214, 128);
	EXPECT_EQ(storage.get_binned_likelihood(), 128);
	EXPECT_EQ(storage.get_molecule_index(), 16777214);
}

TEST(TFeatureLikelihood_Tests, nested_vector) {
	TFeatureLikelihood fl_1(16777214, 128);

	TFeatureLikelihood fl_2(245, 120);

	std::vector<TFeatureLikelihood> vec_1 = {fl_1};
	std::vector<TFeatureLikelihood> vec_2 = {fl_1, fl_2};
	TMassSpecRun nested_vector;
	nested_vector.add_likelihood_vector(vec_1);
	nested_vector.add_likelihood_vector(vec_2);

	EXPECT_EQ(nested_vector.get_likelihoods_for_feature(0).size(), 1);
	EXPECT_EQ(nested_vector.get_likelihoods_for_feature(1).size(), 2);

	EXPECT_EQ(nested_vector.get_likelihoods_for_feature(1).at(0).get_binned_likelihood(), 120);

	coretools::TNestedVector<TMassSpecRun> msms_data;
	msms_data.push_back({nested_vector});
	msms_data.push_back({nested_vector, nested_vector, nested_vector});
	msms_data.push_back({});
	EXPECT_EQ(msms_data.at(0).size(), 1);
	EXPECT_EQ(msms_data.at(1).size(), 3);
	EXPECT_EQ(msms_data.at(2).size(), 0);
	EXPECT_THROW(msms_data.at(3), std::out_of_range);
}
