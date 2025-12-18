#include "TMassSpecStorage.h"
#include "gtest/gtest.h"

TEST(TFeatureLikelihood_Tests, set_linear_index_1) {
	TFeatureLikelihood storage;
	storage.set_molecule_index(16777215);
	EXPECT_EQ(storage.get_molecule_index(), 16777215);
}

TEST(TFeatureLikelihood_Tests, set_linear_index_2) {
	TFeatureLikelihood storage;
	storage.set_molecule_index(0);
	EXPECT_EQ(storage.get_molecule_index(), 0);
}

TEST(TFeatureLikelihood_Tests, set_linear_index_3) {
	TFeatureLikelihood storage;
	storage.set_molecule_index(1);
	EXPECT_EQ(storage.get_molecule_index(), 1);
}

TEST(TFeatureLikelihood_Tests, set_binned_likelihood_1) {
	TFeatureLikelihood storage;
	storage.set_binned_likelihood(255);
	EXPECT_EQ(storage.get_binned_likelihood(), 255);
}

TEST(TFeatureLikelihood_Tests, set_binned_likelihood_2) {
	TFeatureLikelihood storage;
	storage.set_binned_likelihood(0);
	EXPECT_EQ(storage.get_binned_likelihood(), 0);
}

TEST(TFeatureLikelihood_Tests, set_binned_likelihood_3) {
	TFeatureLikelihood storage;
	storage.set_binned_likelihood(128);
	EXPECT_EQ(storage.get_binned_likelihood(), 128);
}

TEST(TFeatureLikelihood_Tests, combined) {
	TFeatureLikelihood storage;
	storage.set_binned_likelihood(128);
	storage.set_molecule_index(16777215);
	EXPECT_EQ(storage.get_binned_likelihood(), 128);
	EXPECT_EQ(storage.get_molecule_index(), 16777215);
}

TEST(TFeatureLikelihood_Tests, nested_vector) {
	TFeatureLikelihood fl_1;
	fl_1.set_binned_likelihood(128);
	fl_1.set_molecule_index(16777215);

	TFeatureLikelihood fl_2;
	fl_2.set_binned_likelihood(120);
	fl_2.set_molecule_index(245);

	std::vector<TFeatureLikelihood> vec_1 = {fl_1};
	std::vector<TFeatureLikelihood> vec_2 = {fl_1, fl_2};
	TMassSpecRun nested_vector;
	nested_vector.add_likelihood_vector(vec_1);
	nested_vector.add_likelihood_vector(vec_2);

	EXPECT_EQ(nested_vector.get_likelihoods_for_feature(0).size(), 1);
	EXPECT_EQ(nested_vector.get_likelihoods_for_feature(1).size(), 2);
	for (const auto &fl : nested_vector) {
		for (const auto &f : fl) {
			std::cout << f.get_molecule_index() << " " << f.get_binned_likelihood() << std::endl;
		}
	};

	EXPECT_EQ(nested_vector.get_likelihoods_for_feature(1).at(0).get_binned_likelihood(), 120);
}
