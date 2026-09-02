//
// Where a chain's leaf layer starts, over every pairing of the two storages.
//
// The start holds no tree and draws nothing. It reads the LOTUS records and writes the field and
// both tree fields, so it is asserted against storages alone. The shapes and the pairings come
// from backend_pairings.h.
//

#include "constants.h"
#include "backend_pairings.h"
#include "field/leaf_layer_start.h"
#include "field/link_backend.h"
#include "gtest/gtest.h"

#include <cstddef>
#include <string>

namespace {

using backends::AllBackends;
using backends::field_shape;
using backends::make_storage;
using backends::molecule_shape;
using backends::seed_ones;
using backends::species_shape;
using backends::tree_pairs;

/// The records the start reads: every third cell, so that the smallest container space the tree
/// pairs ask for still holds both a record and a cell without one. The sparse backends hold only
/// what is put in them, so a cell left out is one the start has to insert rather than write in
/// place.
template<typename Field> void seed_records(Field &records) {
	for (size_t i = 0; i < records.total_size_of_container_space(); ++i) {
		if (i % 3U == 0U) { records.insert_one(i); }
	}
}

/// coretools' dev error is a CTAD class template, so it cannot be named in an exception
/// declaration and EXPECT_THROW cannot be used on it. Catching the err::TError base and asserting
/// it *is* a dev error is what "this is a programming error" means here.
template<typename F> void expect_dev_error(F &&f) {
	try {
		f();
	} catch (coretools::err::TError &e) {
		EXPECT_TRUE(e.isDevError()) << "expected a dev error, got a user error: " << e.what();
		return;
	}
	ADD_FAILURE() << "expected a dev error, but nothing was thrown";
}

/// Whether the records hold this leaf pair, read the way the field's own index reads it.
template<typename Field> bool records_hold(const Field &records, size_t s, size_t m) {
	return records.is_one(records.get_linear_index_in_container_space(IndexArray{s, m}));
}

// -------------------------------------------------------------------------
// The suite, over all four storage pairings
// -------------------------------------------------------------------------

template<typename Backends> class LeafLayerStart : public ::testing::Test {
public:
	using Field     = typename Backends::field;
	using NodeState = typename Backends::node_state;
};

TYPED_TEST_SUITE(LeafLayerStart, AllBackends);

/// The field reads one exactly where a record exists, and both tree fields read the field.
TYPED_TEST(LeafLayerStart, starts_all_three_at_the_records) {
	using Field     = typename TestFixture::Field;
	using NodeState = typename TestFixture::NodeState;

	for (const auto &pair : tree_pairs()) {
		SCOPED_TRACE(pair.name);
		auto records    = make_storage<Field>(field_shape(pair));
		auto Y          = make_storage<Field>(field_shape(pair));
		auto Z_species  = make_storage<NodeState>(species_shape(pair));
		auto Z_molecule = make_storage<NodeState>(molecule_shape(pair));
		seed_records(records);
		// States the start has to overwrite: a tree field that read one where no record does is
		// the case a start that only ever inserts ones would leave behind.
		seed_ones(Z_species, 2);
		seed_ones(Z_molecule, 3);

		leaf_layer_start::start_the_field_at(records, Y);
		const auto counters = leaf_layer_start::hold_tree_fields_at_the_field<TLinkPolicy>(
		    Y, Z_species, Z_molecule);
		EXPECT_EQ(counters.total(), Y.total_size_of_container_space());

		for (size_t s = 0; s < pair.species.n_leaves(); ++s) {
			for (size_t m = 0; m < pair.molecule.n_leaves(); ++m) {
				SCOPED_TRACE("leaf pair " + std::to_string(s) + "," + std::to_string(m));
				const bool reported = records_hold(records, s, m);
				EXPECT_EQ(Y.is_one(Y.get_linear_index_in_container_space(IndexArray{s, m})),
				          reported);
				EXPECT_EQ(Z_species.is_one(
				              Z_species.get_linear_index_in_container_space(IndexArray{s, m})),
				          reported);
				EXPECT_EQ(Z_molecule.is_one(
				              Z_molecule.get_linear_index_in_container_space(IndexArray{s, m})),
				          reported);
			}
		}
	}
}

/// With no record anywhere the start is all zeros -- the field, both tree fields, and whatever
/// they held before it ran.
TYPED_TEST(LeafLayerStart, with_no_records_everything_starts_at_zero) {
	using Field     = typename TestFixture::Field;
	using NodeState = typename TestFixture::NodeState;

	for (const auto &pair : tree_pairs()) {
		SCOPED_TRACE(pair.name);
		const auto records = make_storage<Field>(field_shape(pair));
		auto Y             = make_storage<Field>(field_shape(pair));
		auto Z_species     = make_storage<NodeState>(species_shape(pair));
		auto Z_molecule    = make_storage<NodeState>(molecule_shape(pair));
		seed_ones(Z_species, 2);
		seed_ones(Z_molecule, 3);

		leaf_layer_start::start_the_field_at(records, Y);
		const auto counters = leaf_layer_start::hold_tree_fields_at_the_field<TLinkPolicy>(
		    Y, Z_species, Z_molecule);

		for (size_t s = 0; s < pair.species.n_leaves(); ++s) {
			for (size_t m = 0; m < pair.molecule.n_leaves(); ++m) {
				SCOPED_TRACE("leaf pair " + std::to_string(s) + "," + std::to_string(m));
				EXPECT_FALSE(Y.is_one(Y.get_linear_index_in_container_space(IndexArray{s, m})));
				EXPECT_FALSE(Z_species.is_one(
				    Z_species.get_linear_index_in_container_space(IndexArray{s, m})));
				EXPECT_FALSE(Z_molecule.is_one(
				    Z_molecule.get_linear_index_in_container_space(IndexArray{s, m})));
			}
		}
		// Every leaf pair is a zero in bucket 0, which is the degenerate tally in its emptiest
		// form.
		EXPECT_EQ(counters.count(0, false), Y.total_size_of_container_space());
		EXPECT_EQ(counters.total(), Y.total_size_of_container_space());
	}
}

/// The six counters the start leaves behind are degenerate: a record puts a leaf pair in bucket 2
/// with the field at one, and no record puts it in bucket 0 with the field at zero. Bucket 1 holds
/// nothing, so the AND diagnostic can say nothing about them. The first block update replaces the
/// whole tally.
TYPED_TEST(LeafLayerStart, tallies_the_counters_it_leaves_degenerate) {
	using Field     = typename TestFixture::Field;
	using NodeState = typename TestFixture::NodeState;

	for (const auto &pair : tree_pairs()) {
		SCOPED_TRACE(pair.name);
		auto records    = make_storage<Field>(field_shape(pair));
		auto Y          = make_storage<Field>(field_shape(pair));
		auto Z_species  = make_storage<NodeState>(species_shape(pair));
		auto Z_molecule = make_storage<NodeState>(molecule_shape(pair));
		seed_records(records);

		leaf_layer_start::start_the_field_at(records, Y);
		const auto counters = leaf_layer_start::hold_tree_fields_at_the_field<TLinkPolicy>(
		    Y, Z_species, Z_molecule);

		size_t n_reported = 0;
		for (size_t s = 0; s < pair.species.n_leaves(); ++s) {
			for (size_t m = 0; m < pair.molecule.n_leaves(); ++m) {
				n_reported += static_cast<size_t>(records_hold(records, s, m));
			}
		}
		const size_t n_cells = Y.total_size_of_container_space();
		ASSERT_GT(n_reported, 0U); // otherwise the case below is the empty one, not this one

		EXPECT_EQ(counters.count(2, true), n_reported);
		EXPECT_EQ(counters.count(0, false), n_cells - n_reported);
		EXPECT_EQ(counters.count(2, false), 0U);
		EXPECT_EQ(counters.count(0, true), 0U);
		EXPECT_EQ(counters.count(1, false), 0U);
		EXPECT_EQ(counters.count(1, true), 0U);
	}
}

/// Records of another shape are a programming error, not a state the start makes the best of.
TYPED_TEST(LeafLayerStart, refuses_records_of_another_shape) {
	using Field = typename TestFixture::Field;

	const auto &pair = tree_pairs().front();
	const auto records =
	    make_storage<Field>(IndexArray{pair.species.n_leaves() + 1, pair.molecule.n_leaves()});
	auto Y = make_storage<Field>(field_shape(pair));

	expect_dev_error([&] { leaf_layer_start::start_the_field_at(records, Y); });
}

/// The start writes ones and never zeros, so it reads one *exactly* where a record exists only on
/// an empty field. A field that already holds states is refused rather than merged into.
TYPED_TEST(LeafLayerStart, refuses_a_field_that_already_holds_states) {
	using Field = typename TestFixture::Field;

	const auto &pair = tree_pairs().front();
	auto records     = make_storage<Field>(field_shape(pair));
	auto Y           = make_storage<Field>(field_shape(pair));
	seed_records(records);
	seed_ones(Y, 7);

	expect_dev_error([&] { leaf_layer_start::start_the_field_at(records, Y); });
}

} // namespace
