//
// The field a simulated chain draws, over every pairing of the two storages.
//
// The draw reads the two tree fields and the link, and writes the field and the six counters. It
// holds no tree, so it is asserted against storages alone. The shapes and the pairings come from
// backend_pairings.h.
//
// Two properties carry the model. Every leaf pair is drawn from `P(Y = 1 | Z_s, Z_m, omega)` and
// from nothing else, which is checked against the link table cell by cell. And the counters that
// come back are the tally of the field the draw left behind, which is checked by recounting it --
// the counters are what the error probability's whole likelihood reads, so a tally that drifts
// from the field would be silent (ADR-0005).
//

#include "backend_pairings.h"
#include "constants.h"
#include "field/link_backend.h"
#include "field/simulate_field.h"
#include "written_uniforms.h"
#include "gtest/gtest.h"

#include <cmath>
#include <cstddef>
#include <random>
#include <vector>

namespace {

using backends::AllBackends;
using backends::field_shape;
using backends::make_storage;
using backends::molecule_shape;
using backends::seed_ones;
using backends::species_shape;
using backends::tree_pairs;
using simulate_field::draw_from_the_tree_fields;

/// coretools' dev error is a CTAD class template, so it cannot be named in an exception
/// declaration and EXPECT_THROW cannot be used on it.
template<typename F> void expect_dev_error(F &&f) {
	try {
		f();
	} catch (coretools::err::TError &e) {
		EXPECT_TRUE(e.isDevError()) << "expected a dev error, got a user error: " << e.what();
		return;
	}
	ADD_FAILURE() << "expected a dev error, but nothing was thrown";
}

/// One tree field cell, read the way the field's own index reads it.
template<typename NodeState>
bool tree_field_at(const NodeState &tree_field, const IndexArray &cell) {
	return tree_field.is_one(tree_field.get_linear_index_in_container_space(cell));
}

// -------------------------------------------------------------------------
// The suite, over all four storage pairings
// -------------------------------------------------------------------------

template<typename Backends> class SimulateField : public ::testing::Test {
public:
	using Field     = typename Backends::field;
	using NodeState = typename Backends::node_state;
};

TYPED_TEST_SUITE(SimulateField, AllBackends);

/// Every leaf pair reads one exactly when its own uniform falls below the link's probability for
/// the two tree field states there.
///
/// The threshold is written out from the error probability rather than taken from the policy the
/// draw itself calls, so this asserts the loop -- every cell visited once, against the two tree
/// field cells at its own subscript -- and not the table. The table has its own suite.
TYPED_TEST(SimulateField, draws_each_cell_from_the_link_and_from_nothing_else) {
	using Field                        = typename TestFixture::Field;
	using NodeState                    = typename TestFixture::NodeState;
	constexpr double ERROR_PROBABILITY = 0.1;
	const field_math::TErrorProbability omega(ERROR_PROBABILITY);
	const auto corrupted_read = [](bool state) {
		return state ? 1.0 - ERROR_PROBABILITY : ERROR_PROBABILITY;
	};

	for (const auto &pair : tree_pairs()) {
		SCOPED_TRACE(pair.name);
		auto Y          = make_storage<Field>(field_shape(pair));
		auto Z_species  = make_storage<NodeState>(species_shape(pair));
		auto Z_molecule = make_storage<NodeState>(molecule_shape(pair));
		seed_ones(Z_species, 1);
		seed_ones(Z_molecule, 2);

		uniforms::TWrittenUniforms uniforms(Y.total_size_of_container_space());
		std::mt19937_64 rng(20260907);
		uniforms.fill_from(rng);

		const auto counters =
		    draw_from_the_tree_fields<TLinkPolicy>(Y, Z_species, Z_molecule, omega, uniforms);

		for (size_t index = 0; index < Y.total_size_of_container_space(); ++index) {
			const IndexArray cell = Y.get_multi_dimensional_index(index);
			const bool z_s        = tree_field_at(Z_species, cell);
			const bool z_m        = tree_field_at(Z_molecule, cell);
			const double expected = corrupted_read(z_s) * corrupted_read(z_m);
			EXPECT_EQ(Y.is_one(index), uniforms.at(index) < expected) << "cell " << index;
		}
		EXPECT_EQ(counters.total(), Y.total_size_of_container_space());
	}
}

/// The counters that come back are the tally of the field the draw left behind. They are what the
/// error probability's likelihood reads, so a tally that drifts from the field would be silent.
TYPED_TEST(SimulateField, hands_back_the_counters_of_the_field_it_drew) {
	using Field     = typename TestFixture::Field;
	using NodeState = typename TestFixture::NodeState;
	const field_math::TErrorProbability omega(0.2);

	for (const auto &pair : tree_pairs()) {
		SCOPED_TRACE(pair.name);
		auto Y          = make_storage<Field>(field_shape(pair));
		auto Z_species  = make_storage<NodeState>(species_shape(pair));
		auto Z_molecule = make_storage<NodeState>(molecule_shape(pair));
		seed_ones(Z_species, 3);
		seed_ones(Z_molecule, 4);

		uniforms::TWrittenUniforms uniforms(Y.total_size_of_container_space());
		std::mt19937_64 rng(20260908);
		uniforms.fill_from(rng);

		const auto counters =
		    draw_from_the_tree_fields<TLinkPolicy>(Y, Z_species, Z_molecule, omega, uniforms);

		field_math::TLinkCounters recounted;
		for (size_t index = 0; index < Y.total_size_of_container_space(); ++index) {
			const IndexArray cell = Y.get_multi_dimensional_index(index);
			recounted.add(TLinkPolicy::bucket(tree_field_at(Z_species, cell),
			                                  tree_field_at(Z_molecule, cell)),
			              Y.is_one(index));
		}
		for (size_t bucket = 0; bucket < field_math::TLinkCounters::n_buckets; ++bucket) {
			for (const bool y : {false, true}) {
				EXPECT_EQ(counters.count(bucket, y), recounted.count(bucket, y))
				    << "bucket " << bucket << ", field state " << y;
			}
		}
	}
}

/// The AND, at the two ends of the error probability's range. Just above 0 the field is the AND of
/// the two tree fields; just below 0.5 it is a coin flip whatever they say.
TYPED_TEST(SimulateField, is_the_and_of_the_tree_fields_when_the_error_probability_vanishes) {
	using Field     = typename TestFixture::Field;
	using NodeState = typename TestFixture::NodeState;

	for (const auto &pair : tree_pairs()) {
		SCOPED_TRACE(pair.name);
		auto Y          = make_storage<Field>(field_shape(pair));
		auto Z_species  = make_storage<NodeState>(species_shape(pair));
		auto Z_molecule = make_storage<NodeState>(molecule_shape(pair));
		seed_ones(Z_species, 5);
		seed_ones(Z_molecule, 6);

		// Mid-range uniforms: at an error probability of 1e-9 the link's four probabilities are
		// either within 1e-9 of 1 or within 1e-9 of 0, so 0.5 lands on the AND at every cell.
		const uniforms::TWrittenUniforms uniforms(Y.total_size_of_container_space(), 0.5);
		const auto counters = draw_from_the_tree_fields<TLinkPolicy>(
		    Y, Z_species, Z_molecule, field_math::TErrorProbability(1e-9), uniforms);

		for (size_t index = 0; index < Y.total_size_of_container_space(); ++index) {
			const IndexArray cell = Y.get_multi_dimensional_index(index);
			EXPECT_EQ(Y.is_one(index),
			          tree_field_at(Z_species, cell) && tree_field_at(Z_molecule, cell))
			    << "cell " << index;
		}
		// Only bucket 2 holds a one, which is what the AND says.
		EXPECT_EQ(counters.count(0, true), 0u);
		EXPECT_EQ(counters.count(1, true), 0u);
	}
}

/// The field a draw writes into must be empty, or the states it already held would survive a draw
/// that never writes a zero.
TYPED_TEST(SimulateField, refuses_a_field_that_already_holds_states) {
	using Field     = typename TestFixture::Field;
	using NodeState = typename TestFixture::NodeState;

	const auto &pair = tree_pairs().front();
	auto Y           = make_storage<Field>(field_shape(pair));
	auto Z_species   = make_storage<NodeState>(species_shape(pair));
	auto Z_molecule  = make_storage<NodeState>(molecule_shape(pair));
	Y.insert_one(0);

	const uniforms::TWrittenUniforms uniforms(Y.total_size_of_container_space(), 0.5);
	expect_dev_error([&] {
		(void)draw_from_the_tree_fields<TLinkPolicy>(Y, Z_species, Z_molecule,
		                                             field_math::TErrorProbability(0.1), uniforms);
	});
}

} // namespace
