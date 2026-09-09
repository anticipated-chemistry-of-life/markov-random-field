//
// The field update's loop, over every pairing of the two storages.
//
// The link is pure and tested against brute force (TFieldMath_Tests.cpp), and the storages are
// conformance-tested against each other (TStorageConformance_Tests.cpp). What is left belongs to
// the loop alone, and that is what this file asserts: every cell visited exactly once, the two tree
// field cells read being the ones at that leaf pair, the data terms landing on the cell they were
// scored for, the write reaching the storage, the same chain whatever the thread count, and a tally
// that matches a naive recount of the configuration the pass left.
//
// The loop is asked those questions through a model of its own rather than through the data
// sources, because a model the test writes can be *driven*: a data likelihood of zero for one field
// state leaves the draw one state to take. Every write is then a known value at a known cell, and a
// write that lands on the wrong cell is a wrong value rather than a coincidence. A stub model also
// runs in a build that compiled no data source in.
//
// Every body is instantiated over all four field/node-state pairings. Continuous integration gates
// two of them (`just parity`), so this is where the other two are exercised at all.
//

#include "backend_pairings.h"
#include "cli.h"
#include "constants.h"
#include "coretools/Main/TError.h"
#include "coretools/Types/probability.h"
#include "field/TFieldMath.h"
#include "field/field_update.h"
#include "field/link_backend.h"
#include "random/TCellUniforms.h"
#include "storages/y_storage/TStorageYSparse.h"
#include "storages/z_storage/TStorageZSparse.h"
#include "tree/TPhylogeny.h"
#include "written_uniforms.h"
#include "gtest/gtest.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

namespace {

/// Small enough that the link never runs degenerate, and recognisable in a failure message.
constexpr double OMEGA = 0.125;

using backends::AllBackends;
using backends::field_shape;
using backends::make_storage;
using backends::molecule_shape;
using backends::seed_ones;
using backends::species_shape;
using backends::tree_pairs;
using backends::TTreePair;

/// Every cell of a container space, as states, so two runs can be compared without either backend
/// having to say which cells it holds.
template<typename Storage> std::vector<uint8_t> states_of(const Storage &storage) {
	std::vector<uint8_t> states;
	states.reserve(storage.total_size_of_container_space());
	for (size_t i = 0; i < storage.total_size_of_container_space(); ++i) {
		states.push_back(static_cast<uint8_t>(storage.is_one(i)));
	}
	return states;
}

// -------------------------------------------------------------------------
// The models the loop is run against
// -------------------------------------------------------------------------

/// The field state one cell is driven to. An arbitrary but fixed function of the leaf pair, so that
/// a write landing on the wrong cell writes the wrong value.
bool target_at(size_t species_leaf, size_t molecule_leaf) {
	const uint64_t bits = (6364136223846793005ULL * (species_leaf + 1)) ^
	                      (1442695040888963407ULL * (molecule_leaf + 3));
	return ((bits >> 17U) & 1U) != 0U;
}

/// A model that drives every cell to `target_at`, and records what the loop asked it.
///
/// The forcing is what makes the recording checkable. A data likelihood of 0 for one field state
/// leaves that state with no mass, so the draw has one state left to take whatever uniform it is
/// given.
class TForcingModel {
public:
	/// What the loop did at one cell.
	struct TVisit {
		size_t n_asked        = 0; ///< how often `factors` was asked about this cell
		size_t n_recorded     = 0; ///< how often `record` was told about it
		/// Whether `record` was handed back the factors `factors` scored this cell with.
		bool kept_its_factors = false;
		bool drawn            = false;
	};

	TForcingModel(size_t n_species_leaves, size_t n_molecule_leaves)
	    : _n_molecule_leaves(n_molecule_leaves), _visits(n_species_leaves * n_molecule_leaves) {}

	[[nodiscard]] field_update::TCellFactors factors(const IndexArray &cell) {
		TVisit &recorded = visit(cell[0], cell[1]);
		++recorded.n_asked;

		const bool target = target_at(cell[0], cell[1]);
		return {.lotus = {coretools::P(target ? 0.0 : 1.0), coretools::P(target ? 1.0 : 0.0)},
		        .simple_error = {coretools::P(1.0), coretools::P(1.0)}};
	}

	void record(const IndexArray &cell, const field_update::TCellFactors &factors, bool drawn) {
		TVisit &recorded = visit(cell[0], cell[1]);
		++recorded.n_recorded;
		recorded.drawn            = drawn;
		const bool target         = target_at(cell[0], cell[1]);
		recorded.kept_its_factors = factors.lotus[1].get() == (target ? 1.0 : 0.0);
	}

	[[nodiscard]] TVisit &visit(size_t species_leaf, size_t molecule_leaf) {
		return _visits[species_leaf * _n_molecule_leaves + molecule_leaf];
	}

private:
	size_t _n_molecule_leaves;
	std::vector<TVisit> _visits;
};

static_assert(field_update::FieldModel<TForcingModel>,
              "The forcing model must answer what the field update asks a model.");

/// A model that says nothing about any cell, so the draw comes from the link alone.
class TNeutralModel {
public:
	[[nodiscard]] static field_update::TCellFactors factors(const IndexArray &) { return {}; }
	static void record(const IndexArray &, const field_update::TCellFactors &, bool) {}
};

static_assert(field_update::FieldModel<TNeutralModel>,
              "The neutral model must answer what the field update asks a model.");

/// A model that leaves the draw a real choice at every cell, and keeps no state of its own.
///
/// The forcing model above pins every cell, which would let a wrong uniform pass unnoticed. Here
/// both field states carry mass, so the state a cell ends in depends on the uniform it drew --
/// which is what makes a chain comparable between two thread counts. Being stateless is what makes
/// it safe to run on many threads.
class TFreeModel {
public:
	[[nodiscard]] static field_update::TCellFactors factors(const IndexArray &cell) {
		// A value that depends on the whole leaf pair, so that a read of the wrong cell moves the
		// chain.
		const double drift = 0.1 * static_cast<double>((cell[0] + cell[1]) % 4U);
		return {.lotus        = {coretools::P(0.4), coretools::P(0.6 - drift)},
		        .simple_error = {coretools::P(0.55), coretools::P(0.45 + drift)}};
	}

	static void record(const IndexArray &, const field_update::TCellFactors &, bool) {}
};

static_assert(field_update::FieldModel<TFreeModel>,
              "The free model must answer what the field update asks a model.");

/// The six counters recomputed from a whole configuration, sharing nothing with the tally the pass
/// kept as it went.
template<typename Field, typename NodeState>
field_math::TLinkCounters recount(const Field &Y, const NodeState &Z_species,
                                  const NodeState &Z_molecule, const TTreePair &pair) {
	field_math::TLinkCounters counters;
	for (size_t s = 0; s < pair.species.n_leaves(); ++s) {
		for (size_t m = 0; m < pair.molecule.n_leaves(); ++m) {
			const bool z_s =
			    Z_species.is_one(Z_species.get_linear_index_in_container_space({s, m}));
			const bool z_m =
			    Z_molecule.is_one(Z_molecule.get_linear_index_in_container_space({s, m}));
			counters.add(TLinkPolicy::bucket(z_s, z_m),
			             Y.is_one(Y.get_linear_index_in_container_space({s, m})));
		}
	}
	return counters;
}

/// The tallies of one pass, merged the way the caller of a field update merges them.
field_math::TLinkCounters merged(const std::vector<field_math::TLinkCounters> &tallies) {
	field_math::TLinkCounters counters;
	for (const auto &tally : tallies) { counters.merge(tally); }
	return counters;
}

/// Sets the thread count for one test and puts it back afterwards. The count is a global, and a
/// test that left it raised would change what every later test runs on.
class TThreadCount {
private:
	size_t _previous;

public:
	explicit TThreadCount(size_t n_threads) : _previous(ProgramOptions::NUMBER_OF_THREADS) {
		ProgramOptions::NUMBER_OF_THREADS = n_threads;
	}
	~TThreadCount() { ProgramOptions::NUMBER_OF_THREADS = _previous; }

	TThreadCount(const TThreadCount &)            = delete;
	TThreadCount &operator=(const TThreadCount &) = delete;
	TThreadCount(TThreadCount &&)                 = delete;
	TThreadCount &operator=(TThreadCount &&)      = delete;
};

// -------------------------------------------------------------------------
// The suite, over all four storage pairings
// -------------------------------------------------------------------------

template<typename Backends> class FieldUpdate : public ::testing::Test {
public:
	using Field     = typename Backends::field;
	using NodeState = typename Backends::node_state;
};

TYPED_TEST_SUITE(FieldUpdate, AllBackends);

/// Every cell is asked about once and told what it was given once, and no cell is missed. Each is
/// handed back the factors it was scored with, which is what lets a data source keep its own
/// likelihood bookkeeping.
TYPED_TEST(FieldUpdate, visits_every_cell_exactly_once) {
	using Field     = typename TestFixture::Field;
	using NodeState = typename TestFixture::NodeState;

	for (const auto &pair : tree_pairs()) {
		SCOPED_TRACE(pair.name);
		auto Y          = make_storage<Field>(field_shape(pair));
		auto Z_species  = make_storage<NodeState>(species_shape(pair));
		auto Z_molecule = make_storage<NodeState>(molecule_shape(pair));
		seed_ones(Y, 1);
		seed_ones(Z_species, 2);
		seed_ones(Z_molecule, 3);

		TForcingModel model(pair.species.n_leaves(), pair.molecule.n_leaves());
		std::vector<field_math::TLinkCounters> tallies(ProgramOptions::NUMBER_OF_THREADS);
		const TCellUniforms uniforms(4242, TCellStream::field, 0);
		field_update::run<TLinkPolicy>(Y, Z_species, Z_molecule,
		                               field_math::TErrorProbability(OMEGA), model, uniforms,
		                               tallies);

		for (size_t s = 0; s < pair.species.n_leaves(); ++s) {
			for (size_t m = 0; m < pair.molecule.n_leaves(); ++m) {
				SCOPED_TRACE("leaf pair " + std::to_string(s) + "," + std::to_string(m));
				EXPECT_EQ(model.visit(s, m).n_asked, 1U);
				EXPECT_EQ(model.visit(s, m).n_recorded, 1U);
				EXPECT_TRUE(model.visit(s, m).kept_its_factors);
			}
		}
	}
}

/// The state the draw assigned is in the field afterwards, at the cell it was drawn for -- through
/// an in-place write where the backend held the cell, and through the deferred insert where it did
/// not.
TYPED_TEST(FieldUpdate, writes_the_drawn_state_back) {
	using Field     = typename TestFixture::Field;
	using NodeState = typename TestFixture::NodeState;

	for (const auto &pair : tree_pairs()) {
		SCOPED_TRACE(pair.name);
		auto Y          = make_storage<Field>(field_shape(pair));
		auto Z_species  = make_storage<NodeState>(species_shape(pair));
		auto Z_molecule = make_storage<NodeState>(molecule_shape(pair));
		seed_ones(Y, 1);
		seed_ones(Z_species, 2);
		seed_ones(Z_molecule, 3);

		const std::vector<uint8_t> species_before  = states_of(Z_species);
		const std::vector<uint8_t> molecule_before = states_of(Z_molecule);

		TForcingModel model(pair.species.n_leaves(), pair.molecule.n_leaves());
		std::vector<field_math::TLinkCounters> tallies(ProgramOptions::NUMBER_OF_THREADS);
		const TCellUniforms uniforms(4242, TCellStream::field, 0);
		field_update::run<TLinkPolicy>(Y, Z_species, Z_molecule,
		                               field_math::TErrorProbability(OMEGA), model, uniforms,
		                               tallies);

		for (size_t s = 0; s < pair.species.n_leaves(); ++s) {
			for (size_t m = 0; m < pair.molecule.n_leaves(); ++m) {
				SCOPED_TRACE("leaf pair " + std::to_string(s) + "," + std::to_string(m));
				// the draw had one state to take, so the model was driven to the target
				EXPECT_EQ(model.visit(s, m).drawn, target_at(s, m));
				// and the target is what the field holds at that cell
				EXPECT_EQ(Y.is_one(IndexArray{s, m}), target_at(s, m));
			}
		}

		// The pass draws the field alone. Both tree fields are their own tree's to draw.
		EXPECT_EQ(states_of(Z_species), species_before);
		EXPECT_EQ(states_of(Z_molecule), molecule_before);
	}
}

/// The two tree field cells the loop reads are the ones at its own leaf pair.
///
/// With no data source saying anything and a uniform of exactly one half, the link decides the
/// cell: two tree fields at one put the field at one, and any other pair puts it at zero. So the
/// field the pass leaves is the AND of the two tree fields, and a read of the wrong cell shows up
/// as a wrong state.
TYPED_TEST(FieldUpdate, reads_the_tree_field_cells_of_its_own_leaf_pair) {
	using Field     = typename TestFixture::Field;
	using NodeState = typename TestFixture::NodeState;

	for (const auto &pair : tree_pairs()) {
		SCOPED_TRACE(pair.name);
		auto Y          = make_storage<Field>(field_shape(pair));
		auto Z_species  = make_storage<NodeState>(species_shape(pair));
		auto Z_molecule = make_storage<NodeState>(molecule_shape(pair));
		seed_ones(Y, 1);
		seed_ones(Z_species, 2);
		seed_ones(Z_molecule, 3);

		// (1 - omega)^2 is above one half and both other buckets are below it, so one half splits
		// the three buckets exactly where the AND does.
		ASSERT_GT(TLinkPolicy::prob_for_bucket(2, field_math::TErrorProbability(OMEGA)), 0.5);
		ASSERT_LT(TLinkPolicy::prob_for_bucket(1, field_math::TErrorProbability(OMEGA)), 0.5);

		TNeutralModel model;
		std::vector<field_math::TLinkCounters> tallies(ProgramOptions::NUMBER_OF_THREADS);
		const uniforms::TWrittenUniforms uniforms(Y.total_size_of_container_space(), 0.5);
		field_update::run<TLinkPolicy>(Y, Z_species, Z_molecule,
		                               field_math::TErrorProbability(OMEGA), model, uniforms,
		                               tallies);

		for (size_t s = 0; s < pair.species.n_leaves(); ++s) {
			for (size_t m = 0; m < pair.molecule.n_leaves(); ++m) {
				SCOPED_TRACE("leaf pair " + std::to_string(s) + "," + std::to_string(m));
				const bool z_s =
				    Z_species.is_one(Z_species.get_linear_index_in_container_space({s, m}));
				const bool z_m =
				    Z_molecule.is_one(Z_molecule.get_linear_index_in_container_space({s, m}));
				EXPECT_EQ(Y.is_one(IndexArray{s, m}), z_s && z_m);
			}
		}
	}
}

/// One thread and many give the same field and the same six counters. A cell's uniform is hashed
/// from its position (ADR-0007), so the thread that reaches it does not decide what it gets.
TYPED_TEST(FieldUpdate, gives_the_same_chain_at_any_thread_count) {
	using Field     = typename TestFixture::Field;
	using NodeState = typename TestFixture::NodeState;

	for (const auto &pair : tree_pairs()) {
		SCOPED_TRACE(pair.name);

		// Both runs start from the same configuration and draw from the same stream.
		const auto run_once = [&pair](size_t n_threads) {
			const TThreadCount threads(n_threads);
			auto Y          = make_storage<Field>(field_shape(pair));
			auto Z_species  = make_storage<NodeState>(species_shape(pair));
			auto Z_molecule = make_storage<NodeState>(molecule_shape(pair));
			seed_ones(Y, 1);
			seed_ones(Z_species, 2);
			seed_ones(Z_molecule, 3);

			TFreeModel model;
			std::vector<field_math::TLinkCounters> tallies(n_threads);
			const TCellUniforms uniforms(4242, TCellStream::field, 7);
			field_update::run<TLinkPolicy>(Y, Z_species, Z_molecule,
			                               field_math::TErrorProbability(OMEGA), model, uniforms,
			                               tallies);
			return std::tuple{states_of(Y), merged(tallies)};
		};

		const auto [one_Y, one_counters]   = run_once(1);
		const auto [many_Y, many_counters] = run_once(4);

		EXPECT_EQ(one_Y, many_Y);
		for (size_t bucket = 0; bucket < field_math::TLinkCounters::n_buckets; ++bucket) {
			for (const bool y : {false, true}) {
				EXPECT_EQ(one_counters.count(bucket, y), many_counters.count(bucket, y));
			}
		}
	}
}

/// The counters the pass accumulated are the tally of the configuration it left behind, and they
/// count every cell once.
TYPED_TEST(FieldUpdate, counters_tally_the_configuration_it_left) {
	using Field     = typename TestFixture::Field;
	using NodeState = typename TestFixture::NodeState;

	for (const auto &pair : tree_pairs()) {
		SCOPED_TRACE(pair.name);
		const TThreadCount threads(3);
		auto Y          = make_storage<Field>(field_shape(pair));
		auto Z_species  = make_storage<NodeState>(species_shape(pair));
		auto Z_molecule = make_storage<NodeState>(molecule_shape(pair));
		seed_ones(Y, 1);
		seed_ones(Z_species, 2);
		seed_ones(Z_molecule, 3);

		TFreeModel model;
		std::vector<field_math::TLinkCounters> tallies(ProgramOptions::NUMBER_OF_THREADS);
		const TCellUniforms uniforms(4242, TCellStream::field, 11);
		field_update::run<TLinkPolicy>(Y, Z_species, Z_molecule,
		                               field_math::TErrorProbability(OMEGA), model, uniforms,
		                               tallies);

		const field_math::TLinkCounters kept     = merged(tallies);
		const field_math::TLinkCounters expected = recount(Y, Z_species, Z_molecule, pair);

		EXPECT_EQ(kept.total(), pair.species.n_leaves() * pair.molecule.n_leaves());
		for (size_t bucket = 0; bucket < field_math::TLinkCounters::n_buckets; ++bucket) {
			for (const bool y : {false, true}) {
				SCOPED_TRACE("bucket " + std::to_string(bucket) + ", field state " +
				             std::to_string(static_cast<int>(y)));
				EXPECT_EQ(kept.count(bucket, y), expected.count(bucket, y));
			}
		}
	}
}

/// The retally a held field gets counts the configuration it was handed, and moves no cell.
///
/// It is the same tally the drawing pass leaves, so it is asserted against the same naive recount
/// -- and against the pass itself, which is what makes the two agree for a run that turns `--fix_Y`
/// on halfway through its thinking.
TYPED_TEST(FieldUpdate, tally_counts_a_configuration_without_drawing) {
	using Field     = typename TestFixture::Field;
	using NodeState = typename TestFixture::NodeState;

	for (const auto &pair : tree_pairs()) {
		SCOPED_TRACE(pair.name);
		auto Y          = make_storage<Field>(field_shape(pair));
		auto Z_species  = make_storage<NodeState>(species_shape(pair));
		auto Z_molecule = make_storage<NodeState>(molecule_shape(pair));
		seed_ones(Y, 1);
		seed_ones(Z_species, 2);
		seed_ones(Z_molecule, 3);

		const std::vector<uint8_t> field_before    = states_of(Y);
		const std::vector<uint8_t> species_before  = states_of(Z_species);
		const std::vector<uint8_t> molecule_before = states_of(Z_molecule);

		const field_math::TLinkCounters kept =
		    field_update::tally<TLinkPolicy>(Y, Z_species, Z_molecule);
		const field_math::TLinkCounters expected = recount(Y, Z_species, Z_molecule, pair);

		EXPECT_EQ(kept.total(), pair.species.n_leaves() * pair.molecule.n_leaves());
		for (size_t bucket = 0; bucket < field_math::TLinkCounters::n_buckets; ++bucket) {
			for (const bool y : {false, true}) {
				SCOPED_TRACE("bucket " + std::to_string(bucket) + ", field state " +
				             std::to_string(static_cast<int>(y)));
				EXPECT_EQ(kept.count(bucket, y), expected.count(bucket, y));
			}
		}

		// Nothing moved. A held field is held, and the two tree fields are their trees' to draw.
		EXPECT_EQ(states_of(Y), field_before);
		EXPECT_EQ(states_of(Z_species), species_before);
		EXPECT_EQ(states_of(Z_molecule), molecule_before);

		// The retally splits the cells over the threads, so the six numbers must not move with the
		// thread count. Merging shares is exact integer addition, which is what makes that true.
		for (const size_t n_threads : {size_t{1}, size_t{4}}) {
			const TThreadCount threads(n_threads);
			const field_math::TLinkCounters got =
			    field_update::tally<TLinkPolicy>(Y, Z_species, Z_molecule);
			for (size_t bucket = 0; bucket < field_math::TLinkCounters::n_buckets; ++bucket) {
				for (const bool y : {false, true}) {
					EXPECT_EQ(got.count(bucket, y), expected.count(bucket, y))
					    << n_threads << " threads, bucket " << bucket << ", field state " << y;
				}
			}
		}
	}
}

/// The tally the drawing pass keeps and the tally the retally builds are the same six numbers over
/// one configuration. So a chain that holds its field is scored against what a chain that draws it
/// would have been scored against.
TYPED_TEST(FieldUpdate, tally_agrees_with_the_tally_the_pass_keeps) {
	using Field     = typename TestFixture::Field;
	using NodeState = typename TestFixture::NodeState;

	for (const auto &pair : tree_pairs()) {
		SCOPED_TRACE(pair.name);
		auto Y          = make_storage<Field>(field_shape(pair));
		auto Z_species  = make_storage<NodeState>(species_shape(pair));
		auto Z_molecule = make_storage<NodeState>(molecule_shape(pair));
		seed_ones(Y, 1);
		seed_ones(Z_species, 2);
		seed_ones(Z_molecule, 3);

		TFreeModel model;
		std::vector<field_math::TLinkCounters> tallies(ProgramOptions::NUMBER_OF_THREADS);
		const TCellUniforms uniforms(4242, TCellStream::field, 17);
		field_update::run<TLinkPolicy>(Y, Z_species, Z_molecule,
		                               field_math::TErrorProbability(OMEGA), model, uniforms,
		                               tallies);

		const field_math::TLinkCounters kept = merged(tallies);
		const field_math::TLinkCounters retallied =
		    field_update::tally<TLinkPolicy>(Y, Z_species, Z_molecule);
		for (size_t bucket = 0; bucket < field_math::TLinkCounters::n_buckets; ++bucket) {
			for (const bool y : {false, true}) {
				SCOPED_TRACE("bucket " + std::to_string(bucket) + ", field state " +
				             std::to_string(static_cast<int>(y)));
				EXPECT_EQ(kept.count(bucket, y), retallied.count(bucket, y));
			}
		}
	}
}

/// A tally short of one per thread is rejected before the region opens.
///
/// A thread writes the tally at its own index, so a short vector would be a write past the end from
/// inside a parallel region. The caller sizes the vector, because the caller merges it.
TYPED_TEST(FieldUpdate, rejects_fewer_tallies_than_threads) {
	using Field     = typename TestFixture::Field;
	using NodeState = typename TestFixture::NodeState;

	const auto &pair = tree_pairs().front();
	const TThreadCount threads(3);
	auto Y          = make_storage<Field>(field_shape(pair));
	auto Z_species  = make_storage<NodeState>(species_shape(pair));
	auto Z_molecule = make_storage<NodeState>(molecule_shape(pair));

	TFreeModel model;
	std::vector<field_math::TLinkCounters> tallies(ProgramOptions::NUMBER_OF_THREADS - 1);
	const TCellUniforms uniforms(4242, TCellStream::field, 0);
	EXPECT_THROW(field_update::run<TLinkPolicy>(Y, Z_species, Z_molecule,
	                                            field_math::TErrorProbability(OMEGA), model,
	                                            uniforms, tallies),
	             coretools::err::TError);
}

/// The two backends leave the same field and the same six counters.
///
/// The whole-binary gate (`just parity`) asserts this of a chain. Here it is asserted of one pass,
/// where a failure names the loop rather than the run that diverged from it. The sparse storages
/// hold only the cells they were given, so this is also where a deferred insert is compared against
/// the write the dense storages take in place.
TEST(FieldUpdate, gives_the_same_chain_under_both_backends) {
	for (const auto &pair : tree_pairs()) {
		SCOPED_TRACE(pair.name);

		const auto run_once = [&pair]<typename Field, typename NodeState>() {
			auto Y          = make_storage<Field>(field_shape(pair));
			auto Z_species  = make_storage<NodeState>(species_shape(pair));
			auto Z_molecule = make_storage<NodeState>(molecule_shape(pair));
			seed_ones(Y, 1);
			seed_ones(Z_species, 2);
			seed_ones(Z_molecule, 3);

			TFreeModel model;
			std::vector<field_math::TLinkCounters> tallies(ProgramOptions::NUMBER_OF_THREADS);
			const TCellUniforms uniforms(4242, TCellStream::field, 13);
			field_update::run<TLinkPolicy>(Y, Z_species, Z_molecule,
			                               field_math::TErrorProbability(OMEGA), model, uniforms,
			                               tallies);
			return std::tuple{states_of(Y), merged(tallies)};
		};

		const auto [dense_Y, dense_counters] =
		    run_once.template operator()<TStorageYDense, TStorageZDense>();
		const auto [sparse_Y, sparse_counters] =
		    run_once.template operator()<TStorageYSparse, TStorageZSparse>();

		EXPECT_EQ(dense_Y, sparse_Y);
		for (size_t bucket = 0; bucket < field_math::TLinkCounters::n_buckets; ++bucket) {
			for (const bool y : {false, true}) {
				EXPECT_EQ(dense_counters.count(bucket, y), sparse_counters.count(bucket, y));
			}
		}
	}
}

} // namespace
