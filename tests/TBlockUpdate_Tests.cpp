//
// The block update's loop, over every pairing of the two storages.
//
// The kernel is pure and tested against brute force (TFieldMath_Tests.cpp), and the windows are
// conformance-tested against each other (TStorageConformance_Tests.cpp). What is left belongs to
// the loop alone, and that is what this file asserts: every leaf pair visited exactly once, the two
// tree parents read being the right cells, the two data terms landing on the leaf pair they were
// scored for, the writes reaching the storage, and the same chain whatever the thread count.
//
// The loop is asked those questions through a model of its own rather than through the trees and
// the data sources, because a model the test writes can be *driven*: giving one tree field state no
// mass at all, and the other field state a data likelihood of zero, leaves the draw exactly one
// state to take. Every write is then a known value at a known cell, and a write that lands on the
// wrong leaf pair is a wrong value rather than a coincidence.
//
// Every body is instantiated over all four field/node-state pairings. Continuous integration gates
// two of them (`just parity`), so this is where the other two are exercised at all.
//

#include "constants.h"
#include "cli.h"
#include "coretools/Types/probability.h"
#include "field/TBlockUpdate.h"
#include "field/TFieldMath.h"
#include "field/link_backend.h"
#include "phylogeny_generators.h"
#include "random/TCellUniforms.h"
#include "storages/storage_concepts.h"
#include "storages/y_storage/TStorageYDense.h"
#include "storages/y_storage/TStorageYMatrix.h"
#include "storages/z_storage/TStorageZDense.h"
#include "storages/z_storage/TStorageZMatrix.h"
#include "tree/TPhylogeny.h"
#include "tree/node_state_shape.h"
#include "gtest/gtest.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

namespace {

constexpr size_t N_ITERATIONS = 300;
/// Small enough that the block never sees a degenerate link, and small enough to be recognisable
/// in a failure message.
constexpr double OMEGA = 0.125;

// -------------------------------------------------------------------------
// The shapes the loop is run over
// -------------------------------------------------------------------------

/// One pairing of a species tree with a molecule tree. Both trees are real phylogenies rather than
/// a fixture, because the shape is what decides whether a property is a property: a chain has one
/// leaf, and a star has one internal node.
struct TTreePair {
	std::string name;
	TPhylogeny species;
	TPhylogeny molecule;
};

const std::vector<TTreePair> &tree_pairs() {
	static const std::vector<TTreePair> pairs = [] {
		std::mt19937_64 rng(20260902);
		std::vector<TTreePair> built;
		built.push_back({"forest_x_forest", build_phylogeny(phylo::random_forest(rng, 18, 2)),
		                 build_phylogeny(phylo::random_forest(rng, 14, 1))});
		// a chain has a single leaf, so it makes a container one cell wide in that dimension
		built.push_back({"star_x_chain", build_phylogeny(phylo::star(9)),
		                 build_phylogeny(phylo::chain(6))});
		built.push_back({"chain_x_star", build_phylogeny(phylo::chain(6)),
		                 build_phylogeny(phylo::star(9))});
		return built;
	}();
	return pairs;
}

/// The three container spaces a pair of trees implies. Taken from `node_state_dimensions`, the
/// same function TTree sizes the real node state with.
IndexArray field_shape(const TTreePair &pair) {
	return IndexArray{pair.species.n_leaves(), pair.molecule.n_leaves()};
}
IndexArray species_shape(const TTreePair &pair) {
	return node_state_dimensions(field_shape(pair), 0, pair.species);
}
IndexArray molecule_shape(const TTreePair &pair) {
	return node_state_dimensions(field_shape(pair), 1, pair.molecule);
}

/// Both storages take their dimensions at construction; a field takes the chain length too, because
/// that is what its counter is sized from.
template<typename Storage> Storage make_storage(const IndexArray &dimensions) {
	if constexpr (FieldStorage<Storage>) {
		return Storage(N_ITERATIONS, dimensions);
	} else {
		return Storage(dimensions);
	}
}

/// A scattering of ones over a container space. The sparse backends hold only what is put in them,
/// so a cell left out is a cell a later write has to defer -- which is the path this seeds.
template<typename Storage> void seed_ones(Storage &storage, uint64_t salt) {
	const size_t total = storage.total_size_of_container_space();
	for (size_t i = 0; i < total; ++i) {
		if (((i * 2654435761ULL + salt * 40503ULL) >> 5U) % 3U == 0U) { storage.insert_one(i); }
	}
}

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

/// The states one leaf pair is driven to. An arbitrary but fixed function of the leaf pair, so that
/// a write landing on the wrong cell writes the wrong value.
field_math::TBlockStates target_at(size_t species_leaf, size_t molecule_leaf) {
	const uint64_t bits = (6364136223846793005ULL * (species_leaf + 1)) ^
	                      (1442695040888963407ULL * (molecule_leaf + 3));
	return {.y   = ((bits >> 17U) & 1U) != 0U,
	        .z_s = ((bits >> 29U) & 1U) != 0U,
	        .z_m = ((bits >> 41U) & 1U) != 0U};
}

/// A model that drives every leaf pair to `target_at`, and records what the loop asked it.
///
/// The forcing is what makes the recording checkable. A tree factor of 0 leaves that tree field
/// state with no mass, and a data likelihood of 0 for one field state does the same for the field,
/// so the eight-state draw has exactly one state left to take whatever uniform it is given.
class TForcingModel {
public:
	/// What the loop did at one leaf pair.
	struct TVisit {
		size_t n_asked   = 0; ///< how often `factors` was asked about this leaf pair
		size_t n_recorded = 0; ///< how often `record` was told about it
		bool species_parent  = false;
		bool molecule_parent = false;
		field_math::TBlockStates drawn;
	};

	class TRow {
	private:
		TForcingModel *_model;
		size_t _species_leaf;

	public:
		TRow(TForcingModel &model, size_t species_leaf)
		    : _model(&model), _species_leaf(species_leaf) {}

		[[nodiscard]] block_update::TLeafPairFactors factors(size_t molecule_leaf, bool species_parent,
		                                                 bool molecule_parent) const {
			TVisit &visit         = _model->visit(_species_leaf, molecule_leaf);
			++visit.n_asked;
			visit.species_parent  = species_parent;
			visit.molecule_parent = molecule_parent;

			const auto target = target_at(_species_leaf, molecule_leaf);
			return {.prob_z_s_is_one = coretools::P(target.z_s ? 1.0 : 0.0),
			        .prob_z_m_is_one = coretools::P(target.z_m ? 1.0 : 0.0),
			        .lotus = {coretools::P(target.y ? 0.0 : 1.0), coretools::P(target.y ? 1.0 : 0.0)},
			        .simple_error = {coretools::P(1.0), coretools::P(1.0)}};
		}

		void record(size_t molecule_leaf, const block_update::TLeafPairFactors &,
		            const field_math::TBlockStates &drawn) const {
			TVisit &visit = _model->visit(_species_leaf, molecule_leaf);
			++visit.n_recorded;
			visit.drawn = drawn;
		}
	};

	TForcingModel(size_t n_species_leaves, size_t n_molecule_leaves)
	    : _n_molecule_leaves(n_molecule_leaves), _visits(n_species_leaves * n_molecule_leaves) {}

	[[nodiscard]] TRow open_row(size_t species_leaf) { return TRow(*this, species_leaf); }

	[[nodiscard]] TVisit &visit(size_t species_leaf, size_t molecule_leaf) {
		return _visits[species_leaf * _n_molecule_leaves + molecule_leaf];
	}

private:
	size_t _n_molecule_leaves;
	std::vector<TVisit> _visits;
};

static_assert(block_update::BlockModel<TForcingModel>,
              "The forcing model must answer what the block update asks a model.");

/// A model that leaves the draw a real choice at every leaf pair, and keeps no state of its own.
///
/// The forcing model above pins every cell, which would let a wrong uniform pass unnoticed. Here
/// all eight states carry mass, so the state a leaf pair ends in depends on the uniform it drew --
/// which is what makes a chain comparable between two thread counts. Being stateless is what makes
/// it safe to run on many threads.
class TFreeModel {
public:
	class TRow {
	private:
		size_t _species_leaf;

	public:
		explicit TRow(size_t species_leaf) : _species_leaf(species_leaf) {}

		[[nodiscard]] block_update::TLeafPairFactors factors(size_t molecule_leaf, bool species_parent,
		                                                 bool molecule_parent) const {
			// Values that depend on both the leaf pair and the parents, so that a read of the
			// wrong parent, or of the wrong cell, moves the chain.
			const double drift = 0.1 * static_cast<double>((_species_leaf + molecule_leaf) % 4U);
			return {.prob_z_s_is_one = coretools::P(species_parent ? 0.7 : 0.2 + drift),
			        .prob_z_m_is_one = coretools::P(molecule_parent ? 0.65 : 0.15 + drift),
			        .lotus           = {coretools::P(0.4), coretools::P(0.6 - drift)},
			        .simple_error    = {coretools::P(0.55), coretools::P(0.45)}};
		}

		void record(size_t, const block_update::TLeafPairFactors &,
		            const field_math::TBlockStates &) const {}
	};

	[[nodiscard]] static TRow open_row(size_t species_leaf) { return TRow(species_leaf); }
};

static_assert(block_update::BlockModel<TFreeModel>,
              "The free model must answer what the block update asks a model.");

/// The six counters recomputed from a whole configuration, sharing nothing with the tally the
/// update kept as it went.
template<typename Field, typename NodeState>
field_math::TLinkCounters recount(const Field &Y, const NodeState &Z_species,
                                  const NodeState &Z_molecule, const TTreePair &pair) {
	field_math::TLinkCounters counters;
	for (size_t s = 0; s < pair.species.n_leaves(); ++s) {
		for (size_t m = 0; m < pair.molecule.n_leaves(); ++m) {
			const bool z_s = Z_species.is_one(Z_species.get_linear_index_in_container_space({s, m}));
			const bool z_m =
			    Z_molecule.is_one(Z_molecule.get_linear_index_in_container_space({s, m}));
			counters.add(TLinkPolicy::bucket(z_s, z_m),
			             Y.is_one(Y.get_linear_index_in_container_space({s, m})));
		}
	}
	return counters;
}

/// The tallies of one run, merged the way the caller of a block update merges them.
field_math::TLinkCounters merged(const std::vector<block_update::TThreadTally> &tallies) {
	field_math::TLinkCounters counters;
	for (const auto &tally : tallies) { counters.merge(tally.counters); }
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

template<typename Field, typename NodeState> struct TBackends {
	using field      = Field;
	using node_state = NodeState;
};

using AllBackends = ::testing::Types<TBackends<TStorageYDense, TStorageZDense>,
                                     TBackends<TStorageYDense, TStorageZMatrix>,
                                     TBackends<TStorageYMatrix, TStorageZDense>,
                                     TBackends<TStorageYMatrix, TStorageZMatrix>>;

template<typename Backends> class BlockUpdate : public ::testing::Test {
public:
	using Field     = typename Backends::field;
	using NodeState = typename Backends::node_state;
};

TYPED_TEST_SUITE(BlockUpdate, AllBackends);

/// Every leaf pair is asked about once and told what it was given once, and no leaf pair is missed.
TYPED_TEST(BlockUpdate, visits_every_leaf_pair_exactly_once) {
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
		std::vector<block_update::TThreadTally> tallies(ProgramOptions::NUMBER_OF_THREADS);
		const TCellUniforms uniforms(4242, TCellStream::field, 0);
		block_update::run<TLinkPolicy>(Y, Z_species, Z_molecule, pair.species, pair.molecule,
		                               field_math::TErrorProbability(OMEGA), model, uniforms,
		                               tallies);

		for (size_t s = 0; s < pair.species.n_leaves(); ++s) {
			for (size_t m = 0; m < pair.molecule.n_leaves(); ++m) {
				SCOPED_TRACE("leaf pair " + std::to_string(s) + "," + std::to_string(m));
				EXPECT_EQ(model.visit(s, m).n_asked, 1U);
				EXPECT_EQ(model.visit(s, m).n_recorded, 1U);
			}
		}
	}
}

/// The three states the draw assigned are in the three containers afterwards, at the leaf pair's
/// own cell -- through an in-place write where the backend held the cell, and through the deferred
/// insert where it did not.
TYPED_TEST(BlockUpdate, writes_the_drawn_states_back) {
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
		std::vector<block_update::TThreadTally> tallies(ProgramOptions::NUMBER_OF_THREADS);
		const TCellUniforms uniforms(4242, TCellStream::field, 0);
		block_update::run<TLinkPolicy>(Y, Z_species, Z_molecule, pair.species, pair.molecule,
		                               field_math::TErrorProbability(OMEGA), model, uniforms,
		                               tallies);

		for (size_t s = 0; s < pair.species.n_leaves(); ++s) {
			for (size_t m = 0; m < pair.molecule.n_leaves(); ++m) {
				SCOPED_TRACE("leaf pair " + std::to_string(s) + "," + std::to_string(m));
				const auto target = target_at(s, m);
				// the draw had one state to take, so the model was driven to the target
				EXPECT_EQ(model.visit(s, m).drawn.y, target.y);
				EXPECT_EQ(model.visit(s, m).drawn.z_s, target.z_s);
				EXPECT_EQ(model.visit(s, m).drawn.z_m, target.z_m);
				// and the target is what each container holds at that leaf pair
				EXPECT_EQ(Y.is_one(IndexArray{s, m}), target.y);
				EXPECT_EQ(Z_species.is_one(Z_species.get_linear_index_in_container_space({s, m})),
				          target.z_s);
				EXPECT_EQ(Z_molecule.is_one(Z_molecule.get_linear_index_in_container_space({s, m})),
				          target.z_m);
			}
		}
	}
}

/// The two tree parents the loop reads are the cells the two topologies name: the species parent's
/// cell in the same column, and the molecule parent's cell in the same row.
TYPED_TEST(BlockUpdate, reads_the_tree_parent_of_each_leaf_pair) {
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

		// What the parents read before the update. The update writes leaves only, and a parent is
		// an internal node, so these are what the loop had to see.
		std::vector<uint8_t> species_before  = states_of(Z_species);
		std::vector<uint8_t> molecule_before = states_of(Z_molecule);

		TForcingModel model(pair.species.n_leaves(), pair.molecule.n_leaves());
		std::vector<block_update::TThreadTally> tallies(ProgramOptions::NUMBER_OF_THREADS);
		const TCellUniforms uniforms(4242, TCellStream::field, 0);
		block_update::run<TLinkPolicy>(Y, Z_species, Z_molecule, pair.species, pair.molecule,
		                               field_math::TErrorProbability(OMEGA), model, uniforms,
		                               tallies);

		for (size_t s = 0; s < pair.species.n_leaves(); ++s) {
			for (size_t m = 0; m < pair.molecule.n_leaves(); ++m) {
				SCOPED_TRACE("leaf pair " + std::to_string(s) + "," + std::to_string(m));
				const size_t species_parent = Z_species.get_linear_index_in_container_space(
				    IndexArray{pair.species.parent_of(s), m});
				const size_t molecule_parent = Z_molecule.get_linear_index_in_container_space(
				    IndexArray{s, pair.molecule.parent_of(m)});
				EXPECT_EQ(model.visit(s, m).species_parent, species_before[species_parent] != 0);
				EXPECT_EQ(model.visit(s, m).molecule_parent, molecule_before[molecule_parent] != 0);
			}
		}
	}
}

/// One thread and many give the same three containers and the same six counters. A cell's uniform
/// is hashed from its position (ADR-0007), so the thread that reaches it does not decide what it
/// gets.
TYPED_TEST(BlockUpdate, gives_the_same_chain_at_any_thread_count) {
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
			std::vector<block_update::TThreadTally> tallies(n_threads);
			const TCellUniforms uniforms(4242, TCellStream::field, 7);
			block_update::run<TLinkPolicy>(Y, Z_species, Z_molecule, pair.species, pair.molecule,
			                               field_math::TErrorProbability(OMEGA), model, uniforms,
			                               tallies);
			return std::tuple{states_of(Y), states_of(Z_species), states_of(Z_molecule),
			                  merged(tallies)};
		};

		const auto [one_Y, one_species, one_molecule, one_counters]    = run_once(1);
		const auto [many_Y, many_species, many_molecule, many_counters] = run_once(4);

		EXPECT_EQ(one_Y, many_Y);
		EXPECT_EQ(one_species, many_species);
		EXPECT_EQ(one_molecule, many_molecule);
		for (size_t bucket = 0; bucket < field_math::TLinkCounters::n_buckets; ++bucket) {
			for (const bool y : {false, true}) {
				EXPECT_EQ(one_counters.count(bucket, y), many_counters.count(bucket, y));
			}
		}
	}
}

/// The counters the update accumulated are the tally of the configuration it left behind, and they
/// count every leaf pair once.
TYPED_TEST(BlockUpdate, counters_tally_the_configuration_it_left) {
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
		std::vector<block_update::TThreadTally> tallies(ProgramOptions::NUMBER_OF_THREADS);
		const TCellUniforms uniforms(4242, TCellStream::field, 11);
		block_update::run<TLinkPolicy>(Y, Z_species, Z_molecule, pair.species, pair.molecule,
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

} // namespace
