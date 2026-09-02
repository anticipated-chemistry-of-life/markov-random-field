//
// The four storage pairings, and the container shapes a pair of trees asks them to hold.
//
// A suite that runs against the storages rather than against one build's aliases needs the same
// three things every time: the shapes, a way to build a storage of either kind, and the list of
// pairings to instantiate over. They live here so that two suites assert over one set of shapes.
// A second set would be a second thing to keep honest, and the shapes are the point: a chain has
// a single leaf, a star has a single internal node, and a container one cell wide is where an
// index property stops being a coincidence.
//
// Continuous integration gates two of the four pairings (`just parity`), so a suite instantiated
// over this list is where the other two are exercised at all.
//

#pragma once

#include "constants.h"
#include "phylogeny_generators.h"
#include "storages/storage_concepts.h"
#include "storages/y_storage/TStorageYDense.h"
#include "storages/y_storage/TStorageYMatrix.h"
#include "storages/z_storage/TStorageZDense.h"
#include "storages/z_storage/TStorageZMatrix.h"
#include "tree/TPhylogeny.h"
#include "tree/node_state_shape.h"
#include "gtest/gtest.h"

#include <cstddef>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

namespace backends {

/// The chain length a field is sized with. Long enough that a thinned counter is exercised, and
/// short enough that a suite runs in milliseconds.
constexpr size_t N_ITERATIONS = 300;

/// One pairing of a species tree with a molecule tree. Both are real phylogenies rather than a
/// fixture, because the shape is what decides whether a property is a property.
struct TTreePair {
	std::string name;
	TPhylogeny species;
	TPhylogeny molecule;
};

/// Built once: a phylogeny is immutable, and every suite asserts over the same three pairings.
inline const std::vector<TTreePair> &tree_pairs() {
	static const std::vector<TTreePair> pairs = [] {
		std::mt19937_64 rng(20260902);
		std::vector<TTreePair> built;
		built.push_back({"forest_x_forest", build_phylogeny(phylo::random_forest(rng, 18, 2)),
		                 build_phylogeny(phylo::random_forest(rng, 14, 1))});
		// a chain has a single leaf, so it makes a container one cell wide in that dimension
		built.push_back(
		    {"star_x_chain", build_phylogeny(phylo::star(9)), build_phylogeny(phylo::chain(6))});
		built.push_back(
		    {"chain_x_star", build_phylogeny(phylo::chain(6)), build_phylogeny(phylo::star(9))});
		return built;
	}();
	return pairs;
}

/// The three container spaces a pair of trees implies. The two node-state shapes come from
/// `node_state_dimensions`, the same function TTree sizes the real node state with, so changing
/// that rule changes what these suites are asserted over.
inline IndexArray field_shape(const TTreePair &pair) {
	return IndexArray{pair.species.n_leaves(), pair.molecule.n_leaves()};
}
inline IndexArray species_shape(const TTreePair &pair) {
	return node_state_dimensions(field_shape(pair), 0, pair.species);
}
inline IndexArray molecule_shape(const TTreePair &pair) {
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
/// so a cell left out is a cell a later write has to defer.
template<typename Storage> void seed_ones(Storage &storage, uint64_t salt) {
	const size_t total = storage.total_size_of_container_space();
	for (size_t i = 0; i < total; ++i) {
		if (((i * 2654435761ULL + salt * 40503ULL) >> 5U) % 3U == 0U) { storage.insert_one(i); }
	}
}

/// One field storage against one node-state storage, as a typed test takes them.
template<typename Field, typename NodeState> struct TBackends {
	using field      = Field;
	using node_state = NodeState;
};

/// Every pairing of the two field storages with the two node-state storages.
using AllBackends = ::testing::Types<
    TBackends<TStorageYDense, TStorageZDense>, TBackends<TStorageYDense, TStorageZMatrix>,
    TBackends<TStorageYMatrix, TStorageZDense>, TBackends<TStorageYMatrix, TStorageZMatrix>>;

} // namespace backends
