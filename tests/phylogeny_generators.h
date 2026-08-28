//
// The forests the property tests are run over.
//
// Extracted from TPhylogeny_Tests.cpp, where it started, so that the storage suite asserts its
// properties over the same shapes rather than over a fixture of its own. A second generator would
// be a second set of shapes to keep honest, and the shapes are the whole point: a balanced fixture
// satisfies an index property by accident, where a tree with a single leaf or a single internal
// node does not.
//

#pragma once

#include "tree/TPhylogeny.h"

#include <algorithm>
#include <cstddef>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace phylo {

inline TEdge edge(std::string child, std::string parent, double length = 1.0) {
	return TEdge{std::move(child), std::move(parent), length};
}

/// A random forest, as an edge list in arbitrary order. Every root is given a child up front, so
/// the "a root must have at least one child" post-condition is never the thing under test here.
inline std::vector<TEdge> random_forest(std::mt19937_64 &rng, size_t n_nodes, size_t n_roots) {
	std::vector<TEdge> edges;
	edges.reserve(n_nodes - n_roots); // one branch per non-root node
	auto name = [](size_t i) { return "n" + std::to_string(i); };

	// Nodes [0, n_roots) are roots; the next n_roots nodes are one guaranteed child each.
	for (size_t i = 0; i < n_roots; ++i) { edges.push_back(edge(name(n_roots + i), name(i))); }
	// Everything after that hangs off a uniformly chosen node that already exists.
	for (size_t i = 2 * n_roots; i < n_nodes; ++i) {
		std::uniform_int_distribution<size_t> pick(0, i - 1);
		edges.push_back(edge(name(i), name(pick(rng))));
	}
	std::shuffle(edges.begin(), edges.end(), rng);
	return edges;
}

/// A deep chain of `n_nodes` nodes, each the only child of the one before it: one root, one leaf,
/// and everything between them internal. The shape a balanced fixture never produces.
inline std::vector<TEdge> chain(size_t n_nodes) {
	std::vector<TEdge> edges;
	edges.reserve(n_nodes);
	for (size_t i = 1; i < n_nodes; ++i) {
		edges.push_back(edge("n" + std::to_string(i), "n" + std::to_string(i - 1)));
	}
	return edges;
}

/// A wide shallow star: one root with `n_leaves` children, and no node between them. The chain's
/// mirror image -- one internal node, and every other node a leaf.
inline std::vector<TEdge> star(size_t n_leaves) {
	std::vector<TEdge> edges;
	edges.reserve(n_leaves);
	for (size_t i = 1; i <= n_leaves; ++i) {
		edges.push_back(edge("n" + std::to_string(i), "root"));
	}
	return edges;
}

} // namespace phylo
