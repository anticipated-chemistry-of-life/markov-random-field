//
// Runs a kernel over every clique of a tree, in parallel, one clique to a thread's share of work.
//
// TTree's own passes over its cliques -- the node-state walk, the alpha and nu moves, the
// branch-length likelihood, the density, the from-children start -- are five hand-written copies
// of the same shape: open a parallel region, build a clique view, run some per-clique work, drain
// the view's deferred inserts, and gather whatever the work returned in clique order. This is that
// shape, named once.
//
// It knows nothing of TTree. A view factory builds the clique view for one clique number, and a
// kernel does that clique's own work against it -- both supplied by the caller, so this file names
// no tree, no clique space and no parameter. That is what makes it constructible in a test from a
// bare phylogeny and node state, the way node_state_walk.h and node_state_density.h already are.
//
// The branch-length likelihood is the one pass that does not fit: its result is reduced across
// fixed slots rather than gathered one entry per clique, on a schedule this file does not use. It
// stays on its own hand-written loop.
//

#pragma once

#include "cli.h"

#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

namespace clique_traversal {

/// What a kernel needs to do one clique's share of the work: take the clique's number and the view
/// its caller built for it, and return whatever this pass gathers -- or nothing, for a pass that
/// only writes through the view. `View` is the view type the caller's view factory hands back, so
/// this is checked against the same type `run` deduces it from.
template<typename Kernel, typename View>
concept CliqueKernel = requires(Kernel &kernel, size_t clique, View &view) {
	kernel(clique, view);
};

/// Runs `kernel` over every one of `n_cliques` cliques, in parallel.
///
/// `view_factory(clique)` builds the view `kernel(clique, view)` is then handed -- the same
/// `TNodeStateCliqueView` a tree's own passes already open, or a test's own view over a bare
/// phylogeny and node state. Every view's deferred inserts are collected into one list per clique,
/// returned to the caller: a mutating caller feeds them to the node state's own insert; a read-only
/// caller asserts the whole result is empty once, rather than per clique inside the loop.
///
/// `kernel`'s return type decides what else comes back. A value-returning kernel has its results
/// gathered into a vector sized to `n_cliques`, in clique order -- the caller reads `.first`
/// alongside the inserts in `.second`. A `void` kernel gathers nothing, and `run` returns only the
/// inserts: the from-children start writes through its view and returns no value, so nothing here
/// pays for a gather vector it would never read.
///
/// `schedule(dynamic)` and `ProgramOptions::NUMBER_OF_THREADS`, `default(none)`: the schedule every
/// clique-parallel pass but the branch-length likelihood already agrees on. A tree's cliques can
/// cost wildly different amounts -- a clique's likelihood walk is bounded by that tree's depth, not
/// by a constant -- so a static split would leave some threads idle while one carries a deep
/// subtree. Nothing gathered here is summed across cliques inside this loop, so no ordering
/// discipline is owed to the schedule; a pass that does sum across cliques (the branch-length
/// likelihood) owns that discipline itself, on its own fixed-slot schedule, outside this file.
template<typename ViewFactory, typename Kernel>
[[nodiscard]] auto run(size_t n_cliques, ViewFactory &view_factory, Kernel &kernel) {
	using View = std::invoke_result_t<ViewFactory &, size_t>;
	static_assert(CliqueKernel<Kernel, View>,
	              "The kernel must be callable as kernel(clique, view).");
	using Result = std::invoke_result_t<Kernel &, size_t, View &>;

	std::vector<std::vector<size_t>> inserts(n_cliques);

	if constexpr (std::is_void_v<Result>) {
#pragma omp parallel for num_threads(ProgramOptions::NUMBER_OF_THREADS) schedule(dynamic)          \
    default(none) shared(n_cliques, view_factory, kernel, inserts)
		for (size_t clique = 0; clique < n_cliques; ++clique) {
			auto view = view_factory(clique);
			kernel(clique, view);
			inserts[clique] = view.take_deferred_inserts();
		}
		return inserts;
	} else {
		std::vector<Result> gathered(n_cliques);
#pragma omp parallel for num_threads(ProgramOptions::NUMBER_OF_THREADS) schedule(dynamic)          \
    default(none) shared(n_cliques, view_factory, kernel, inserts, gathered)
		for (size_t clique = 0; clique < n_cliques; ++clique) {
			auto view        = view_factory(clique);
			gathered[clique] = kernel(clique, view);
			inserts[clique]  = view.take_deferred_inserts();
		}
		return std::pair<std::vector<Result>, std::vector<std::vector<size_t>>>{std::move(gathered),
		                                                                       std::move(inserts)};
	}
}

} // namespace clique_traversal
