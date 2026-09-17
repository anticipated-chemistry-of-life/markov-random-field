//
// A test's own thread count, restored when the test ends.
//
// Extracted from TBlockUpdate_Tests.cpp, where it started, so that TCliqueTraversal_Tests.cpp
// runs the same guard rather than a second copy of it. Small, but the way
// tests/phylogeny_generators.h and tests/written_uniforms.h are already shared.
//

#pragma once

#include "cli.h"

#include <cstddef>

namespace threads {

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

} // namespace threads
