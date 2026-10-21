#pragma once
//
// The CPU profile of a whole run, in a build configured with ACOL_PROFILING (the `profile-*`
// presets, `pixi run infer profile-lm ...`). Anywhere else `TProgramProfile` does nothing.
//
// `main` holds one for the whole run, so the profile covers all of it: reading the data, the
// burn-in, the chain start, the chain, and writing the posteriors. It is written to
// <out>_profile.prof when the run names its output prefix on the command line (`--out X` or
// `out=X`), and to acol_profile.prof otherwise -- the profile starts before coretools reads its
// parameters, so a prefix given in a parameter file does not reach it.
//
// gperftools samples on SIGPROF, so it needs neither root nor perf_event_paranoid. It reads its
// sampling rate from CPUPROFILE_FREQUENCY when the binary starts, in samples per CPU-second
// (default 100); CPUPROFILE_REALTIME=1 samples wall-clock time instead. CPUPROFILE must stay unset:
// it starts a second profile of its own.
//
// Read the profile with pprof (github.com/google/pprof):
//     pprof -top    build/profile-lm/acol <out>_profile.prof
//     pprof -http=: build/profile-lm/acol <out>_profile.prof
//

#include <cstdio>
#include <string>
#include <string_view>

#ifdef ACOL_PROFILING
#include <gperftools/profiler.h>
#endif

namespace profiling {

/// Where the profile of a run with these arguments goes: `<out>_profile.prof`, reading `out` the
/// two ways coretools accepts it on a command line, or `acol_profile.prof` when it is not there.
[[nodiscard]] inline std::string profile_file(int argc, char *argv[]) {
	constexpr std::string_view spaced = "--out";
	constexpr std::string_view equals = "out=";
	for (int i = 1; i < argc; ++i) {
		const std::string_view arg = argv[i];
		if (arg == spaced && i + 1 < argc) { return std::string(argv[i + 1]) + "_profile.prof"; }
		if (arg.starts_with(equals)) { return std::string(arg.substr(equals.size())) + "_profile.prof"; }
	}
	return "acol_profile.prof";
}

/// Profiles the process from construction to destruction. `main` holds one around everything it
/// does, so the profile is written however the run ends, a thrown error included.
class TProgramProfile {
public:
	TProgramProfile([[maybe_unused]] int argc, [[maybe_unused]] char *argv[]) {
#ifdef ACOL_PROFILING
		const std::string file = profile_file(argc, argv);
		if (ProfilerStart(file.c_str()) == 0) {
			// Before coretools has a log to write to, and not a reason to stop the run.
			std::fprintf(stderr,
			             "acol: cannot profile into '%s': the file cannot be written, or a profile "
			             "is already running (is CPUPROFILE set?). Running without a profile.\n",
			             file.c_str());
			return;
		}
		std::fprintf(stderr, "acol: profiling the whole run into '%s'.\n", file.c_str());
#endif
	}

	~TProgramProfile() {
#ifdef ACOL_PROFILING
		ProfilerStop();
#endif
	}

	TProgramProfile(const TProgramProfile &)            = delete;
	TProgramProfile &operator=(const TProgramProfile &) = delete;
};

} // namespace profiling
