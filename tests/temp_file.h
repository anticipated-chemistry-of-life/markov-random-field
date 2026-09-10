//
// A file that lives for the duration of one test.
//
// Extracted from TPhylogeny_Tests.cpp, where it started, so that the node-state file suite writes
// and reads its fixtures the same way. Both suites want the same two things from it: a path
// outside the source tree, and a file that is gone when the test ends.
//

#pragma once

#include <atomic>
#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>

/// A file that only exists for the duration of one test, and never inside the source tree: a
/// failing test must not leave a stray fixture behind for the next run to trip over.
///
/// The two-argument form writes the content up front, which is what a reader under test needs.
/// The one-argument form names a path and creates nothing, which is what a writer under test
/// needs; `content` then reads back what the writer put there.
class TTempFile {
private:
	std::filesystem::path _path;

	static std::filesystem::path _unique_path(const std::string &name) {
		// A per-process salt and a counter, because two test binaries running side by side would
		// otherwise write the same path and read each other's content.
		static const auto salt = std::to_string(std::random_device{}());
		static std::atomic<unsigned> counter{0};
		return std::filesystem::temp_directory_path() /
		       ("acol_" + salt + "_" + std::to_string(counter++) + "_" + name);
	}

public:
	explicit TTempFile(const std::string &name) : _path(_unique_path(name)) {}

	TTempFile(const std::string &name, const std::string &content) : TTempFile(name) {
		std::ofstream out(_path);
		// Without this, a test that could not write its own fixture reads an absent file, the
		// reader throws, and the test passes for entirely the wrong reason. gtest assertions are
		// not usable in a constructor, so this throws instead.
		if (!out.is_open()) { throw std::runtime_error("could not create temp file " + path()); }
		out << content;
	}
	~TTempFile() {
		std::error_code ignored;
		std::filesystem::remove(_path, ignored);
	}
	TTempFile(const TTempFile &)            = delete;
	TTempFile &operator=(const TTempFile &) = delete;

	[[nodiscard]] std::string path() const { return _path.string(); }

	/// What the file holds now. Empty when nothing has written it.
	[[nodiscard]] std::string content() const {
		std::ifstream in(_path);
		std::ostringstream text;
		text << in.rdbuf();
		return text.str();
	}
};
