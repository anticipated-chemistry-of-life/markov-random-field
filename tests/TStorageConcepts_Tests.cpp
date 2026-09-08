
#include "storages/storage_backend.h"
#include <cstddef>
#include <vector>

// The storage concepts are checked entirely at compile time: TStorageYMatrix.h and
// TStorageZMatrix.h assert that the sparse pair conforms, and storage_backend.h asserts it of the
// pair its two aliases select. What those cannot say is that the concepts *reject* anything
// -- a concept whose requires-expression named no member at all would satisfy them just as well.
//
// So this file is the other half: types that deliberately fall short, asserted not to conform.
// There is nothing to run, which is why there is no TEST here.

// The stub bodies below are never executed -- only their signatures are looked at -- so
// clang-tidy's "could be static" is answering a question this file is not asking.
// NOLINTBEGIN(readability-convert-member-functions-to-static)
namespace {

/// Everything StorageWindow asks for.
struct Window {
	[[nodiscard]] size_t size() const { return 0; }
	[[nodiscard]] bool is_one(size_t) const { return false; }
	[[nodiscard]] size_t linear_index(size_t) const { return 0; }
	void set_state(size_t, bool) {}
	std::vector<size_t> take_buffered_inserts() { return {}; }
	void close() {}
};
static_assert(StorageWindow<Window>, "the full window surface must conform");

/// The same window, with no way to close it. A window that never closes never commits the writes
/// it buffered, so the sparse storage would lose every insert.
struct UnclosableWindow {
	[[nodiscard]] size_t size() const { return 0; }
	[[nodiscard]] bool is_one(size_t) const { return false; }
	[[nodiscard]] size_t linear_index(size_t) const { return 0; }
	void set_state(size_t, bool) {}
	std::vector<size_t> take_buffered_inserts() { return {}; }
};
static_assert(!StorageWindow<UnclosableWindow>, "a window that cannot close must not conform");

/// Everything the shared concept asks for except `remove_zeros`.
struct AlmostAStorage {
	[[nodiscard]] bool is_one(size_t) const { return false; }
	void set_state(size_t, bool) {}
	void insert_one(size_t) {}
	void insert_zero(size_t) {}
	[[nodiscard]] size_t total_size_of_container_space() const { return 0; }
	[[nodiscard]] bool empty() const { return true; }
	[[nodiscard]] size_t get_linear_index_in_container_space(const IndexArray &) const { return 0; }
	[[nodiscard]] IndexArray get_multi_dimensional_index(size_t) const { return {}; }
};
static_assert(!BinaryStorage<AlmostAStorage>, "a missing member must not conform");

struct FullStorage : AlmostAStorage {
	void remove_zeros() {}
};
static_assert(BinaryStorage<FullStorage>, "the full shared surface must conform");
static_assert(!FieldStorage<FullStorage>, "no counter means it is not a field");
static_assert(!WindowedStorage<FullStorage>, "a storage with no window is not a windowed one");

/// The same storage, with the window the update loops reach the field and the node state through.
struct FullWindowedStorage : FullStorage {
	using TWindow = Window;
	TWindow open_window(const IndexArray &, size_t, size_t) { return {}; }
};
static_assert(WindowedStorage<FullWindowedStorage>, "the full windowed surface must conform");

/// `is_one` returning something convertible to bool is not the same as returning bool: a storage
/// answering with a count would read as "is one" for every stored cell.
struct WrongReturnType : FullStorage {
	[[nodiscard]] int is_one(size_t) const { return 0; }
};
static_assert(!BinaryStorage<WrongReturnType>, "the return types are part of the interface");

/// The window a storage hands out is part of what makes it a *windowed* storage: the update loops
/// read and write through it, so a storage whose window falls short is one they cannot use. It is
/// no part of being a storage, because a data source opens no window at all.
struct StorageWithAWindowThatFallsShort : FullStorage {
	using TWindow = UnclosableWindow;
	TWindow open_window(const IndexArray &, size_t, size_t) { return {}; }
};
static_assert(BinaryStorage<StorageWithAWindowThatFallsShort>,
              "a window is not what makes a storage a storage");
static_assert(!WindowedStorage<StorageWithAWindowThatFallsShort>,
              "the window is part of the windowed storage interface");

} // namespace
// NOLINTEND(readability-convert-member-functions-to-static)

// The refinement is a real one in both directions: the node state carries no posterior counter
// and must not pass for a field, while the field satisfies both.
static_assert(BinaryStorage<TNodeStateStorage>);
static_assert(!FieldStorage<TNodeStateStorage>);
static_assert(BinaryStorage<TFieldStorage>);
static_assert(FieldStorage<TFieldStorage>);

// The observations are storages and no more. Neither opens a window, and neither carries a counter.
static_assert(BinaryStorage<TBinaryStorage>);
static_assert(!FieldStorage<TBinaryStorage>);
static_assert(BinaryStorage<TSparseBinaryArray>);
static_assert(!WindowedStorage<TSparseBinaryArray>);
static_assert(!FieldStorage<TSparseBinaryArray>);
