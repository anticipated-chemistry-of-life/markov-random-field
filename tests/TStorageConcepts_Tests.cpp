
#include "storages/TDense.h"
#include "storages/storage_backend.h"
#include "storages/y_storage/TStorageYSparse.h"
#include "storages/z_storage/TStorageZSparse.h"
#include <cstddef>
#include <cstdint>
#include <type_traits>

// The storage concepts are checked entirely at compile time: every storage header asserts that its
// own type conforms, and storage_backend.h asserts it of the pair its two aliases select. What
// those cannot say is that the concepts *reject* anything -- a concept whose requires-expression
// named no member at all would satisfy them just as well.
//
// So this file is the other half: types that deliberately fall short, asserted not to conform.
// There is nothing to run, which is why there is no TEST here.

// The stub bodies below are never executed -- only their signatures are looked at -- so
// clang-tidy's "could be static" is answering a question this file is not asking.
// NOLINTBEGIN(readability-convert-member-functions-to-static)
namespace {

/// Everything the shared concept asks for except `remove_zeros`.
struct AlmostAStorage {
	using TCell = uint8_t;
	[[nodiscard]] bool is_one(size_t) const { return false; }
	void set_state(size_t, bool) {}
	void insert_one(size_t) {}
	void insert_zero(size_t) {}
	[[nodiscard]] size_t total_size_of_container_space() const { return 0; }
	[[nodiscard]] bool empty() const { return true; }
	[[nodiscard]] size_t get_linear_index_in_container_space(const IndexArray &) const { return 0; }
	[[nodiscard]] IndexArray get_multi_dimensional_index(size_t) const { return {}; }
	IsOneResult<TCell> locate(size_t) { return {}; }
	IsOneResult<TCell> locate(const IndexArray &) { return {}; }
};
static_assert(!BinaryStorage<AlmostAStorage>, "a missing member must not conform");

struct FullStorage : AlmostAStorage {
	void remove_zeros() {}
};
static_assert(BinaryStorage<FullStorage>, "the full shared surface must conform");
static_assert(!FieldStorage<FullStorage>, "no counter means it is not a field");

/// A storage that cannot point an updater at a cell. Every storage has to: an update writes
/// through a handle, and there is no second way in.
struct WithoutAHandle {
	[[nodiscard]] bool is_one(size_t) const { return false; }
	void set_state(size_t, bool) {}
	void insert_one(size_t) {}
	void insert_zero(size_t) {}
	void remove_zeros() {}
	[[nodiscard]] size_t total_size_of_container_space() const { return 0; }
	[[nodiscard]] bool empty() const { return true; }
	[[nodiscard]] size_t get_linear_index_in_container_space(const IndexArray &) const { return 0; }
	[[nodiscard]] IndexArray get_multi_dimensional_index(size_t) const { return {}; }
};
static_assert(!BinaryStorage<WithoutAHandle>,
              "a storage that cannot point at a cell is not a storage");

/// The same storage, locating by linear index alone. An updater that holds a multidimensional
/// index would have to convert it first, and that arithmetic is what `locate` keeps in one place.
struct LocatableByLinearIndexAlone : WithoutAHandle {
	using TCell = uint8_t;
	IsOneResult<TCell> locate(size_t) { return {}; }
};
static_assert(!BinaryStorage<LocatableByLinearIndexAlone>,
              "both index forms are part of the interface");

/// `is_one` returning something convertible to bool is not the same as returning bool: a storage
/// answering with a count would read as "is one" for every stored cell.
struct WrongReturnType : FullStorage {
	[[nodiscard]] int is_one(size_t) const { return 0; }
};
static_assert(!BinaryStorage<WrongReturnType>, "the return types are part of the interface");

} // namespace
// NOLINTEND(readability-convert-member-functions-to-static)

// The refinement is a real one in both directions: the node state carries no posterior counter
// and must not pass for a field, while the field satisfies both.
static_assert(BinaryStorage<TNodeStateStorage>);
static_assert(!FieldStorage<TNodeStateStorage>);
static_assert(BinaryStorage<TFieldStorage>);
static_assert(FieldStorage<TFieldStorage>);

// Every storage points an updater at one of its cells, so `BinaryStorage` above already says it of
// all six. What is left to check is that the two cells they point at are the two that exist: a
// bare state, and the state packed with the posterior counter beside it.
static_assert(std::is_same_v<TDenseBinary::TCell, uint8_t>);
static_assert(std::is_same_v<TSparseBinary::TCell, uint8_t>);
static_assert(std::is_same_v<TStorageZDense::TCell, uint8_t>);
static_assert(std::is_same_v<TStorageZSparse::TCell, uint8_t>);
static_assert(std::is_same_v<TStorageYDense::TCell, TStorageY>);
static_assert(std::is_same_v<TStorageYSparse::TCell, TStorageY>);

// The observations are storages and no more. Neither carries a counter.
static_assert(BinaryStorage<TBinaryStorage>);
static_assert(!FieldStorage<TBinaryStorage>);
static_assert(BinaryStorage<TSparseBinary>);
static_assert(!FieldStorage<TSparseBinary>);
