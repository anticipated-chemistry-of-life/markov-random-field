
#include "storages/storage_backend.h"
#include <cstddef>
#include <cstdint>

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

/// A storage that answers `locate` as well: the state, whether it holds the cell, the linear index
/// and where the cell is.
struct FullLocatableStorage : FullStorage {
	using TCell = uint8_t;
	IsOneResult<TCell> locate(size_t) { return {}; }
	IsOneResult<TCell> locate(const IndexArray &) { return {}; }
};
static_assert(LocatableStorage<FullLocatableStorage>, "the full locatable surface must conform");
static_assert(!LocatableStorage<FullStorage>,
              "a storage that cannot point at a cell is not locatable");

/// The same storage, locating by linear index alone. An updater that holds a multidimensional
/// index would have to convert it first, and that arithmetic is what `locate` keeps in one place.
struct LocatableByLinearIndexAlone : FullStorage {
	using TCell = uint8_t;
	IsOneResult<TCell> locate(size_t) { return {}; }
};
static_assert(!LocatableStorage<LocatableByLinearIndexAlone>,
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

// The four storages that point an updater at one of their cells, and the two that do not. A
// sorted-vector matrix keeps every cell twice, once in its row and once in its column, so it has no
// single cell to point at.
static_assert(LocatableStorage<TDenseStateArray>);
static_assert(LocatableStorage<TSparseBinaryArray>);
static_assert(LocatableStorage<TStorageYDense>);
static_assert(LocatableStorage<TStorageZDense>);
static_assert(!LocatableStorage<TStorageYMatrix>);
static_assert(!LocatableStorage<TStorageZMatrix>);

// The observations are storages and no more. Neither carries a counter.
static_assert(BinaryStorage<TBinaryStorage>);
static_assert(!FieldStorage<TBinaryStorage>);
static_assert(BinaryStorage<TSparseBinaryArray>);
static_assert(!FieldStorage<TSparseBinaryArray>);
