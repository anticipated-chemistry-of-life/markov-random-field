
#include "storages/storage_backend.h"
#include <cstddef>
#include <cstdint>
#include <vector>

// The storage concepts are checked entirely at compile time: TStorageYMatrix.h and
// TStorageZMatrix.h assert that the sparse pair conforms, and storage_backend.h asserts it of
// whichever pair the build selected. What those cannot say is that the concepts *reject* anything
// -- a concept whose requires-expression named no member at all would satisfy them just as well.
//
// So this file is the other half: types that deliberately fall short, asserted not to conform.
// There is nothing to run, which is why there is no TEST here.

// The stub bodies below are never executed -- only their signatures are looked at -- so
// clang-tidy's "could be static" is answering a question this file is not asking.
// NOLINTBEGIN(readability-convert-member-functions-to-static)
namespace {

/// Everything the shared concept asks for except `remove_zeros`.
struct AlmostBinaryStorage {
	[[nodiscard]] bool is_one(size_t) const { return false; }
	void set_state(size_t, bool) {}
	void insert_one(size_t) {}
	void insert_zero(size_t) {}
	[[nodiscard]] size_t total_size_of_container_space() const { return 0; }
	[[nodiscard]] bool empty() const { return true; }
	[[nodiscard]] size_t get_linear_index_in_container_space(const IndexArray &) const { return 0; }
	[[nodiscard]] IndexArray get_multi_dimensional_index(size_t) const { return {}; }
	void fill_current_state(const IndexArray &, size_t, size_t, std::vector<uint8_t> &,
	                        std::vector<uint8_t> &, std::vector<size_t> &) const {}
};
static_assert(!BinaryFieldStorage<AlmostBinaryStorage>, "a missing member must not conform");

struct BinaryStorage : AlmostBinaryStorage {
	void remove_zeros() {}
};
static_assert(BinaryFieldStorage<BinaryStorage>, "the full shared surface must conform");
static_assert(!FieldStorage<BinaryStorage>, "no counter means it is not a field");

/// `is_one` returning something convertible to bool is not the same as returning bool: a storage
/// answering with a count would read as "is one" for every stored cell.
struct WrongReturnType : BinaryStorage {
	[[nodiscard]] int is_one(size_t) const { return 0; }
};
static_assert(!BinaryFieldStorage<WrongReturnType>, "the return types are part of the interface");

} // namespace
// NOLINTEND(readability-convert-member-functions-to-static)

// The refinement is a real one in both directions: the internal state carries no posterior counter
// and must not pass for a field, while the field satisfies both.
static_assert(BinaryFieldStorage<TInternalStateStorage>);
static_assert(!FieldStorage<TInternalStateStorage>);
static_assert(BinaryFieldStorage<TFieldStorage>);
static_assert(FieldStorage<TFieldStorage>);
