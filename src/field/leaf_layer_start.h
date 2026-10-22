//
// Where a chain's leaf layer starts.
//
// The block update draws the field and both tree fields together, so the three of them start
// together too. The start is two steps. The field takes the LOTUS records, and both tree fields
// take the field. CONTEXT.md names the concept under *Chain start*, and ADR-0005 argues it.
//
// Both steps read and write through the storage concepts and hold no tree. One suite therefore
// asserts them over every backend pairing, and not over the one pairing a chain is built with.
//

#pragma once

#include "constants.h"
#include "coretools/Main/TError.h"
#include "field/TFieldMath.h"
#include "storages/cell_write.h"
#include "storages/storage_concepts.h"
#include <algorithm>
#include <cstddef>
#include <vector>

namespace leaf_layer_start {

/// Puts the field at the records. It then reads one exactly where a record exists, and zero
/// elsewhere.
///
/// The field must be empty. This writes ones and never zeros, because a cell no record holds
/// already reads as zero, and inserting that zero would grow a sparse field for nothing. A field
/// that already held states would keep them, so the postcondition above rests on the field being
/// empty rather than on the caller.
///
/// The records and the field are addressed the same way -- both are indexed in leaf space, one
/// dimension per tree -- so a record's linear index is already the field's. They are two types all
/// the same: the records carry no posterior counter, and are a binary storage rather than a field.
template<BinaryStorage Records, FieldStorage Field>
void start_the_field_at(const Records &records, Field &field) {
	if (records.total_size_of_container_space() != field.total_size_of_container_space()) {
		throw coretools::TDevError("Cannot start the field at the records: the records hold ",
		                           records.total_size_of_container_space(),
		                           " cells and the field holds ",
		                           field.total_size_of_container_space(), ".");
	}
	if (!field.empty()) {
		throw coretools::TDevError(
		    "Cannot start the field at the records: the field already holds states.");
	}
	for (size_t index = 0; index < field.total_size_of_container_space(); ++index) {
		if (records.is_one(index)) { field.insert_one(index); }
	}
}

/// How many species leaves the copy below takes at a time. See `hold_tree_fields_at_the_field`.
inline constexpr size_t BAND_ROWS = 64;

/// Puts both tree fields at the field, and returns the six counters over the whole of it.
///
/// The tally is degenerate by construction. Both tree fields read the field, so every leaf pair
/// falls in bucket 0 or in bucket 2, and bucket 1 holds nothing. The AND diagnostic therefore says
/// nothing about it. That is the truth of this configuration and not an approximation of it. What
/// a caller does about it is the caller's business.
///
/// Every cell is set once and the tally is a sum of integers, so the order the cells are visited
/// in is free, and it is chosen for the species node state: that keeps its cliques contiguous
/// (ADR-0013), so a band of species leaves at one molecule leaf is one cache line of it, where a
/// field row read it a line per cell. Cells are named by subscript, which a node state laid out
/// either way turns into a position with no division.
template<typename Link, FieldStorage Field, BinaryStorage NodeState>
[[nodiscard]] field_math::TLinkCounters hold_tree_fields_at_the_field(const Field &field,
                                                                      NodeState &species_field,
                                                                      NodeState &molecule_field) {
	static_assert(NUMBER_OF_TREES == 2, "the field is a species leaf by a molecule leaf");
	// The leaf block of each node state is the field: rows are species leaves, columns molecule
	// leaves. Each node state knows one of the two counts; the other dimension is its own tree's.
	const size_t n_rows    = molecule_field.dimensions()[0];
	const size_t n_columns = species_field.dimensions()[1];
	if (n_rows * n_columns != field.total_size_of_container_space()) {
		throw coretools::TDevError(
		    "Cannot hold the tree fields at the field: their leaf blocks are ", n_rows, " x ",
		    n_columns, " and the field holds ", field.total_size_of_container_space(), " cells.");
	}

	field_math::TLinkCounters counters;
	std::vector<size_t> species_inserts;
	std::vector<size_t> molecule_inserts;
	for (size_t first_row = 0; first_row < n_rows; first_row += BAND_ROWS) {
		const size_t end_row = std::min(first_row + BAND_ROWS, n_rows);
		for (size_t column = 0; column < n_columns; ++column) {
			for (size_t row = first_row; row < end_row; ++row) {
				const bool y = field.is_one(row * n_columns + column);
				const IndexArray cell{row, column};
				// Only a cell that disagrees is written. A sparse node state does not hold every
				// cell, and inserting one that already reads as zero would grow it for nothing.
				if (species_field.is_one(cell) != y) {
					(void)write_or_defer(species_field, cell, y, species_inserts);
				}
				if (molecule_field.is_one(cell) != y) {
					(void)write_or_defer(molecule_field, cell, y, molecule_inserts);
				}
				counters.add(Link::bucket(y, y), y);
			}
		}
	}
	// The cells a sparse node state did not hold. No cell is visited twice above, so none was
	// read after it was deferred.
	for (const size_t index : species_inserts) { species_field.insert_one(index); }
	for (const size_t index : molecule_inserts) { molecule_field.insert_one(index); }
	return counters;
}

} // namespace leaf_layer_start
