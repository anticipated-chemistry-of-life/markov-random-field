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
#include "storages/storage_concepts.h"
#include <cstddef>

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
/// dimension per tree -- so a record's linear index is already the field's.
template<FieldStorage Field> void start_the_field_at(const Field &records, Field &field) {
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

/// Puts both tree fields at the field, and returns the six counters over the whole of it.
///
/// The tally is degenerate by construction. Both tree fields read the field, so every leaf pair
/// falls in bucket 0 or in bucket 2, and bucket 1 holds nothing. The AND diagnostic therefore says
/// nothing about it. That is the truth of this configuration and not an approximation of it. What
/// a caller does about it is the caller's business.
template<typename Link, FieldStorage Field, BinaryFieldStorage NodeState>
[[nodiscard]] field_math::TLinkCounters hold_tree_fields_at_the_field(const Field &field,
                                                                      NodeState &species_field,
                                                                      NodeState &molecule_field) {
	field_math::TLinkCounters counters;
	for (size_t index = 0; index < field.total_size_of_container_space(); ++index) {
		const bool y          = field.is_one(index);
		const IndexArray cell = field.get_multi_dimensional_index(index);
		// Only a cell that disagrees is written. A sparse node state does not hold every cell, and
		// inserting one that already reads as zero would grow it for nothing.
		for (auto *tree_field : {&species_field, &molecule_field}) {
			const size_t index_in_tree_field =
			    tree_field->get_linear_index_in_container_space(cell);
			if (tree_field->is_one(index_in_tree_field) == y) { continue; }
			if (y) {
				tree_field->insert_one(index_in_tree_field);
			} else {
				tree_field->insert_zero(index_in_tree_field);
			}
		}
		counters.add(Link::bucket(y, y), y);
	}
	return counters;
}

} // namespace leaf_layer_start
