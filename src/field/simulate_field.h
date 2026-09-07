//
// The field a simulated chain draws, from the two tree fields and the link.
//
// Under this model the field has a forward draw: each tree draws its own node state top-down
// (tree/node_state_draw.h), and the field is then a noisy AND of the two tree fields, cell by cell
// (ADR-0005). Nothing iterates.
//
// The six counters come back with the field, because the draw already knows every leaf pair's
// bucket and field state. Recounting them afterwards would be a second walk over the same cells,
// and a second place the bucket arithmetic is written.
//
// Like field/leaf_layer_start.h this reads and writes through the storage concepts and holds no
// tree, so one suite asserts it over every backend pairing.
//

#pragma once

#include "constants.h"
#include "coretools/Main/TError.h"
#include "field/TFieldMath.h"
#include "random/TCellUniforms.h"
#include "storages/storage_concepts.h"

#include <concepts>
#include <cstddef>

namespace simulate_field {

/// Draws every field cell from the two tree fields, and returns the six counters over the whole
/// field.
///
/// The field must be empty. This writes ones and never zeros, because a cell the link left at zero
/// already reads as zero, and inserting that zero would grow a sparse field for nothing. A field
/// that already held states would keep them, so the postcondition rests on the field being empty
/// rather than on the caller.
///
/// A tree field and the field are addressed at the same `(row, column)` for a given leaf pair, so
/// each tree field cell is found by the field's own subscript rather than by a conversion.
template<typename Link, FieldStorage Field, BinaryFieldStorage NodeState, CellUniforms Uniforms>
[[nodiscard]] field_math::TLinkCounters
draw_from_the_tree_fields(Field &field, const NodeState &species_field,
                          const NodeState &molecule_field,
                          const field_math::TErrorProbability &omega, const Uniforms &uniforms) {
	if (!field.empty()) {
		throw coretools::TDevError(
		    "Cannot draw the field from the tree fields: the field already holds states.");
	}

	field_math::TLinkCounters counters;
	for (size_t index = 0; index < field.total_size_of_container_space(); ++index) {
		const IndexArray cell = field.get_multi_dimensional_index(index);
		const bool z_s =
		    species_field.is_one(species_field.get_linear_index_in_container_space(cell));
		const bool z_m =
		    molecule_field.is_one(molecule_field.get_linear_index_in_container_space(cell));

		const bool y = uniforms.at(index) < Link::prob_y_is_one(z_s, z_m, omega);
		if (y) { field.insert_one(index); }
		counters.add(Link::bucket(z_s, z_m), y);
	}
	return counters;
}

} // namespace simulate_field
