//
// Which link between the tree fields and the field this build uses.
//

#pragma once

#include "field/TFieldMath.h"

// One alias, mirroring storages/storage_backend.h: the relationship between the two tree fields
// and the field is named in exactly one place, so swapping it is a one-line change and the seam is
// visible without reading the sampler.
//
// Unlike the storage there is only one implementation today, so this is a plain alias rather than
// a compile-definition switch. When a second link exists this is where the `#if` goes.
using TLinkPolicy = field_math::TAndLinkPolicy;

static_assert(field_math::LinkPolicy<TLinkPolicy>,
              "The selected link does not implement the link policy interface.");
