//
// Which storage backs the field, and which backs the node state.
//

#pragma once

#include "storages/storage_concepts.h"

#include "storages/y_storage/TStorageYDense.h"
#include "storages/y_storage/TStorageYMatrix.h"
#include "storages/z_storage/TStorageZDense.h"
#include "storages/z_storage/TStorageZMatrix.h"

// The field and the node state choose their storage independently. Each choice is one alias below.
// Changing one is an edit to one line. No build option and no build directory take part.
//
// The choice is a type alias and not a runtime switch, so every storage access inlines. Nothing
// here is reached through a virtual call. storage_concepts.h states the interface an alias has to
// satisfy, and the foot of this file asserts it.
//
// An external define wins over the alias it guards. Passing
// `-DACOL_FIELD_STORAGE=TStorageYDense` to the compiler is how `just parity`
// (tests/backend_parity/) builds two binaries from one source tree.
//
// Both defaults are dense for now, which is one of the two pairs the gate covers. ADR-0006 argues
// for a sparse field against a dense node state, on fill rather than size. That pairing is one
// line away when the runs need it.

#ifdef ACOL_FIELD_STORAGE
using TFieldStorage = ACOL_FIELD_STORAGE;
#else
using TFieldStorage = TStorageYDense;
#endif

#ifdef ACOL_NODE_STATE_STORAGE
using TNodeStateStorage = ACOL_NODE_STATE_STORAGE;
#else
using TNodeStateStorage = TStorageZDense;
#endif

// All four pairings compile. Continuous integration gates two of them:
//
//     field   node state   gated
//     ------  -----------  ---------------------
//     sparse  sparse       yes, by `just parity`
//     dense   dense        yes, by `just parity`  <- the default
//     sparse  dense        no
//     dense   sparse       no
//
// The two gated pairs between them exercise both storages. An ungated pair is untested, not
// unsupported. ADR-0006 gives the argument.

static_assert(FieldStorage<TFieldStorage>,
              "The selected field storage does not implement the field storage interface.");
static_assert(BinaryFieldStorage<TNodeStateStorage>,
              "The selected node-state storage does not implement the binary storage interface.");
