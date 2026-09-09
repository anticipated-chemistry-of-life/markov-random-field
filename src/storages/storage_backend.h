//
// Which storage backs the field, which backs the node state, and which backs the observed data.
//

#pragma once

#include "storages/storage_concepts.h"

#include "storages/TDenseStateArray.h"
#include "storages/TSparseBinaryArray.h"
#include "storages/y_storage/TStorageYDense.h"
#include "storages/y_storage/TStorageYMatrix.h"
#include "storages/z_storage/TStorageZDense.h"
#include "storages/z_storage/TStorageZMatrix.h"

#include <type_traits>

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

/// The storage the simple error model data takes.
///
/// It makes no choice of its own. The field's flag decides both, because the two are the same
/// shape and are read cell for cell against each other. A run that wants the sparse field wants
/// this sparse too. A third define would also be a third pairing for the parity gate to cover.
using TBinaryStorage = std::conditional_t<std::is_same_v<TFieldStorage, TStorageYDense>,
                                          TDenseStateArray, TSparseBinaryArray>;

// The alias above reads the field's choice as a type, so it has to know every type that choice can
// be. A define naming a third field storage would fall to the sparse side without saying so.
static_assert(std::is_same_v<TFieldStorage, TStorageYDense> ||
                  std::is_same_v<TFieldStorage, TStorageYMatrix>,
              "The binary storage follows the field, so the field has to be one of the two "
              "storages it knows. Add the new one to TBinaryStorage before selecting it here.");

// The LOTUS records take no alias at all. They are always a TSparseBinaryArray, whatever the build
// selects above. A run reads them in once and never writes them again. They hold a few ones per
// thousand cells under any backend.

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

// The counted storage is the field and nothing else. Everything else the sampler holds is a plain
// binary storage: a state per cell, and no posterior counter that is never incremented and never
// read.
static_assert(BinaryStorage<TNodeStateStorage> && !FieldStorage<TNodeStateStorage>,
              "The selected node-state storage must be a binary storage and not a field.");
static_assert(BinaryStorage<TBinaryStorage> && !FieldStorage<TBinaryStorage>,
              "The selected binary storage must be a binary storage and not a field.");
// The LOTUS records are the third, and TSparseBinaryArray.h asserts them where the type is
// defined, because they are pinned to that type rather than selected here.
