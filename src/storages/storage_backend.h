//
// Which implementation of the field and of the internal state this build uses.
//

#pragma once

#include "storages/storage_concepts.h"

// The build sets exactly one ACOL_STORAGE_<BACKEND> macro (cmake option ACOL_STORAGE_BACKEND,
// which defaults to `sparse`). Selecting a backend is a compile definition rather than a runtime
// switch so that every storage access inlines: nothing here is ever reached through a virtual call.
//
// A second backend is one branch below plus one value accepted by the cmake option.
#if defined(ACOL_STORAGE_SPARSE)

#include "storages/y_storage/TStorageYMatrix.h"
#include "storages/z_storage/TStorageZMatrix.h"

using TFieldStorage         = TStorageYMatrix;
using TInternalStateStorage = TStorageZMatrix;

#else
#error                                                                                             \
    "No storage backend was selected. Configure with -DACOL_STORAGE_BACKEND=sparse (the default), which defines ACOL_STORAGE_SPARSE."
#endif

static_assert(FieldStorage<TFieldStorage>,
              "The selected field backend does not implement the field storage interface.");
static_assert(
    BinaryFieldStorage<TInternalStateStorage>,
    "The selected internal-state backend does not implement the binary storage interface.");
