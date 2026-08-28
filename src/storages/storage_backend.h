//
// Which implementation of the field and of the internal state this build uses.
//

#pragma once

#include "storages/storage_concepts.h"

// The build sets exactly one ACOL_STORAGE_<BACKEND> macro (cmake option ACOL_STORAGE_BACKEND,
// which defaults to `sparse`). Selecting a backend is a compile definition rather than a runtime
// switch so that every storage access inlines: nothing here is ever reached through a virtual call.
//
// The two branches below are the whole of the choice. Because it is made at compile time, checking
// that the backends agree is a comparison of two *binaries* rather than of two runtime modes --
// which is what `just parity` (tests/backend_parity/) does, and it therefore exercises exactly
// what ships.
#if defined(ACOL_STORAGE_SPARSE)

#include "storages/y_storage/TStorageYMatrix.h"
#include "storages/z_storage/TStorageZMatrix.h"

using TFieldStorage         = TStorageYMatrix;
using TInternalStateStorage = TStorageZMatrix;

#elif defined(ACOL_STORAGE_DENSE)

#include "storages/y_storage/TStorageYDense.h"
#include "storages/z_storage/TStorageZDense.h"

using TFieldStorage         = TStorageYDense;
using TInternalStateStorage = TStorageZDense;

#else
#error                                                                                             \
    "No storage backend was selected. Configure with -DACOL_STORAGE_BACKEND=sparse (the default) or =dense, which define ACOL_STORAGE_SPARSE and ACOL_STORAGE_DENSE."
#endif

static_assert(FieldStorage<TFieldStorage>,
              "The selected field backend does not implement the field storage interface.");
static_assert(
    BinaryFieldStorage<TInternalStateStorage>,
    "The selected internal-state backend does not implement the binary storage interface.");
