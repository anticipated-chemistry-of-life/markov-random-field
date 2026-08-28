//
// The dense internal state.
//

#pragma once

#include "storages/TDenseStateArray.h"

/// The dense internal state is the dense state array and nothing else: `Z` carries one binary
/// state per cell and no posterior counter, so there is nothing for this to add. The name exists
/// so that the two backends spell their internal state the same way -- `TStorageZMatrix` beside
/// `TStorageZDense` -- and so that the sharing the dense field relies on is visible here rather
/// than buried in a member declaration.
///
/// It is the same type, so TDenseStateArray.h's own assertion that it satisfies
/// `BinaryFieldStorage` is this file's assertion too; repeating it here would check nothing new.
using TStorageZDense = TDenseStateArray;
