//
// The data sources that take part in the block update through a pass of their own, as one
// compile-time list (ADR-0017).
//
// The block update reads which passes a build carries from this list alone. TDataModel names the
// data source that hands out each pass. The fold over the list is field/TBlockData.h. The list
// order is part of the chain (ADR-0017): LOTUS, then the simple error model, then mass
// spectrometry. LOTUS and mass spectrometry have no pass yet, so the block model multiplies them
// in at their own positions (field/TBlockModel.h).
//

#pragma once

#include "data_sources/data_source_pass.h"
#include "simple_error_model/block_pass.h"

/// A compile-time list of passes.
template<DataSourcePass... Passes> struct TPassList {};

#ifdef USE_SIMPLE_ERROR_MODEL
/// Every compiled-in pass, in the order the block update multiplies their cell factors.
using TBlockSources = TPassList<simple_error_model::TBlockPass>;
#else
using TBlockSources = TPassList<>;
#endif
