//
// The data sources that take part in the block update through a pass of their own, as one
// compile-time list (ADR-0017).
//
// The block update reads which passes a build carries from this list alone. TDataModel names the
// data source that hands out each pass. The fold over the list is field/TBlockData.h. The list
// order is part of the chain (ADR-0017): LOTUS, then the simple error model, then mass
// spectrometry.
//

#pragma once

#include "data_sources/data_source_pass.h"
#include "lotus/block_pass.h"
#include "mass_spec/block_pass.h"
#include "simple_error_model/block_pass.h"

/// A compile-time list of passes.
template<DataSourcePass... Passes> struct TPassList {};

/// The passes of every list, one list after the other.
template<typename... Lists> struct TPassListConcat;
template<DataSourcePass... Passes> struct TPassListConcat<TPassList<Passes...>> {
	using type = TPassList<Passes...>;
};
template<DataSourcePass... First, DataSourcePass... Second, typename... Rest>
struct TPassListConcat<TPassList<First...>, TPassList<Second...>, Rest...>
    : TPassListConcat<TPassList<First..., Second...>, Rest...> {};

#ifdef USE_LOTUS
using TLotusPasses = TPassList<lotus_math::TBlockPass>;
#else
using TLotusPasses = TPassList<>;
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
using TSimpleErrorPasses = TPassList<simple_error_model::TBlockPass>;
#else
using TSimpleErrorPasses = TPassList<>;
#endif
#ifdef USE_MS_DATA
using TMassSpecPasses = TPassList<mass_spec_math::TBlockPass>;
#else
using TMassSpecPasses = TPassList<>;
#endif

/// Every compiled-in pass, in the order the block update multiplies their cell factors.
using TBlockSources =
    typename TPassListConcat<TLotusPasses, TSimpleErrorPasses, TMassSpecPasses>::type;
