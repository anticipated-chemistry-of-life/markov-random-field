//
// Created by madleina on 22.10.24.
//

#ifndef ACOL_TYPES_H
#define ACOL_TYPES_H

#include "coretools/Types/commonWeakTypes.h"
#include "stattools/ParametersObservations/THash.h"
#include "stattools/ParametersObservations/TObservation.h"
#include "stattools/ParametersObservations/spec.h"
#include "stattools/Priors/TPriorUniform.h"
#include <stdint.h>

// Parameter types
using TypeGamma = coretools::Positive;
using TypeDelta = coretools::Positive;

// Gamma
using PriorOnGamma = stattools::prior::TUniformFixed<TypeGamma>;
using SpecGamma    = stattools::ParamSpec<TypeGamma, stattools::name("gamma"), PriorOnGamma>;

// Delta
using PriorOnDelta = stattools::prior::TUniformFixed<TypeDelta>;
using SpecDelta    = stattools::ParamSpec<TypeDelta, stattools::name("delta"), PriorOnDelta>;

// Observation: Lotus
class TLotus; // forward declaration to avoid circular inclusion
using TypeLotus                     = coretools::Boolean;
constexpr static size_t NumDimLotus = 2;
using StorageLotus                  = coretools::TMultiDimensionalStorage<TypeLotus, NumDimLotus>;
using SpecLotus                     = stattools::TObservation<TypeLotus, stattools::name("lotus"), NumDimLotus, TLotus>;

/**
 * Type for the number of bins for the branches.
 */
using TypeBinBranches = uint8_t;

#endif // ACOL_TYPES_H
