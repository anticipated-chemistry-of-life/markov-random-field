//
// Created by madleina on 22.10.24.
//

#ifndef ACOL_TYPES_H
#define ACOL_TYPES_H

#include "coretools/Types/commonWeakTypes.h"
#include "coretools/Types/probability.h"
#include "stattools/ParametersObservations/THash.h"
#include "stattools/ParametersObservations/TObservation.h"
#include "stattools/ParametersObservations/spec.h"
#include "stattools/Priors/TPriorExponential.h"
#include "stattools/Priors/TPriorNormal.h"
#include "stattools/Priors/TPriorUniform.h"
#include <cstdint>

class TTree; // forward declaration to avoid circular inclusion
inline size_t NUMBER_OF_THREADS;
inline bool SIMULATION_NO_Z_INITIALIZATION;
inline bool SIMULATION_NO_Y_INITIALIZATION;
inline bool WRITE_Y;
inline bool WRITE_Y_TRACE;
inline bool WRITE_Z;
inline bool WRITE_Z_TRACE;
inline bool WRITE_BRANCH_LENGTHS;
inline bool WRITE_JOINT_LOG_PROB_DENSITY;

// use simple error model for Lotus?
#ifdef USE_LOTUS
constexpr static bool UseSimpleErrorModel = false;
#else
constexpr static bool UseSimpleErrorModel = true;
#endif

// Parameter types
using TypeGamma               = coretools::Positive;
using TypeErrorRate           = coretools::Positive;
using TypeAlpha               = coretools::Probability;
using TypeLogNu               = coretools::Unbounded;
using TypeNu                  = coretools::StrictlyPositive;
using TypeMeanLogNu           = coretools::Unbounded;
using TypeVarLogNu            = coretools::StrictlyPositive;
using TypeBinnedBranchLengths = coretools::UnsignedInt8WithMax<0>;

// Gamma
using PriorOnGamma = stattools::prior::TUniformFixed<TypeGamma>;
using SpecGamma    = stattools::ParamSpec<TypeGamma, stattools::name("gamma"), PriorOnGamma>;

// Epsilon
using PriorOnErrorRate = stattools::prior::TExponentialFixed<TypeErrorRate>;
using SpecErrorRate    = stattools::ParamSpec<TypeErrorRate, stattools::name("epsilon"), PriorOnErrorRate,
                                              stattools::EnforceUniqueHash<false>>;

// Alpha
using PriorOnAlpha = stattools::prior::TUniformFixed<TypeAlpha>;
using SpecAlpha =
    stattools::ParamSpec<TypeAlpha, stattools::name("alpha"), PriorOnAlpha, stattools::EnforceUniqueHash<false>>;

// Mean Nu
using PriorOnMeanLogNu = stattools::prior::TUniformFixed<TypeMeanLogNu>;
using SpecMeanLogNu    = stattools::ParamSpec<TypeMeanLogNu, stattools::name("mean_log_nu"), PriorOnMeanLogNu,
                                              stattools::EnforceUniqueHash<false>>;

// Var Nu
using PriorOnVarLogNu = stattools::prior::TExponentialFixed<TypeVarLogNu>;
using SpecVarLogNu    = stattools::ParamSpec<TypeVarLogNu, stattools::name("var_log_nu"), PriorOnVarLogNu,
                                             stattools::EnforceUniqueHash<false>>;

// Log Nu
using PriorOnLogNu = stattools::prior::TNormalInferred<SpecMeanLogNu, SpecVarLogNu, TypeLogNu>;
using SpecLogNu =
    stattools::ParamSpec<TypeLogNu, stattools::name("log_nu"), PriorOnLogNu, stattools::EnforceUniqueHash<false>>;

// binned branch lengths
using PriorOnBinnedBranches = stattools::prior::TUniformFixed<TypeBinnedBranchLengths>;
using SpecBinnedBranches    = stattools::ParamSpec<TypeBinnedBranchLengths, stattools::name("bin_branch"),
                                                   PriorOnBinnedBranches, stattools::EnforceUniqueHash<false>>;

// Markov Field (only needed for stattools purposes)
using TypeMarkovField                     = coretools::Boolean;
constexpr static size_t NumDimMarkovField = 1; // note: only for stattools, actually not known at compile time
using PriorOnMarkovField                  = TTree;
using SpecMarkovField = stattools::ParamSpec<TypeMarkovField, stattools::name("MRF"), PriorOnMarkovField,
                                             stattools::EnforceUniqueHash<false>>;

// Observation: Lotus
class TLotus; // forward declaration to avoid circular inclusion
using TypeLotus                     = coretools::Boolean;
constexpr static size_t NumDimLotus = 2;
using StorageLotus                  = coretools::TMultiDimensionalStorage<TypeLotus, NumDimLotus>;
using SpecLotus = stattools::TObservation<TypeLotus, stattools::name("lotus_obs"), NumDimLotus, TLotus>;

// Type for calculating the number of 1's per clique
using TypeCounter1 = uint32_t;

#endif // ACOL_TYPES_H
