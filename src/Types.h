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
#include "stattools/Priors/TPriorBeta.h"
#include "stattools/Priors/TPriorExponential.h"
#include "stattools/Priors/TPriorGamma.h"
#include "stattools/Priors/TPriorNormal.h"
#include "stattools/Priors/TPriorUniform.h"
#include <cstdint>

class TTree; // forward declaration to avoid circular inclusion

// use simple error model for Lotus?
#ifdef USE_LOTUS
constexpr static bool UseSimpleErrorModel = false;
#else
constexpr static bool UseSimpleErrorModel = true;
#endif

// Parameter types
using TypeGamma               = coretools::StrictlyPositive;
using TypeErrorRate           = coretools::ZeroOneOpen;
using TypeAlpha               = coretools::Probability;
using TypeLogNu               = coretools::Unbounded;
using TypeNu                  = coretools::StrictlyPositive;
using TypeMeanLogNu           = coretools::Unbounded;
using TypeVarLogNu            = coretools::StrictlyPositive;
using TypeBinnedBranchLengths = coretools::UnsignedInt8WithMax<0>;

// Gamma
// Weakly-informative Gamma(shape=alpha, rate=beta) prior on the detection rate.
// Counts enter as log(paper_count+1) (median c̄≈1.61), so gamma*=ln2/c̄≈0.43 gives
// detection 1-exp(-gamma·c̄)=0.5 at the median. Gamma(2, 4.6) -> mean=alpha/beta≈0.43,
// mode=(alpha-1)/beta≈0.22, with negligible mass on the saturation corner (gamma≫1)
// that lets the field collapse against detection. Hyperparameters set via
// --gamma.priorParameters "<alpha>,<beta>".
// TODO : verify this Claude bullshit
using PriorOnGamma = stattools::prior::TGammaFixed<TypeGamma>;
using SpecGamma    = stattools::ParamSpec<TypeGamma, stattools::name("gamma"), PriorOnGamma>;

// Epsilon
using PriorOnErrorRate = stattools::prior::TBetaFixed<TypeErrorRate>;
using SpecErrorRate = stattools::ParamSpec<TypeErrorRate, stattools::name("epsilon"),
                                           PriorOnErrorRate, stattools::EnforceUniqueHash<false>>;

// Alpha
using PriorOnAlpha = stattools::prior::TUniformFixed<TypeAlpha>;
using SpecAlpha    = stattools::ParamSpec<TypeAlpha, stattools::name("alpha"), PriorOnAlpha,
                                          stattools::EnforceUniqueHash<false>>;

// Mean Nu
// Weakly-informative Normal(mean, var) prior on mean_log_nu. Branch lengths are
// normalized to mean 1 and the likelihood depends on the product nu·t, so nu≈O(1)
// (mean_log_nu≈0) is the natural scale: the field is neither frozen (nu→0, constant
// cliques) nor maximally noisy. var=3 (sd≈1.73) keeps it weak while excluding the
// nu→0 collapse (mean_log_nu≈-9) at ~5 sd. Hyperparameters set via
// --<tree>_mean_log_nu.priorParameters "<mean>,<var>".
// TODO : verify this Claude bullshit
using PriorOnMeanLogNu = stattools::prior::TNormalFixed<TypeMeanLogNu>;
using SpecMeanLogNu = stattools::ParamSpec<TypeMeanLogNu, stattools::name("mean_log_nu"),
                                           PriorOnMeanLogNu, stattools::EnforceUniqueHash<false>>;

// Var Nu
using PriorOnVarLogNu = stattools::prior::TExponentialFixed<TypeVarLogNu>;
using SpecVarLogNu    = stattools::ParamSpec<TypeVarLogNu, stattools::name("var_log_nu"),
                                             PriorOnVarLogNu, stattools::EnforceUniqueHash<false>>;

// Log Nu
using PriorOnLogNu = stattools::prior::TNormalInferred<SpecMeanLogNu, SpecVarLogNu, TypeLogNu>;
using SpecLogNu    = stattools::ParamSpec<TypeLogNu, stattools::name("log_nu"), PriorOnLogNu,
                                          stattools::EnforceUniqueHash<false>>;

// binned branch lengths
using PriorOnBinnedBranches = stattools::prior::TUniformFixed<TypeBinnedBranchLengths>;
using SpecBinnedBranches =
    stattools::ParamSpec<TypeBinnedBranchLengths, stattools::name("bin_branch"),
                         PriorOnBinnedBranches, stattools::EnforceUniqueHash<false>>;

// Markov Field (only needed for stattools purposes)
using TypeMarkovField = coretools::Boolean;
constexpr static size_t NumDimMarkovField =
    1; // note: only for stattools, actually not known at compile time
using PriorOnMarkovField = TTree;
using SpecMarkovField =
    stattools::ParamSpec<TypeMarkovField, stattools::name("MRF"), PriorOnMarkovField,
                         stattools::EnforceUniqueHash<false>>;

// Observation: Lotus
class TLotus; // forward declaration to avoid circular inclusion
using TypeLotus                     = coretools::Boolean;
constexpr static size_t NumDimLotus = 2;
using StorageLotus                  = coretools::TMultiDimensionalStorage<TypeLotus, NumDimLotus>;
using SpecLotus =
    stattools::TObservation<TypeLotus, stattools::name("lotus_obs"), NumDimLotus, TLotus>;

// Type for calculating the number of 1's per clique
using TypeCounter1 = uint32_t;

#endif // ACOL_TYPES_H
