//
// Created by madleina on 22.10.24.
//

#ifndef ACOL_TYPES_H
#define ACOL_TYPES_H

#include "coretools/Types/TStringHash.h"
#include "coretools/Types/commonWeakTypes.h"
#include "coretools/Types/probability.h"
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
using TypeGamma                    = coretools::StrictlyPositive;
using TypeErrorRate                = coretools::ZeroOneOpen;
using TypeAlpha                    = coretools::Probability;
using TypeLogNu                    = coretools::Unbounded;
using TypeNu                       = coretools::StrictlyPositive;
using TypeMeanLogNu                = coretools::Unbounded;
using TypeVarLogNu                 = coretools::StrictlyPositive;
using TypeBinnedBranchLengths      = coretools::UnsignedInt8WithMax<0>;
using TypeFilterProbability        = coretools::Probability;
using TypeContaminationProbability = coretools::ZeroOneOpen;

// Gamma
// Weakly-informative Gamma(shape=alpha, rate=beta) prior on the detection rate.
// Counts enter as log(paper_count+1) (median c̄≈1.61), so gamma*=ln2/c̄≈0.43 gives
// detection 1-exp(-gamma·c̄)=0.5 at the median. Gamma(2, 4.6) -> mean=alpha/beta≈0.43,
// mode=(alpha-1)/beta≈0.22, with negligible mass on the saturation corner (gamma≫1)
// that lets the field collapse against detection. Hyperparameters set via
// --gamma.priorParameters "<alpha>,<beta>".
// TODO : verify this Claude bullshit
using PriorOnGamma = stattools::prior::TGammaFixed<stattools::TParameterBase, TypeGamma, 1>;
using SpecGamma =
    stattools::ParamSpec<TypeGamma, stattools::Hash<coretools::toHash("gamma")>, PriorOnGamma>;

// Epsilon
using PriorOnErrorRate = stattools::prior::TBetaFixed<stattools::TParameterBase, TypeErrorRate, 1>;
using SpecErrorRate =
    stattools::ParamSpec<TypeErrorRate, stattools::Hash<coretools::toHash("epsilon")>,
                         PriorOnErrorRate>;

// Alpha
using PriorOnAlpha = stattools::prior::TUniformFixed<stattools::TParameterBase, TypeAlpha, 1>;
using SpecAlpha =
    stattools::ParamSpec<TypeAlpha, stattools::Hash<coretools::toHash("alpha")>, PriorOnAlpha>;

// Mean Nu
// Weakly-informative Normal(mean, var) prior on mean_log_nu. Branch lengths are
// normalized to mean 1 and the likelihood depends on the product nu·t, so nu≈O(1)
// (mean_log_nu≈0) is the natural scale: the field is neither frozen (nu→0, constant
// cliques) nor maximally noisy. var=3 (sd≈1.73) keeps it weak while excluding the
// nu→0 collapse (mean_log_nu≈-9) at ~5 sd. Hyperparameters set via
// --<tree>_mean_log_nu.priorParameters "<mean>,<var>".
// TODO : verify this Claude bullshit
using PriorOnMeanLogNu =
    stattools::prior::TNormalFixed<stattools::TParameterBase, TypeMeanLogNu, 1>;
using SpecMeanLogNu =
    stattools::ParamSpec<TypeMeanLogNu, stattools::Hash<coretools::toHash("mean_log_nu")>,
                         PriorOnMeanLogNu>;

// Var Nu
using PriorOnVarLogNu =
    stattools::prior::TExponentialFixed<stattools::TParameterBase, TypeVarLogNu, 1>;
using SpecVarLogNu =
    stattools::ParamSpec<TypeVarLogNu, stattools::Hash<coretools::toHash("var_log_nu")>,
                         PriorOnVarLogNu>;

// Log Nu
using PriorOnLogNu = stattools::prior::TNormalInferred<stattools::TParameterBase, TypeLogNu, 1,
                                                       SpecMeanLogNu, SpecVarLogNu>;
using SpecLogNu =
    stattools::ParamSpec<TypeLogNu, stattools::Hash<coretools::toHash("log_nu")>, PriorOnLogNu>;

// binned branch lengths
using PriorOnBinnedBranches =
    stattools::prior::TUniformFixed<stattools::TParameterBase, TypeBinnedBranchLengths, 1>;
using SpecBinnedBranches =
    stattools::ParamSpec<TypeBinnedBranchLengths, stattools::Hash<coretools::toHash("bin_branch")>,
                         PriorOnBinnedBranches>;

// Probability to pass mass spec filter
using PriorOnMassSpecFilter =
    stattools::prior::TUniformFixed<stattools::TParameterBase, TypeFilterProbability, 1>;
using SpecMassSpecFilter =
    stattools::ParamSpec<TypeFilterProbability, stattools::Hash<coretools::toHash("filter_proba")>,
                         PriorOnMassSpecFilter>;

// Contamination probability in MassSpec : Y = 0 and MSData = 1
using PriorOnContaminationProba =
    stattools::prior::TBetaFixed<stattools::TParameterBase, TypeContaminationProbability, 1>;
using SpecContaminationProba =
    stattools::ParamSpec<TypeContaminationProbability,
                         stattools::Hash<coretools::toHash("contamination_proba")>,
                         PriorOnContaminationProba>;

// Markov Field (only needed for stattools purposes)
using TypeMarkovField = coretools::Boolean;
constexpr static size_t NumDimMarkovField =
    1; // note: only for stattools, actually not known at compile time
using PriorOnMarkovField = TTree;
using SpecMarkovField =
    stattools::ParamSpec<TypeMarkovField, stattools::Hash<coretools::toHash("MRF")>,
                         PriorOnMarkovField>;

// Observation: Lotus
class TLotus; // forward declaration to avoid circular inclusion
using TypeLotus                     = coretools::Boolean;
constexpr static size_t NumDimLotus = 2;
using StorageLotus                  = coretools::TMultiDimensionalStorage<TypeLotus, NumDimLotus>;
using SpecLotus                     = stattools::TObservation<TypeLotus, NumDimLotus, TLotus>;

// Observations: Mass Spec
class TMSMSData;
using TypeMSData                     = coretools::Boolean;
constexpr static size_t NumDimMSData = 1;
using StorageMSData = coretools::TMultiDimensionalStorage<TypeMSData, NumDimMSData>;
using SpecMSData    = stattools::TObservation<TypeMSData, NumDimMSData, TMSMSData>;

// Type for calculating the number of 1's per clique
using TypeCounter1 = uint32_t;

#endif // ACOL_TYPES_H
