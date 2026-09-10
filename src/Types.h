//
// Created by madleina on 22.10.24.
//

#ifndef ACOL_TYPES_H
#define ACOL_TYPES_H

#include "coretools/Types/TStringHash.h"
#include "coretools/Types/commonWeakTypes.h"
#include "coretools/Types/probability.h"
#include "stattools/ParametersObservations/TObservation.h"
#include "stattools/ParametersObservations/TParameter.h"
#include "stattools/ParametersObservations/spec.h"
#include "stattools/Priors/TPriorBeta.h"
#include "stattools/Priors/TPriorExponential.h"
#include "stattools/Priors/TPriorGamma.h"
#include "stattools/Priors/TPriorNormal.h"
#include "stattools/Priors/TPriorUniform.h"
#include <cstdint>
#include <memory>
#include <vector>

class TTree; // forward declaration to avoid circular inclusion

// Which sources of information are compiled in? These are independent switches (the cmake options
// LOTUS and SIMPLE_DATA, i.e. the letters 'l' and 's' passed to `just build`), so a build can use
// either, both, or -- once MS data is wired up -- all three. Note that these constants are only
// for static_asserts and for logging: members and member functions of a data source must be
// guarded with #ifdef, not with `if constexpr`, because a discarded `if constexpr` branch still
// has to name-resolve.
#ifdef USE_LOTUS
constexpr static bool UseLotus = true;
#else
constexpr static bool UseLotus = false;
#endif

#ifdef USE_SIMPLE_ERROR_MODEL
constexpr static bool UseSimpleErrorModel = true;
#else
constexpr static bool UseSimpleErrorModel = false;
#endif

static_assert(UseLotus || UseSimpleErrorModel,
              "No source of data was compiled in: Y could not be informed by anything. Configure "
              "the build with at least one of -DLOTUS=ON or -DSIMPLE_DATA=ON, i.e. pass 'l' "
              "and/or 's' to `just build`.");

// Parameter types
using TypeGamma                    = coretools::StrictlyPositive;
using TypeErrorRate                = coretools::ZeroOneOpen;
using TypeEpsilonSimpleModel       = coretools::ZeroOneOpen;
using TypeAlpha                    = coretools::Probability;
using TypeLogNu                    = coretools::Unbounded;
using TypeNu                       = coretools::StrictlyPositive;
using TypeMeanLogNu                = coretools::Unbounded;
using TypeVarLogNu                 = coretools::StrictlyPositive;
using TypeBinnedBranchLengths      = coretools::UnsignedInt8WithMax<0>;
using TypeFilterProbability        = coretools::UnsignedInt8;
using TypeContaminationProbability = coretools::ZeroOneOpen;

// Gamma
// Weakly-informative Gamma(shape=alpha, rate=beta) prior on the detection rate.
// Counts enter as log(paper_count+1) (median c̄≈1.61), so gamma*=ln2/c̄≈0.43 gives
// detection 1-exp(-gamma·c̄)=0.5 at the median. Gamma(2, 4.6) -> mean=alpha/beta≈0.43,
// mode=(alpha-1)/beta≈0.22, with negligible mass on the saturation corner (gamma≫1)
// that lets the field collapse against detection. Hyperparameters set via
// --gamma.priorParameters "<alpha>,<beta>".
// TODO : verify this Claude bullshit
using PriorOnGamma = stattools::prior::TUniformFixed<stattools::TParameterBase, TypeGamma, 1>;
using SpecGamma =
    stattools::ParamSpec<TypeGamma, stattools::Hash<coretools::toHash("gamma")>, PriorOnGamma>;

// Epsilon
using PriorOnErrorRate =
    stattools::prior::TUniformFixed<stattools::TParameterBase, TypeErrorRate, 1>;
using SpecErrorRate =
    stattools::ParamSpec<TypeErrorRate, stattools::Hash<coretools::toHash("epsilon")>,
                         PriorOnErrorRate>;

// Epsilon of the simple error model
// The rate at which the simple error model data D misreports the latent state Y. Default prior is
// Beta(1, 1), i.e. uniform on (0, 1): the simple error model exists to diagnose the rest of the
// model, so its error rate should be driven by the likelihood alone and not pulled anywhere by the
// prior. Hyperparameters set via --epsilon_simple_model.priorParameters "<alpha>,<beta>".
// Note this is a *different* type from SpecErrorRate even though both are ZeroOneOpen + Beta: the
// name hash is a template parameter, which is what lets TDataModel overload calculateLLRatio on
// each of them.
using PriorOnEpsilonSimpleModel =
    stattools::prior::TUniformFixed<stattools::TParameterBase, TypeEpsilonSimpleModel, 1>;
using SpecEpsilonSimpleModel =
    stattools::ParamSpec<TypeEpsilonSimpleModel,
                         stattools::Hash<coretools::toHash("epsilon_simple_model")>,
                         PriorOnEpsilonSimpleModel>;

// Omega, the error probability
// The rate at which one tree field cell is corrupted before the two are reconciled into the field
// (ADR-0005). One scalar, shared by both trees.
//
// The support is the open interval (0, 0.5), which is a statement about the model rather than a
// range on an argument: at 0 the link is the deterministic AND and the block update takes log(0),
// and at 0.5 and above the tree fields are anti-correlated with the field. No coretools interval
// spells that, so the bounds are set at run time -- see
// TMarkovField::set_error_probability_support, which must run before stattools sizes the parameter.
//
// The prior is an exponential truncated to that interval, so it is concentrated near small values:
// an error probability near a half says the tree fields carry almost nothing. The default rate of
// 20 puts the prior mean at 0.05, and truncation moves it by 2e-5. Set it with
// --error_probability_prior_rate.
using TypeErrorProbability = coretools::MinMaxVariable<coretools::toHash("omega")>;
using PriorOnErrorProbability =
    stattools::prior::TExponentialFixed<stattools::TParameterBase, TypeErrorProbability, 1,
                                        /*TruncatedTop=*/true, /*TruncatedBottom=*/true>;
using SpecErrorProbability =
    stattools::ParamSpec<TypeErrorProbability, stattools::Hash<coretools::toHash("omega")>,
                         PriorOnErrorProbability>;

// Alpha
using PriorOnAlpha = stattools::prior::TUniformFixed<stattools::TParameterBase, TypeAlpha, 1>;
using SpecAlpha =
    stattools::ParamSpec<TypeAlpha, stattools::Hash<coretools::toHash("alpha")>, PriorOnAlpha>;

// Mean Nu
using PriorOnMeanLogNu =
    stattools::prior::TUniformFixed<stattools::TParameterBase, TypeMeanLogNu, 1>;
using SpecMeanLogNu =
    stattools::ParamSpec<TypeMeanLogNu, stattools::Hash<coretools::toHash("mean_log_nu")>,
                         PriorOnMeanLogNu>;

// Var Nu
using PriorOnVarLogNu = stattools::prior::TGammaFixed<stattools::TParameterBase, TypeVarLogNu, 1>;
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

// The likelihood box that anchors the model in the stattools DAG. It owns the Markov field and
// every compiled-in data source; see TDataModel.h.
class TDataModel; // forward declaration to avoid circular inclusion

// Markov Field (only needed for stattools purposes)
using TypeMarkovField = coretools::Boolean;
constexpr static size_t NumDimMarkovField =
    1; // note: only for stattools, actually not known at compile time
using PriorOnMarkovField = TTree;
using SpecMarkovField =
    stattools::ParamSpec<TypeMarkovField, stattools::Hash<coretools::toHash("MRF")>,
                         PriorOnMarkovField>;
using ParamMarkovField  = stattools::TParameter<SpecMarkovField, TDataModel>;
using MarkovFieldParams = std::vector<std::unique_ptr<ParamMarkovField>>;

// Observation anchoring TDataModel in the DAG. This is a "fake" observation: the real data lives
// in the individual data sources, but stattools requires a box to sit above an observation.
using TypeDataObs                     = coretools::Boolean;
constexpr static size_t NumDimDataObs = 2;
using StorageDataObs = coretools::TMultiDimensionalStorage<TypeDataObs, NumDimDataObs>;
using SpecDataObs    = stattools::TObservation<TypeDataObs, NumDimDataObs, TDataModel>;

// Observations: Mass Spec
class TMSMSData;
using TypeMSData                     = coretools::Boolean;
constexpr static size_t NumDimMSData = 1;
using StorageMSData = coretools::TMultiDimensionalStorage<TypeMSData, NumDimMSData>;
using SpecMSData    = stattools::TObservation<TypeMSData, NumDimMSData, TMSMSData>;

#endif // ACOL_TYPES_H
