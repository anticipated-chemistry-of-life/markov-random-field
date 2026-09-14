//
// The seam every per-cell data source satisfies: an independent likelihood term over the field,
// observing Y one cell at a time (ADR-0005's data factor, `log p(L, D | Y)`). LOTUS and the simple
// error model are this shape. Mass spectrometry data is not -- it scores an assignment problem
// over a species' runs, not a per-cell state (see mass_spec/msms_data.h) -- so it is not named by
// ACOL_DATA_SOURCES and does not satisfy DataSource below.
//
// ACOL_DATA_SOURCES(X) lists whichever per-cell sources are compiled in, once, so that
// TDataModel's six uniform per-source sites -- initialize, guess_initial_values, the likelihood
// sum, notifier-stats collection, simulate + write, and the accessor pair -- expand from this one
// list instead of being hand-mirrored at each site (TDataModel.h/.cpp). Adding a per-cell source
// is: write a class that satisfies DataSource, add one row here.
//
// What is *not* generated from this list: each source's own member declaration, its constructor
// call, and the stattools parameter-move dispatch (TDataModel::calculateLLRatio / updateTempVals).
// A source's parameter count and constructor arguments vary source to source, and stattools
// dispatches a move by the parameter's own type, not by its owning source -- both stay hand-written
// in TDataModel.h/.cpp and TCore.h/.cpp, the way they are today.
//
#pragma once

#include "ntfy/TNtfyNotifier.h"
#include "storages/storage_backend.h"
#include <concepts>
#include <string>
#include <vector>

class TDataModel;

/// Everything the three MCMC-progress notifications (one burn-in round finished, burn-in
/// finished, MCMC finished) report, collected once per notification so each compiled-in data
/// source can add to it without the collection site (TDataModel::_collect_notifier_stats) having
/// to know which sources exist.
struct TNotifierStats {
	std::vector<std::string> dim_names;
	std::vector<TNtfyNotifier::ParamStats> gamma_stats;
	std::vector<TNtfyNotifier::NamedStats> scalar_stats;
};

/// A per-cell data source. `T` must, given a box to hang its parameters off, a flag saying whether
/// this run simulates, the field, and an output prefix:
///
/// - size itself and its parameters, and read its data file when not simulating (`initialize`);
/// - set its parameters to their configured starting values and derive whatever cached state its
///   likelihood needs from the current field (`guess_initial_values`);
/// - report its current log-likelihood (`log_likelihood`);
/// - draw its data from a simulated field and write it out (`simulate_from_Y`, `write_simulated`);
/// - append whatever it has to report to a shared stats collector (`contribute_stats`).
///
/// Not part of this concept: the per-cell hook the block update reads (TLotus::calculate_LL_update_Y,
/// TSimpleErrorModel::probabilities_for_Y_update) -- their shapes still differ source to source --
/// and the parameter-move dispatch, which stattools forces to be one overload per parameter type.
template <typename T>
concept DataSource =
    requires(T &s, TDataModel *box, bool simulate, const TFieldStorage &Y,
             const std::string &prefix, TNotifierStats &stats) {
	    { s.initialize(box, simulate) } -> std::same_as<void>;
	    { s.guess_initial_values(Y) } -> std::same_as<void>;
	    { s.log_likelihood() } -> std::convertible_to<double>;
	    { s.simulate_from_Y(Y) } -> std::same_as<void>;
	    { s.write_simulated(prefix) } -> std::same_as<void>;
	    { s.contribute_stats(stats) } -> std::same_as<void>;
    };

// The list of compiled-in per-cell data sources, one X(member_name, Type) row per source. A
// consumer that wants "for each compiled-in source, do X" defines its own X and expands
// ACOL_DATA_SOURCES(X) -- see TDataModel.h/.cpp for the six sites this drives. `member_name` is the
// source's field name without the leading underscore, so a row here names both `TDataModel::_lotus`
// and the accessor `get_lotus()`.
#if defined(USE_LOTUS) && defined(USE_SIMPLE_ERROR_MODEL)
#define ACOL_DATA_SOURCES(X) X(lotus, TLotus) X(simple_error_model, TSimpleErrorModel)
#elif defined(USE_LOTUS)
#define ACOL_DATA_SOURCES(X) X(lotus, TLotus)
#elif defined(USE_SIMPLE_ERROR_MODEL)
#define ACOL_DATA_SOURCES(X) X(simple_error_model, TSimpleErrorModel)
#else
#define ACOL_DATA_SOURCES(X)
#endif
