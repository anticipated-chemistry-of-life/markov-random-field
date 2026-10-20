//
// The tree update's moves, on one tree's stattools parameters.
//
// TTree builds one of these once per iteration and hands it to `tree_update::run`. It fills the
// `TreeMoves` concept (tree/tree_update.h). It adds the prior and Hastings terms that stattools
// owns, in a fixed association: `(ratio) + prior ratio` per clique, and
// `ratio + (prior ratio + prior ratio)` per pair. The chain's output depends on that association.
//
// The parameters' box type is TTree, so TTree is only declared here. The members are defined in
// TStattoolsTreeMoves.cpp, where TTree is complete.
//

#pragma once

#include "Types.h"
#include "stattools/ParametersObservations/TParameter.h"
#include "stattools/Updates/TPairIndexSampler.h"

#include <cstddef>
#include <utility>

class TTree;

/// A `TreeMoves` adapter over the alpha, log nu and binned branch-length parameters of one tree.
/// It owns this iteration's branch pairs, and draws them and the coins from the generator.
///
/// A parameter fixed with `--<name>.update false` needs no special case. stattools makes its
/// proposal a no-op and its decision a rejection.
class TStattoolsTreeMoves {
public:
	using TypeParamAlpha       = stattools::TParameter<SpecAlpha, TTree>;
	using TypeParamLogNu       = stattools::TParameter<SpecLogNu, TTree>;
	using TypeParamBinBranches = stattools::TParameter<SpecBinnedBranches, TTree>;

private:
	TypeParamAlpha *_alpha                = nullptr;
	TypeParamLogNu *_log_nu               = nullptr;
	TypeParamBinBranches *_branch_lengths = nullptr;
	size_t _n_cliques                     = 0;
	stattools::TPairIndexSampler _pairs;

public:
	TStattoolsTreeMoves(TypeParamAlpha &alpha, TypeParamLogNu &log_nu,
	                    TypeParamBinBranches &branch_lengths, size_t n_cliques, size_t n_branches)
	    : _alpha(&alpha), _log_nu(&log_nu), _branch_lengths(&branch_lengths),
	      _n_cliques(n_cliques), _pairs(n_branches) {}

	// stattools holds the proposal in `value()` while it is pending, and the accepted value
	// otherwise. So `proposed_*` and `accepted_*` read the same accessor, at different times.
	void propose_alphas();
	[[nodiscard]] double proposed_alpha(size_t clique) const;
	[[nodiscard]] double accepted_alpha(size_t clique) const;
	/// The ratio comes first and the prior ratio is added to it: `(ratio) + prior ratio`.
	bool decide_alpha(size_t clique, double log_ratio);

	[[nodiscard]] double accepted_log_nu(size_t clique) const;
	void propose_log_nus();
	[[nodiscard]] double proposed_log_nu(size_t clique) const;
	/// The ratio comes first and the prior ratio is added to it: `(ratio) + prior ratio`.
	bool decide_log_nu(size_t clique, double log_ratio);

	void sample_branch_pairs();
	[[nodiscard]] size_t n_branch_pairs() const;
	[[nodiscard]] std::pair<size_t, size_t> branch_pair(size_t pair) const;
	static bool branch_coin();
	[[nodiscard]] size_t current_bin(size_t branch) const;
	/// Sets both proposed bins, then decides the pair. The prior ratios are added together, and
	/// their sum is added to the ratio: `ratio + (prior ratio + prior ratio)`.
	bool decide_branch_pair(size_t pair, size_t first_bin, size_t second_bin, double log_ratio);
};
