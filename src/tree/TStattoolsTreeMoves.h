//
// The tree update's moves, on one tree's stattools parameters.
//
// TTree builds one of these once per iteration and hands it to `tree_update::run`. It fills the
// `TreeMoves` concept (tree/tree_update.h). It adds the prior and Hastings terms that stattools
// owns, with the same association the update used before the extraction, so the chain stays
// byte-identical.
//
// The parameters' box type is TTree, so TTree is only declared here. The members are defined in
// TStattoolsTreeMoves.cpp, where TTree is complete.
//

#pragma once

#include "Types.h"
#include "stattools/ParametersObservations/TParameter.h"

#include <cstddef>

class TTree;

/// A `TreeMoves` adapter over the alpha and log nu parameters of one tree.
///
/// A parameter fixed with `--<name>.update false` needs no special case. stattools makes its
/// proposal a no-op and its decision a rejection.
class TStattoolsTreeMoves {
public:
	using TypeParamAlpha = stattools::TParameter<SpecAlpha, TTree>;
	using TypeParamLogNu = stattools::TParameter<SpecLogNu, TTree>;

private:
	TypeParamAlpha *_alpha  = nullptr;
	TypeParamLogNu *_log_nu = nullptr;
	size_t _n_cliques       = 0;

public:
	TStattoolsTreeMoves(TypeParamAlpha &alpha, TypeParamLogNu &log_nu, size_t n_cliques)
	    : _alpha(&alpha), _log_nu(&log_nu), _n_cliques(n_cliques) {}

	void propose_alphas();
	[[nodiscard]] double proposed_alpha(size_t clique) const;
	[[nodiscard]] double accepted_alpha(size_t clique) const;
	/// The ratio comes first and the prior ratio is added to it: `(ratio) + prior ratio`.
	bool decide_alpha(size_t clique, double log_ratio);

	[[nodiscard]] double accepted_log_nu(size_t clique) const;
};
