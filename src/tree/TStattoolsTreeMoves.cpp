#include "tree/TStattoolsTreeMoves.h"
#include "tree/TTree.h"
#include "tree/tree_update.h"

static_assert(tree_update::TreeMoves<TStattoolsTreeMoves>);

void TStattoolsTreeMoves::propose_alphas() { _alpha->propose(coretools::TRange(0, _n_cliques, 1)); }

double TStattoolsTreeMoves::proposed_alpha(size_t clique) const { return _alpha->value(clique); }

double TStattoolsTreeMoves::accepted_alpha(size_t clique) const { return _alpha->value(clique); }

bool TStattoolsTreeMoves::decide_alpha(size_t clique, double log_ratio) {
	const double log_H = log_ratio + _alpha->getLogDensityRatio(clique);
	return _alpha->acceptOrReject(log_H, coretools::TRange(clique));
}

double TStattoolsTreeMoves::accepted_log_nu(size_t clique) const {
	return _log_nu->value(clique);
}

void TStattoolsTreeMoves::propose_log_nus() {
	_log_nu->propose(coretools::TRange(0, _n_cliques, 1));
}

double TStattoolsTreeMoves::proposed_log_nu(size_t clique) const { return _log_nu->value(clique); }

bool TStattoolsTreeMoves::decide_log_nu(size_t clique, double log_ratio) {
	const double log_H = log_ratio + _log_nu->getLogDensityRatio(clique);
	return _log_nu->acceptOrReject(log_H, coretools::TRange(clique));
}
