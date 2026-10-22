#include "tree/TStattoolsTreeMoves.h"
#include "tree/tree_update.h"

static_assert(tree_update::TreeMoves<TStattoolsTreeMoves>);

void TStattoolsTreeMoves::propose_alphas() { _alpha->propose(coretools::TRange(0, _n_cliques, 1)); }

double TStattoolsTreeMoves::proposed_alpha(size_t clique) const { return _alpha->value(clique); }

double TStattoolsTreeMoves::accepted_alpha(size_t clique) const { return _alpha->value(clique); }

bool TStattoolsTreeMoves::decide_alpha(size_t clique, double log_ratio) {
	const double log_H = log_ratio + _alpha->getLogDensityRatio(clique);
	return _alpha->acceptOrReject(log_H, coretools::TRange(clique));
}

double TStattoolsTreeMoves::accepted_log_nu(size_t clique) const { return _log_nu->value(clique); }

void TStattoolsTreeMoves::propose_log_nus() {
	_log_nu->propose(coretools::TRange(0, _n_cliques, 1));
}

double TStattoolsTreeMoves::proposed_log_nu(size_t clique) const { return _log_nu->value(clique); }

bool TStattoolsTreeMoves::decide_log_nu(size_t clique, double log_ratio) {
	const double log_H = log_ratio + _log_nu->getLogDensityRatio(clique);
	return _log_nu->acceptOrReject(log_H, coretools::TRange(clique));
}

void TStattoolsTreeMoves::sample_branch_pairs() { _pairs.sampleIndices(); }

size_t TStattoolsTreeMoves::n_branch_pairs() const { return _pairs.length(); }

std::pair<size_t, size_t> TStattoolsTreeMoves::branch_pair(size_t pair) const {
	return _pairs.getIndexPair(pair);
}

bool TStattoolsTreeMoves::branch_coin() {
	return coretools::instances::randomGenerator().pickOneOfTwo(coretools::P(0.5));
}

size_t TStattoolsTreeMoves::current_bin(size_t branch) const {
	return static_cast<size_t>(_branch_lengths->value(branch));
}

bool TStattoolsTreeMoves::decide_branch_pair(size_t pair, size_t first_bin, size_t second_bin,
                                             double log_ratio) {
	const auto [p1, p2] = _pairs.getIndexPair(pair);
	_branch_lengths->set(p1, static_cast<int>(first_bin));
	_branch_lengths->set(p2, static_cast<int>(second_bin));
	const double log_prior_ratio =
	    _branch_lengths->getLogDensityRatio(p1) + _branch_lengths->getLogDensityRatio(p2);
	const double log_H = log_ratio + log_prior_ratio;
	return _branch_lengths->acceptOrReject(log_H, coretools::TRange(p1, p2), coretools::TRange());
}
