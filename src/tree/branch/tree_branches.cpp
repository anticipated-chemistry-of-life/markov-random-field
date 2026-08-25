#include "../TTree.h"
#include "cli.h"

void TTree::_initialize_grid_branch_lengths() {
	// read a, b and K from command-line
	const size_t number_of_bins = ProgramOptions::BRANCH_LENGTHS_BINS;

	const size_t max_type =
	    std::numeric_limits<coretools::underlyingType<TypeBinnedBranchLengths>::type>::max();
	if (number_of_bins >= max_type) {
		throw coretools::TUserError("More bins (", number_of_bins, ") required than type allows (",
		                            max_type, ")! Please decrease n_bins or change type of bins.");
	}
	// The one piece of the grid that cannot live in TBinGrid: this configures a stattools type, so
	// it is global rather than per-grid. TBinGrid carries its own n_bins instead of reading these
	// bounds back out, which is what removes the "initialize the grid first" ordering constraint
	// that used to be implicit in every proposal.
	TypeBinnedBranchLengths::setMax(number_of_bins - 1);

	_bin_grid.emplace(number_of_bins);
}

stattools::TPairIndexSampler TTree::_build_pairs_branch_lengths() const {
	const size_t num_branches = _nodes.size() - get_number_of_roots();
	stattools::TPairIndexSampler sampler(num_branches);
	sampler.sampleIndices();

	return sampler;
}

void TTree::_propose_new_branch_lengths(size_t p1, size_t p2, int val) {
	_binned_branch_lengths->set(p1, (int)_binned_branch_lengths->value(p1) + val);
	_binned_branch_lengths->set(p2, (int)_binned_branch_lengths->value(p2) - val);
}

void TTree::_propose_new_branch_lengths(const stattools::TPairIndexSampler &pairs) {
	const auto &grid = _grid();

	for (size_t p = 0; p < pairs.length(); ++p) {
		auto [p1, p2]       = pairs.getIndexPair(p);
		const size_t first  = static_cast<size_t>(_binned_branch_lengths->value(p1));
		const size_t second = static_cast<size_t>(_binned_branch_lengths->value(p2));

		// Short-circuit on purpose: the coin is drawn only when the direction is actually free. At
		// a boundary the direction is forced, and drawing anyway would consume a random number and
		// shift the whole chain.
		const bool decrease =
		    grid.step_is_free(first, second) &&
		    coretools::instances::randomGenerator().pickOneOfTwo(coretools::P(0.5));

		_propose_new_branch_lengths(p1, p2, grid.step_direction(first, second, decrease));
	}
};

double TTree::_calculate_likelihood_ratio_branch_length(size_t index_in_binned_branch_length,
                                                        const TClique &clique,
                                                        const TCurrentState &current_state) const {
	// translate index in binned branch length vector (of size leaves + internal nodes without
	// roots) to index in nodes
	const size_t index_in_tree =
	    _leaves_and_internal_nodes_without_roots[index_in_binned_branch_length];

	// calculate probability of parent to node for old branch length
	double prob_old = clique.calculate_prob_to_parent<false>(
	    index_in_tree, this, _binned_branch_lengths->oldValue(index_in_binned_branch_length),
	    current_state);

	// calculate probability of parent to node for new branch length
	double prob_new = clique.calculate_prob_to_parent<false>(
	    index_in_tree, this, _binned_branch_lengths->value(index_in_binned_branch_length),
	    current_state);

	return prob_new / prob_old;
}

void TTree::_add_to_LL_branch_lengths(size_t c, const TCurrentState &current_state,
                                      std::vector<coretools::TSumLogProbability> &log_sum,
                                      const stattools::TPairIndexSampler &pairs) const {
	const auto &clique = _cliques[c];

	for (size_t p = 0; p < pairs.length(); ++p) { // loop over all possible pairs
		// get index of branches to calculate LL: p1 and p2
		// that index corresponds to the index in fake, concatenated vector of leaves and internal
		// nodes without roots
		auto [p1, p2]   = pairs.getIndexPair(p);
		double ratio_p1 = _calculate_likelihood_ratio_branch_length(p1, clique, current_state);
		double ratio_p2 = _calculate_likelihood_ratio_branch_length(p2, clique, current_state);

		log_sum[p].add(ratio_p1);
		log_sum[p].add(ratio_p2);
	}
}

void TTree::_evalute_update_branch_length(std::vector<coretools::TSumLogProbability> &log_sum,
                                          const stattools::TPairIndexSampler &pairs) {
	for (size_t p = 0; p < pairs.length(); ++p) {
		const double LL              = log_sum[p].getSum();
		auto [p1, p2]                = pairs.getIndexPair(p);
		const double log_prior_ratio = _binned_branch_lengths->getLogDensityRatio(p1) +
		                               _binned_branch_lengths->getLogDensityRatio(p2);
		const double log_H           = LL + log_prior_ratio;
		_binned_branch_lengths->acceptOrReject(log_H, coretools::TRange(p1, p2),
		                                       coretools::TRange());
	}
}

void TTree::_set_initial_branch_lengths(bool is_simulation) {
	// overwrite simulated branch length: use branch lengths from tree
	if (_binned_branch_lengths->hasFixedInitialValue() || is_simulation) { // use from simulation
		// translate bin into actual branch lengths
		std::vector<double> vals(_binned_branch_lengths->size());
		for (size_t i = 0; i < _binned_branch_lengths->size(); ++i) {
			vals[i] = _grid().grid_branch_length(_binned_branch_lengths->value(i));
		}

		// translate back to bin (_bin_branch_lengths normalizes to mean 1 on the way). These are
		// already branch-space values, so nothing needs excluding.
		auto binned_branch_lengths = _bin_branch_lengths(vals, false);

		// set these values (hack stattools to pretend initial values are not fixed)
		_binned_branch_lengths->fixInitialization(false);
		for (size_t i = 0; i < _binned_branch_lengths->size(); ++i) {
			_binned_branch_lengths->set(i, binned_branch_lengths[i]);

			// we have to do it a second else the oldValue is still the one not normalized.
			// Indeed we first propose the update of the branch length so then in the "update_mu"
			// function, we will have to take the oldValue of the binned branch length.
			_binned_branch_lengths->set(i, binned_branch_lengths[i]);
		}
		_binned_branch_lengths->fixInitialization(true);
	} else { // use from tree
		for (size_t i = 0; i < _binned_branch_lengths->size(); ++i) {
			_binned_branch_lengths->set(i, _binned_branch_lengths_from_tree[i]);
		}
	}
}

std::vector<size_t> TTree::_bin_branch_lengths(const std::vector<double> &branch_lengths,
                                               bool exclude_root) const {
	// Which entries are roots is topology, so the filtering happens here. TBinGrid only ever sees
	// branch-space lengths, all strictly positive -- which is the same set either way, since
	// _add_parent writes 0.0 for a root and the reader rejects every other non-positive length.
	std::vector<double> lengths_of_branches;
	lengths_of_branches.reserve(get_number_of_nodes() - get_number_of_roots());
	for (size_t i = 0; i < branch_lengths.size(); ++i) { // loop over all nodes
		if (exclude_root && _nodes[i].is_root()) { continue; }
		lengths_of_branches.push_back(branch_lengths[i]);
	}

	// Normalizes to mean 1, bins, then walks the result onto the branch-length budget.
	return _grid().bins_from_lengths(lengths_of_branches, [](size_t n) {
		return coretools::instances::randomGenerator().getRand<size_t>(0, n - 1);
	});
}

void TTree::_bin_branch_lengths_from_tree(const std::vector<double> &branch_lengths) {
	_initialize_grid_branch_lengths();
	_binned_branch_lengths_from_tree = _bin_branch_lengths(branch_lengths, true);
};
