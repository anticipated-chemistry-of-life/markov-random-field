#include "../TTree.h"
#include "cli.h"

void TTree::_initialize_grid_branch_lengths() {
	// read a, b and K from command-line
	_number_of_bins = ProgramOptions::BRANCH_LENGTHS_BINS;

	const size_t max_type =
	    std::numeric_limits<coretools::underlyingType<TypeBinnedBranchLengths>::type>::max();
	if (_number_of_bins >= max_type) {
		throw coretools::TUserError("More bins (", _number_of_bins, ") required than type allows (",
		                            max_type, ")! Please decrease n_bins or change type of bins.");
	}
	TypeBinnedBranchLengths::setMax(_number_of_bins - 1);

	// calculate Delta
	_delta = 2.0 / ((double)_number_of_bins + 1.0);

	_grid_branch_lengths.resize(_number_of_bins);
	for (size_t k = 0; k < _number_of_bins; ++k) { _grid_branch_lengths[k] = _delta * (double)k; }
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
	for (size_t p = 0; p < pairs.length(); ++p) {
		auto [p1, p2] = pairs.getIndexPair(p);

		const bool both_at_min =
		    _binned_branch_lengths->value(p1) == TypeBinnedBranchLengths::min() &&
		    _binned_branch_lengths->value(p2) == TypeBinnedBranchLengths::min();
		const bool both_at_max =
		    _binned_branch_lengths->value(p1) == TypeBinnedBranchLengths::max() &&
		    _binned_branch_lengths->value(p2) == TypeBinnedBranchLengths::max();
		if (both_at_min || both_at_max) {
			// can not do update, both are at boundary
			_propose_new_branch_lengths(p1, p2, 0);
		} else if (_binned_branch_lengths->value(p1) == TypeBinnedBranchLengths::min()) {
			// p1 is at minimum (at the left-most position) -> can only go to the right
			_propose_new_branch_lengths(p1, p2, 1);
		} else if (_binned_branch_lengths->value(p2) == TypeBinnedBranchLengths::min()) {
			// p2 is at minimum (at the left-most position) -> can only go to the right
			_propose_new_branch_lengths(p1, p2, -1);
		} else if (_binned_branch_lengths->value(p1) == TypeBinnedBranchLengths::max()) {
			// p1 is at maximum (at the right-most position) -> can only go to the left
			_propose_new_branch_lengths(p1, p2, -1);
		} else if (_binned_branch_lengths->value(p2) == TypeBinnedBranchLengths::max()) {
			// p2 is at maximum (at the right-most position) -> can only go to the left
			_propose_new_branch_lengths(p1, p2, 1);
		} else {
			// we can choose to go left or right
			bool left = coretools::instances::randomGenerator().pickOneOfTwo(coretools::P(0.5));
			if (left) {
				_propose_new_branch_lengths(p1, p2, -1);
			} else {
				_propose_new_branch_lengths(p1, p2, 1);
			}
		}
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
			vals[i] = (_delta * (_binned_branch_lengths->value(i) + 0.5));
		}

		// normalize such that the average branch length is 1
		double average = std::reduce(vals.begin(), vals.end(), 0.0) / (double)vals.size();
		for (double &val : vals) { val = val / average; }

		// translate back to bin
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
	std::vector<size_t> binned_branch_lengths;
	binned_branch_lengths.reserve(get_number_of_nodes() - get_number_of_roots());

	size_t sum_index_branches = 0;
	for (size_t i = 0; i < branch_lengths.size(); ++i) { // loop over all nodes
		if (exclude_root && _nodes[i].is_root()) { continue; }
		// find bin
		auto it = std::lower_bound(_grid_branch_lengths.begin(), _grid_branch_lengths.end(),
		                           branch_lengths[i]);

		if (it == _grid_branch_lengths.end()) {
			// last bin
			binned_branch_lengths.push_back(_grid_branch_lengths.size() - 1);
			sum_index_branches += _grid_branch_lengths.size() - 1;
		} else {
			// take the distance between the lower bin and the value and the higher bin and the
			// value and then we take the one that is closer to the value

			auto it_next = it + 1;

			if (it_next == _grid_branch_lengths.end()) {
				// last bin
				binned_branch_lengths.push_back(std::distance(_grid_branch_lengths.begin(), it));
				sum_index_branches += _grid_branch_lengths.size() - 1;
			} else if (std::abs(branch_lengths[i] - *it) < std::abs(branch_lengths[i] - *it_next)) {
				// take the lower bin
				auto index = std::distance(_grid_branch_lengths.begin(), it);
				binned_branch_lengths.push_back(index);
				sum_index_branches += index;
			} else {
				// take the higher bin
				auto index = std::distance(_grid_branch_lengths.begin(), it_next);
				binned_branch_lengths.push_back(index);
				sum_index_branches += index;
			}
		}
	}

	// if total branch length is smaller than one, we randomly sample some branch lengths and
	// increase them of one until we get a total branch length of one Adjust total_branch_length to
	// exactly 1.0
	size_t number_of_branches = get_number_of_nodes() - get_number_of_roots();
	size_t goal               = (number_of_branches * _number_of_bins) / 2;
	while (sum_index_branches != goal) {
		auto idx = coretools::instances::randomGenerator().getRand<size_t>(
		    0, binned_branch_lengths.size() - 1);
		size_t new_bin;
		if (sum_index_branches < goal) {
			// Increase bin index
			if (binned_branch_lengths[idx] >= _grid_branch_lengths.size() - 1) { continue; }
			new_bin = binned_branch_lengths[idx] + 1;
		} else {
			// Decrease bin index
			if (binned_branch_lengths[idx] == 0) { continue; }
			new_bin = binned_branch_lengths[idx] - 1;
		}
		sum_index_branches += (new_bin - binned_branch_lengths[idx]);
		binned_branch_lengths[idx] = new_bin;
	}

	return binned_branch_lengths;
}

void TTree::_bin_branch_lengths_from_tree(std::vector<double> &branch_lengths) {
	_initialize_grid_branch_lengths();
	// normalize such that they sum to one

	double sum = 0.0;
	int count  = 0;
	for (const double &branch_length : branch_lengths) {
		if (branch_length <= 0.0) { continue; }
		sum += branch_length;
		++count;
	}
	double average = sum / (double)count;

	for (double &branch_length : branch_lengths) {
		if (branch_length <= 0.0) { continue; }
		branch_length = branch_length / average;
	}

	_binned_branch_lengths_from_tree = _bin_branch_lengths(branch_lengths, true);
};
