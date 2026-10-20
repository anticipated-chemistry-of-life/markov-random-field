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
		auto binned_branch_lengths = _bin_branch_lengths(vals);

		// set these values (hack stattools to pretend initial values are not fixed)
		_binned_branch_lengths->fixInitialization(false);
		for (size_t i = 0; i < _binned_branch_lengths->size(); ++i) {
			_binned_branch_lengths->set(i, binned_branch_lengths[i]);
		}
		_binned_branch_lengths->fixInitialization(true);
	} else { // use from tree
		for (size_t i = 0; i < _binned_branch_lengths->size(); ++i) {
			_binned_branch_lengths->set(i, _binned_branch_lengths_from_tree[i]);
		}
	}
}

/// Branch-space lengths in, bins out. Roots never appear: a phylogeny holds one length per branch,
/// so there is nothing left to filter -- which is what removes the exclude_root flag this used to
/// need, and with it the question of which index space the caller was in.
std::vector<size_t> TTree::_bin_branch_lengths(const std::vector<double> &branch_lengths) const {
	// Normalizes to mean 1, bins, then walks the result onto the branch-length budget.
	return _grid().bins_from_lengths(branch_lengths, [](size_t n) {
		return coretools::instances::randomGenerator().getRand<size_t>(0, n - 1);
	});
}

void TTree::_bin_branch_lengths_from_tree() {
	_initialize_grid_branch_lengths();
	_binned_branch_lengths_from_tree = _bin_branch_lengths(_topology().branch_lengths());
};
