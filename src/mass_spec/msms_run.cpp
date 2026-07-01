#include "./msms_run.h"

std::vector<TFeatureLikelihood> TMassSpecRun::_free_candidates_of_feature(size_t f) const {
	std::vector<TFeatureLikelihood> free_candidates;
	for (const auto &candidate : _features.at(f)) {
		const uint32_t molecule_idx = candidate.get_molecule_index();
		if (molecule_idx >= TFeatureLikelihood::get_unknown_molecule_index()) { continue; }
		if (!is_molecule_assigned(molecule_idx)) { free_candidates.push_back(candidate); }
	}
	return free_candidates;
}

TAssignmentProposal TMassSpecRun::_propose_to_unknown() const {
	const auto real = _features_assigned_to_real_molecule();
	if (real.empty()) { return {}; }
	const size_t f = real[coretools::instances::randomGenerator().sample(real.size())];
	TAssignmentProposal proposal;
	proposal.type      = AssignmentMoveType::ToUnknown;
	proposal.feature_a = f;
	proposal.old_a     = _current_assignments[f];
	proposal.new_a     = TFeatureLikelihood::new_unknown_molecule(_probabilities_of_unkowns.at(f));
	// Hastings ratio q(reverse)/q(forward). Forward picks f among R real features (the unknown
	// target is forced) -> q_fwd = 1/R. The reverse is FromUnknown picking f among the (U+1)
	// unknown features and then its old molecule among the (free(f)+1) free candidates it now
	// has (the molecule freed by this move re-enters the set) -> q_rev = 1/((U+1)(free(f)+1)).
	const auto R       = (double)real.size();
	const auto U       = (double)(_current_assignments.size() - real.size());
	const auto free_after = (double)(_free_candidates_of_feature(f).size() + 1);
	proposal.log_hastings = std::log(R / ((U + 1.0) * free_after));
	return proposal;
}

TAssignmentProposal TMassSpecRun::_propose_from_unknown() const {
	const auto unknown = _features_assigned_to_unknown();
	if (unknown.empty()) { return {}; }
	const size_t f = unknown[coretools::instances::randomGenerator().sample(unknown.size())];
	const auto free_candidates = _free_candidates_of_feature(f);
	if (free_candidates.empty()) { return {}; }
	TAssignmentProposal proposal;
	proposal.type      = AssignmentMoveType::FromUnknown;
	proposal.feature_a = f;
	proposal.old_a     = _current_assignments[f];
	proposal.new_a =
	    free_candidates[coretools::instances::randomGenerator().sample(free_candidates.size())];
	// Hastings ratio q(reverse)/q(forward), the exact inverse of the ToUnknown case. Forward
	// picks f among U unknown features and a molecule among its free(f) free candidates ->
	// q_fwd = 1/(U*free(f)). The reverse is ToUnknown picking f among the (R+1) real features
	// (target forced) -> q_rev = 1/(R+1).
	const auto U           = (double)unknown.size();
	const auto R           = (double)(_current_assignments.size() - unknown.size());
	const auto free_before = (double)free_candidates.size();
	proposal.log_hastings  = std::log((U * free_before) / (R + 1.0));
	return proposal;
}

TAssignmentProposal TMassSpecRun::_propose_move_to_free() const {
	const auto real = _features_assigned_to_real_molecule();
	if (real.empty()) { return {}; }
	const size_t f             = real[coretools::instances::randomGenerator().sample(real.size())];
	// the feature's own (assigned) molecule is excluded automatically: it is assigned, hence
	// not a free candidate.
	const auto free_candidates = _free_candidates_of_feature(f);
	if (free_candidates.empty()) { return {}; }
	TAssignmentProposal proposal;
	proposal.type      = AssignmentMoveType::MoveToFree;
	proposal.feature_a = f;
	proposal.old_a     = _current_assignments[f];
	proposal.new_a =
	    free_candidates[coretools::instances::randomGenerator().sample(free_candidates.size())];
	return proposal;
}

TAssignmentProposal TMassSpecRun::_propose_swap() const {
	const auto real = _features_assigned_to_real_molecule();
	if (real.size() < 2) { return {}; }
	auto &rng           = coretools::instances::randomGenerator();
	const size_t pick_i = rng.sample(real.size());
	size_t pick_j       = rng.sample(real.size() - 1);
	if (pick_j >= pick_i) { ++pick_j; } // pick a distinct second feature
	const size_t i = real[pick_i];
	const size_t j = real[pick_j];

	const uint32_t molecule_i = _current_assignments[i].get_molecule_index();
	const uint32_t molecule_j = _current_assignments[j].get_molecule_index();
	// the swap is only legal if each molecule is a candidate of the feature it moves onto.
	const auto i_takes_j      = is_molecule_in_feature(i, molecule_j);
	const auto j_takes_i      = is_molecule_in_feature(j, molecule_i);
	if (!i_takes_j.found || !j_takes_i.found) { return {}; }

	TAssignmentProposal proposal;
	proposal.type      = AssignmentMoveType::Swap;
	proposal.feature_a = i;
	proposal.old_a     = _current_assignments[i];
	proposal.new_a     = TFeatureLikelihood(molecule_j, *i_takes_j.binned_likelihood);
	proposal.feature_b = j;
	proposal.old_b     = _current_assignments[j];
	proposal.new_b     = TFeatureLikelihood(molecule_i, *j_takes_i.binned_likelihood);
	return proposal;
}
