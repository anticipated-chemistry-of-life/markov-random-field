#pragma once

/// Class storing a mass spectrometry run. The iterator goes along the
/// molecule dimension. So for each feature, we have a list of molecule indices and binned
/// likelihoods.
#include "./feature_likelihood.h"
#include "coretools/Containers/TNestedVector.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Math/TSumLog.h"
#include "coretools/devtools.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

struct BinarySearchResult {
	bool found{};
	std::optional<uint8_t> binned_likelihood;
	size_t index{};
};

/// The four ways the molecule<->feature assignment of a single run can change in one MCMC step.
/// `Invalid` means no eligible move of the chosen type could be constructed for the current state
/// (the caller should treat it as a no-op / reject).
enum class AssignmentMoveType : uint8_t {
	ToUnknown,   // a feature assigned to a real molecule -> the unknown molecule (frees that
	             // molecule)
	FromUnknown, // a feature assigned to the unknown molecule -> one of its unassigned candidates
	Swap,        // two features exchange their (real) molecules
	MoveToFree,  // a feature assigned to a real molecule -> one of its unassigned candidate
	             // molecules
	Invalid
};

/// A single proposed change to a run's `_current_assignments`. It records, for every feature it
/// touches, the assignment before (`old_*`) and after (`new_*`) the move, so it can be applied,
/// reverted, and scored without re-deriving anything. `feature_b`/`old_b`/`new_b` are only used by
/// `Swap`.
struct TAssignmentProposal {
	AssignmentMoveType type = AssignmentMoveType::Invalid;
	size_t feature_a        = 0;
	TFeatureLikelihood old_a;
	TFeatureLikelihood new_a;
	size_t feature_b = 0;
	TFeatureLikelihood old_b;
	TFeatureLikelihood new_b;

	[[nodiscard]] bool is_valid() const { return type != AssignmentMoveType::Invalid; }
};

class TMassSpecRun {
private:
	coretools::TNestedVector<TFeatureLikelihood> _features;

	/// Each feature also needs a probability for the unknown molecules. In case a feature is
	/// assigned to an unknown molecule, the probability of unknown is used. This vector should be
	/// of size "number of features" in the MSMS run so that each feature can have a different
	/// probability of for the unknown molecule(s).
	std::vector<uint8_t> _probabilities_of_unkowns;

	/// The current assignment of molecules to features. The vector should be of size "number of
	/// features" in the current MSMS run.
	std::vector<TFeatureLikelihood> _current_assignments;

	/// The filter index for the current MSMS run.
	size_t _filter_index = 0;

	// const stattools::TParameter<SpecContaminationProba, TMSMSData> *_proba_of_contamination;
	// const stattools::TParameter<SpecMassSpecFilter, TMSMSData> *_proba_to_pass_filter;

public:
	TMassSpecRun() = default;
	[[nodiscard]] size_t capacity() const { return _features.capacity(); }
	void reserve(size_t n) { _features.reserve(n); }
	[[nodiscard]] bool empty() const { return _features.empty(); }
	void add_empty_feature() { _features.push_back(); }

	[[nodiscard]] coretools::TConstView<TFeatureLikelihood>
	get_likelihoods_for_feature(size_t i) const {
		return _features.at(i);
	}
	[[nodiscard]] size_t size() const { return _features.size(); }
	[[nodiscard]] size_t number_of_features() const { return this->size(); }
	[[nodiscard]] auto begin() const { return _features.begin(); };
	[[nodiscard]] auto end() const { return _features.end(); }
	auto begin() { return _features.begin(); }
	auto end() { return _features.end(); }
	/// Returns the filter index for the current MSMS run.
	[[nodiscard]] size_t filter_index() const { return _filter_index; }
	///
	[[nodiscard]] BinarySearchResult is_molecule_in_feature(size_t feature_idx,
	                                                        uint32_t molecule_index) const {
		const auto &likelihoods = _features.at(feature_idx);
		return is_molecule_in_feature(likelihoods, molecule_index);
	}

	inline static BinarySearchResult
	is_molecule_in_feature(const coretools::TConstView<TFeatureLikelihood> &feature,
	                       uint32_t molecule_index) {
		auto it = std::lower_bound(feature.begin(), feature.end(), molecule_index);
		if (it == feature.end()) { return {false, std::nullopt, feature.size()}; }
		size_t index = std::distance(feature.begin(), it);
		if (it->get_molecule_index() != molecule_index) { return {false, std::nullopt, index}; }
		return {true, it->get_binned_likelihood(), index};
	}

	static double
	get_likelihood_from_binned_value(const std::array<double, 256> &binned_likelihoods,
	                                 uint8_t binned_value) {
		return binned_likelihoods[binned_value];
	}

	/// Adds a vector of molecule likelihoods. This represents one feature with all the molecule
	/// likelihoods that are associated to it. The vector will be sorted before being added to
	/// guarantee the order and be able to do binary search.
	void add_likelihood_vector(std::vector<TFeatureLikelihood> &feature_likelihoods) {
		std::sort(feature_likelihoods.begin(), feature_likelihoods.end());
		_features.push_back(feature_likelihoods);
	}

	[[nodiscard]] bool feature_has_unknown_molecule_assigned(size_t feature_idx) const {
		return _current_assignments.at(feature_idx).get_molecule_index() >=
		       TFeatureLikelihood::get_unknown_molecule_index();
	}

	/// Calculates the probability if a molecule is assigned to a feature based on the
	/// observed MS/MS data and the model's probability estimates.
	/// If the current state in Y for a molecule-species pair is 1 it will check if the
	/// molecule is assigned to the feature and return the corresponding probability.
	///
	/// This value will then be multiplied by the likelihood of the molecule being present
	/// in the feature to get the final probability.
	/// TMassSpecRun — pure: caller passes the filter prob (for this molecule) and contamination
	[[nodiscard]] static double probability_of_assignment(bool y, bool ms,
	                                                      double proba_to_pass_filter,
	                                                      double proba_contamination) {
		if (y && ms) { return proba_to_pass_filter; }        // present & detected -> passed filter
		if (y && !ms) { return 1.0 - proba_to_pass_filter; } // present & not detected
		if (!y && ms) { return proba_contamination; }        // absent  & detected -> contamination
		return 1.0 - proba_contamination;                    // absent  & not detected
	}

	[[nodiscard]] bool is_molecule_assigned(uint32_t molecule_idx) const {
		return std::any_of(_current_assignments.begin(), _current_assignments.end(),
		                   [molecule_idx](const TFeatureLikelihood &a) {
			                   return a.get_molecule_index() == molecule_idx;
		                   });
	}

	//------------------------------------------------------------------------------------------------
	// Assignment latent state: setup
	//------------------------------------------------------------------------------------------------
	void set_filter_index(size_t filter_index) { _filter_index = filter_index; }

	/// One binned probability per feature, used as P(feature | unknown molecule). Must be sized to
	/// `number_of_features()`.
	void set_probabilities_of_unknowns(std::vector<uint8_t> probabilities_of_unknowns) {
		_probabilities_of_unkowns = std::move(probabilities_of_unknowns);
	}
	[[nodiscard]] uint8_t probability_of_unknown(size_t feature_idx) const {
		return _probabilities_of_unkowns.at(feature_idx);
	}

	/// Start every feature off assigned to the unknown molecule (binned prob taken from
	/// `_probabilities_of_unkowns`). Requires `_probabilities_of_unkowns` to already be sized to
	/// the number of features.
	void initialize_assignments_to_unknown() {
		const size_t n = number_of_features();
		if (_probabilities_of_unkowns.size() != n) {
			throw coretools::TDevError("TMassSpecRun: _probabilities_of_unkowns has size ",
			                           _probabilities_of_unkowns.size(), " but the run has ", n,
			                           " features. Set the unknown probabilities first.");
		}
		_current_assignments.assign(n, TFeatureLikelihood{});
		for (size_t f = 0; f < n; ++f) {
			_current_assignments[f] =
			    TFeatureLikelihood::new_unknown_molecule(_probabilities_of_unkowns[f]);
		}
	}

	void set_current_assignment(size_t feature_idx, TFeatureLikelihood assignment) {
		_current_assignments.at(feature_idx) = assignment;
	}
	[[nodiscard]] const TFeatureLikelihood &get_current_assignment(size_t feature_idx) const {
		return _current_assignments.at(feature_idx);
	}
	[[nodiscard]] size_t number_of_assignments() const { return _current_assignments.size(); }

	//------------------------------------------------------------------------------------------------
	// Assignment latent state: MCMC moves
	//
	// A proposal only ever changes one (or, for a swap, two) feature(s) and is built so that the
	// invariant "a real molecule is assigned to at most one feature" is preserved. The unknown
	// molecule is exempt: many features may share it. `apply_move`/`revert_move` are exact
	// inverses, and the proposal carries the old assignments so the likelihood ratio can be scored
	// elsewhere (see TMSMSData::calculate_LL_ratio_for_assignment_move). NOTE: the returned
	// proposals are symmetric only for `Swap`; a correct Metropolis-Hastings step for the other
	// move types still needs the (asymmetric) proposal/Hastings ratio, which is not computed here.
	//------------------------------------------------------------------------------------------------

	/// Pick one of the four move types uniformly at random and try to build a concrete, invariant-
	/// preserving proposal for the current state. Returns an `Invalid` proposal if no such move of
	/// the chosen type exists right now.
	[[nodiscard]] TAssignmentProposal propose_move() const {
		switch (coretools::instances::randomGenerator().sample(4)) {
		case 0: return _propose_to_unknown();
		case 1: return _propose_from_unknown();
		case 2: return _propose_swap();
		default: return _propose_move_to_free();
		}
	}

	void apply_move(const TAssignmentProposal &proposal) {
		_current_assignments.at(proposal.feature_a) = proposal.new_a;
		if (proposal.type == AssignmentMoveType::Swap) {
			_current_assignments.at(proposal.feature_b) = proposal.new_b;
		}
	}

	void revert_move(const TAssignmentProposal &proposal) {
		_current_assignments.at(proposal.feature_a) = proposal.old_a;
		if (proposal.type == AssignmentMoveType::Swap) {
			_current_assignments.at(proposal.feature_b) = proposal.old_b;
		}
	}

private:
	/// Features whose current assignment is a real (non-unknown) molecule.
	[[nodiscard]] std::vector<size_t> _features_assigned_to_real_molecule() const {
		std::vector<size_t> features;
		for (size_t f = 0; f < _current_assignments.size(); ++f) {
			if (!_current_assignments[f].is_unknown_molecule()) { features.push_back(f); }
		}
		return features;
	}

	/// Features currently assigned to the unknown molecule.
	[[nodiscard]] std::vector<size_t> _features_assigned_to_unknown() const {
		std::vector<size_t> features;
		for (size_t f = 0; f < _current_assignments.size(); ++f) {
			if (_current_assignments[f].is_unknown_molecule()) { features.push_back(f); }
		}
		return features;
	}

	/// Candidate molecules of feature `f` (i.e. molecules with a likelihood for it) that are not
	/// currently assigned to any feature -> the molecules `f` could legally move onto.
	[[nodiscard]] std::vector<TFeatureLikelihood> _free_candidates_of_feature(size_t f) const {
		std::vector<TFeatureLikelihood> free_candidates;
		for (const auto &candidate : _features.at(f)) {
			const uint32_t molecule_idx = candidate.get_molecule_index();
			if (molecule_idx >= TFeatureLikelihood::get_unknown_molecule_index()) { continue; }
			if (!is_molecule_assigned(molecule_idx)) { free_candidates.push_back(candidate); }
		}
		return free_candidates;
	}

	[[nodiscard]] TAssignmentProposal _propose_to_unknown() const {
		const auto real = _features_assigned_to_real_molecule();
		if (real.empty()) { return {}; }
		const size_t f = real[coretools::instances::randomGenerator().sample(real.size())];
		TAssignmentProposal proposal;
		proposal.type      = AssignmentMoveType::ToUnknown;
		proposal.feature_a = f;
		proposal.old_a     = _current_assignments[f];
		proposal.new_a = TFeatureLikelihood::new_unknown_molecule(_probabilities_of_unkowns.at(f));
		return proposal;
	}

	[[nodiscard]] TAssignmentProposal _propose_from_unknown() const {
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
		return proposal;
	}

	[[nodiscard]] TAssignmentProposal _propose_move_to_free() const {
		const auto real = _features_assigned_to_real_molecule();
		if (real.empty()) { return {}; }
		const size_t f = real[coretools::instances::randomGenerator().sample(real.size())];
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

	[[nodiscard]] TAssignmentProposal _propose_swap() const {
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
};
