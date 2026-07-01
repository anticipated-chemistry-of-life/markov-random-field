#pragma once

/// Class storing a mass spectrometry run. The iterator goes along the
/// molecule dimension. So for each feature, we have a list of molecule indices and binned
/// likelihoods.
#include "./feature_likelihood.h"
#include "coretools/Containers/TNestedVector.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/devtools.h"
#include <array>
#include <cmath>
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
	size_t feature_a    = 0;
	size_t feature_b    = 0;
	/// log of the proposal (Hastings) ratio q(old|new)/q(new|old) for this move, computed at
	/// proposal time from the counts the proposal sampled from. Zero for the symmetric moves
	/// (`Swap`, `MoveToFree`); non-zero for `ToUnknown`/`FromUnknown`. The Metropolis-Hastings
	/// acceptance must use `log_likelihood_ratio + log_hastings`.
	double log_hastings = 0.0;
	TFeatureLikelihood old_a;
	TFeatureLikelihood new_a;
	TFeatureLikelihood old_b;
	TFeatureLikelihood new_b;
	AssignmentMoveType type = AssignmentMoveType::Invalid;

	[[nodiscard]] bool is_valid() const { return this->type != AssignmentMoveType::Invalid; }
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

	/// The complete vector of molecules that are currently assigned. This will have length
	/// n_molecules (uint8 that will actually be a bool because bools are uint8).
	std::vector<uint8_t> _molecules_assigned;

	size_t _number_of_molecules = 0;

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
	[[nodiscard]] std::vector<TFeatureLikelihood> _free_candidates_of_feature(size_t f) const;
	[[nodiscard]] TAssignmentProposal _propose_to_unknown() const;
	[[nodiscard]] TAssignmentProposal _propose_from_unknown() const;
	[[nodiscard]] TAssignmentProposal _propose_move_to_free() const;
	[[nodiscard]] TAssignmentProposal _propose_swap() const;

	void _reset_molecules_assigned() {
		_molecules_assigned.clear();
		_molecules_assigned.resize(_number_of_molecules, false);
	}

	void _fill_molecules_assigned() {
		_reset_molecules_assigned();
		for (const auto &assignment : _current_assignments) {
			_molecules_assigned.at(assignment.get_molecule_index()) = true;
		}
	}

	/// This will release the memory of that vector else we just store a vector of capacity N
	/// molecules for every single MSMS run which is way to much memory.
	void _drop_molecules_assigned() {
		_molecules_assigned.clear();
		_molecules_assigned.shrink_to_fit();
	}

public:
	TMassSpecRun() = default;
	void initialize(size_t number_of_molecules) {
		_number_of_molecules = number_of_molecules;
		_reset_molecules_assigned();
	}

	[[nodiscard]] size_t capacity() const { return _features.capacity(); }
	void reserve(size_t n) { _features.reserve(n); }
	[[nodiscard]] bool empty() const { return _features.empty(); }
	void add_empty_feature() { _features.push_back(); }

	[[nodiscard]] coretools::TConstView<TFeatureLikelihood>
	get_likelihoods_for_feature(size_t i) const {
		return _features.at(i);
	}
	[[nodiscard]] coretools::TView<TFeatureLikelihood> get_likelihoods_for_feature(size_t i) {
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
	void update_all_assignments() {
		_fill_molecules_assigned();
		const size_t features_size    = _features.size();
		const size_t assignments_size = _current_assignments.size();
		// check that there are as many features than the length of the current assignments
		if (features_size != assignments_size) {
			throw coretools::TDevError(
			    "Number of features does not match the length of the current assignments");
		}
		// do all the updates
		for (size_t i = 0; i < features_size; ++i) {
			auto &assignment       = _current_assignments[i];
			const auto likelihoods = _features.at(i);
		}
		_drop_molecules_assigned();
	}
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
	[[nodiscard]] static inline double probability_of_assignment(bool y, bool ms,
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
};
