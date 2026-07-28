#pragma once

/// Class storing a mass spectrometry run. The iterator goes along the
/// molecule dimension. So for each feature, we have a list of molecule indices and binned
/// likelihoods.
#include "./feature_likelihood.h"
#include "coretools/Containers/TNestedVector.h"
#include "coretools/Main/TError.h"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

struct BinarySearchResult {
	bool found{};
	std::optional<uint8_t> binned_likelihood;
	size_t index{};
};

/// The five ways the molecule<->feature assignment of a single run can change in one MCMC step.
/// `Invalid` means no eligible move could be constructed for the current state (the caller should
/// treat it as a no-op / reject).
enum class AssignmentMoveType : uint8_t {
	ToUnknown,   // a feature assigned to a real molecule -> the unknown molecule (frees that
	             // molecule)
	FromUnknown, // a feature assigned to the unknown molecule -> one of its unassigned candidates
	Swap,        // two features exchange their (real) molecules
	MoveToFree,  // a feature assigned to a real molecule -> one of its unassigned candidate
	             // molecules
	SwapWithUnknown, // a feature at the unknown molecule takes a molecule from another feature,
	                 // which is pushed to the unknown molecule in exchange
	Invalid
};

/// A single proposed change to a run's `_current_assignments`. It records, for every feature it
/// touches, the assignment before (`old_*`) and after (`new_*`) the move, so it can be applied,
/// reverted, and scored without re-deriving anything. `feature_b`/`old_b`/`new_b` are only used by
/// the two move types for which `touches_two_features()` is true.
struct TAssignmentProposal {
	size_t feature_a    = 0;
	size_t feature_b    = 0;
	/// log of the proposal (Hastings) ratio q(old|new)/q(new|old) for this move, computed at
	/// proposal time from the candidate counts the proposal sampled from. Zero for the symmetric
	/// moves (`Swap`, `MoveToFree`); non-zero for the others. The Metropolis-Hastings acceptance
	/// must use `log_likelihood_ratio + log_hastings`.
	double log_hastings = 0.0;
	TFeatureLikelihood old_a;
	TFeatureLikelihood new_a;
	TFeatureLikelihood old_b;
	TFeatureLikelihood new_b;
	AssignmentMoveType type = AssignmentMoveType::Invalid;

	[[nodiscard]] bool is_valid() const { return this->type != AssignmentMoveType::Invalid; }

	/// True for the move types that change two features, i.e. the ones for which `feature_b`,
	/// `old_b` and `new_b` are meaningful. Every place that has to decide whether to look at the
	/// `*_b` fields must go through this predicate, so that adding a further two-feature move type
	/// only requires changing it here.
	[[nodiscard]] bool touches_two_features() const {
		return this->type == AssignmentMoveType::Swap ||
		       this->type == AssignmentMoveType::SwapWithUnknown;
	}
};

/// Scores a proposed assignment change of a single run: returns log P(new state) - log P(old state)
/// under the full model. `TMassSpecRun` owns the assignments and the proposal mechanics but knows
/// nothing about Y, the mass-spec filter probabilities or the contamination probability, all of
/// which enter that ratio; `TMSMSData` implements this interface to supply them.
class TAssignmentScorer {
public:
	virtual ~TAssignmentScorer() = default;
	[[nodiscard]] virtual double
	log_likelihood_ratio(const TAssignmentProposal &proposal) const = 0;
};

class TMassSpecRun {
public:
	/// Value stored in the holder map for a molecule that is currently not assigned to any feature.
	static constexpr uint32_t NO_HOLDER = std::numeric_limits<uint32_t>::max();

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

	/// Maps a molecule index to the feature that currently holds it, or `NO_HOLDER` if the molecule
	/// is free. Only alive during `update_all_assignments`, where it turns both "is this molecule
	/// free?" and "which feature holds it?" into O(1) lookups; outside the sweep it is dropped
	/// because keeping one vector of n_molecules entries per MSMS run would be far too much memory.
	/// Every access goes through the five `_*_holder*` methods below, so the representation can be
	/// changed (e.g. to a hash map sized on the number of features) without touching the sweep.
	std::vector<uint32_t> _holder_of_molecule;

	size_t _number_of_molecules = 0;

private:
	void _reset_holders() { _holder_of_molecule.assign(_number_of_molecules, NO_HOLDER); }

	/// Builds the holder map from the current assignments. Features at the unknown molecule are
	/// skipped: the unknown molecule is shared by arbitrarily many features and its sentinel index
	/// is out of range of the map. Throws if a real molecule is held twice, which turns a
	/// bookkeeping bug in one sweep into a loud failure at the start of the next one.
	void _fill_holders() {
		_reset_holders();
		for (size_t f = 0; f < _current_assignments.size(); ++f) {
			const auto &assignment = _current_assignments[f];
			if (assignment.is_unknown_molecule()) { continue; }
			const uint32_t molecule_idx = assignment.get_molecule_index();
			if (_holder_of(molecule_idx) != NO_HOLDER) {
				throw coretools::TDevError("TMassSpecRun: molecule ", molecule_idx,
				                           " is assigned to both feature ",
				                           _holder_of(molecule_idx), " and feature ", f, ".");
			}
			_set_holder(molecule_idx, static_cast<uint32_t>(f));
		}
	}

	/// This will release the memory of that vector else we just store a vector of capacity N
	/// molecules for every single MSMS run which is way to much memory.
	void _drop_holders() {
		_holder_of_molecule.clear();
		_holder_of_molecule.shrink_to_fit();
	}

	[[nodiscard]] uint32_t _holder_of(uint32_t molecule_idx) const {
		return _holder_of_molecule.at(molecule_idx);
	}

	void _set_holder(uint32_t molecule_idx, uint32_t feature_idx) {
		_holder_of_molecule.at(molecule_idx) = feature_idx;
	}

	/// Builds the proposal for a feature that is currently assigned to the unknown molecule.
	[[nodiscard]] TAssignmentProposal _propose_for_unknown_feature(size_t f, double beta) const;
	/// Builds the proposal for a feature that is currently assigned to a real molecule.
	[[nodiscard]] TAssignmentProposal _propose_for_assigned_feature(size_t f, double beta) const;

	[[nodiscard]] TAssignmentProposal _make_from_unknown(size_t f, TFeatureLikelihood target,
	                                                     size_t c_f, double beta) const;
	[[nodiscard]] TAssignmentProposal _make_to_unknown(size_t f, size_t c_f, double beta) const;
	[[nodiscard]] TAssignmentProposal _make_move_to_free(size_t f, TFeatureLikelihood target) const;
	[[nodiscard]] TAssignmentProposal _make_swap(size_t f, size_t other,
	                                             TFeatureLikelihood target) const;
	[[nodiscard]] TAssignmentProposal
	_make_swap_with_unknown(size_t f, size_t other, TFeatureLikelihood target, size_t c_f) const;

	/// Brings the holder map back in sync after `apply_move` accepted `proposal`.
	void _sync_holders_after_accept(const TAssignmentProposal &proposal);

	/// Checks everything `update_all_assignments` relies on and returns false if the run has
	/// nothing to update.
	[[nodiscard]] bool _assignments_are_updatable() const;

public:
	TMassSpecRun() = default;

	/// Size of the molecule space this run's candidate indices live in. Must be set before the
	/// first call to `update_all_assignments`, which needs it to size the holder map.
	void set_number_of_molecules(size_t number_of_molecules) {
		_number_of_molecules = number_of_molecules;
	}
	[[nodiscard]] size_t number_of_molecules() const { return _number_of_molecules; }

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

	/// Number of candidate molecules of feature `i`, i.e. `c_f` in the Hastings ratios.
	[[nodiscard]] size_t number_of_candidates(size_t i) const { return _features.size(i); }

	/// Adds a vector of molecule likelihoods. This represents one feature with all the molecule
	/// likelihoods that are associated to it. The vector will be sorted before being added to
	/// guarantee the order and be able to do binary search.
	void add_likelihood_vector(std::vector<TFeatureLikelihood> &feature_likelihoods);

	/// One Metropolis-Hastings sweep over every feature of this run, see `msms_run.cpp` for the
	/// move types and their Hastings ratios. `beta` is the probability of proposing to move an
	/// already assigned feature back to the unknown molecule.
	void update_all_assignments(const TAssignmentScorer &scorer, double beta);

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
			// bin 0 maps to probability 0 (see PHRED_LIKE_PROBA), which would make the initial
			// state impossible and every likelihood ratio involving it undefined.
			if (_probabilities_of_unkowns[f] == 0) {
				throw coretools::TUserError("TMassSpecRun: the probability of the unknown molecule "
				                            "for feature ",
				                            f, " is 0, which makes that feature impossible.");
			}
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
	// A proposal only ever changes one (or, for the two swap types, two) feature(s) and is built so
	// that the invariant "a real molecule is assigned to at most one feature" is preserved. The
	// unknown molecule is exempt: many features may share it. `apply_move`/`revert_move` are exact
	// inverses, and the proposal carries the old assignments so the likelihood ratio can be scored
	// elsewhere (see TMSMSData::calculate_LL_ratio_for_assignment_move). The proposals are
	// asymmetric in general, so acceptance must use `log_likelihood_ratio + log_hastings`.
	//------------------------------------------------------------------------------------------------

	void apply_move(const TAssignmentProposal &proposal) {
		_current_assignments.at(proposal.feature_a) = proposal.new_a;
		if (proposal.touches_two_features()) {
			_current_assignments.at(proposal.feature_b) = proposal.new_b;
		}
	}

	void revert_move(const TAssignmentProposal &proposal) {
		_current_assignments.at(proposal.feature_a) = proposal.old_a;
		if (proposal.touches_two_features()) {
			_current_assignments.at(proposal.feature_b) = proposal.old_b;
		}
	}
};
