#include "./msms_data.h"
#include "TMarkovField.h"
#include "Types.h"
#include "cli.h"
#include "constants.h"
#include "coretools/Main/TError.h"
#include "storages/y_storage/TStorageYMatrix.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

void TMSMSData::initialize() {
	// build molecule leaf names in leaf-index order
	const size_t n = _molecules_tree->get_number_of_leaves();
	std::vector<std::string> molecule_names;
	molecule_names.reserve(n);
	for (size_t i = 0; i < n; ++i) {
		molecule_names.push_back(
		    _molecules_tree->get_node_id(_molecules_tree->get_node_index_from_leaf_index(i)));
	}

	_proba_to_pass_filter->initStorage(
	    this, {n * _number_of_filters},
	    {std::make_shared<coretools::TNamesStrings>(molecule_names)});

	// scalar — default 1-slot names are correct
	_proba_contamination->initStorage(this, {1});
}

void TMSMSData::guessInitialValues() {
	for (size_t i = 0; i < _molecules_tree->get_number_of_leaves(); ++i) {
		_proba_to_pass_filter->set(i, ProgramOptions::INDEX_PROBA_TO_PASS_MS_FILTER);
	}

	_proba_contamination->set(ProgramOptions::PROBA_OF_MS_CONTAMINATION);
}

TMSMSData::TMSMSData(
    const std::vector<std::unique_ptr<TTree>> &trees, const TMarkovField &markov_field,
    size_t number_of_filters,
    const std::vector<std::unique_ptr<stattools::TParameter<SpecMarkovField, TLotus>>>
        &markov_field_stattools_param,
    TypeParamMassSpecFilter *filter_proba, TypeParamContamination *contamination_proba)
    : _number_of_filters(number_of_filters), _markov_field(markov_field),
      _markov_field_stattools_param(markov_field_stattools_param) {
	for (size_t d = 0; d < trees.size(); ++d) {
		const auto &tree = trees[d];
		if (tree->get_tree_name() == "species") {
			_species_tree = tree.get();
			_species_dim  = d;
		} else if (tree->get_tree_name() == "molecules") {
			_molecules_tree = tree.get();
			_molecule_dim   = d;
		}
	}
	this->_check_trees_exist();
	const size_t number_of_molecules = _molecules_tree->get_number_of_leaves();
	if (number_of_molecules >= MAX_NUMBER_OF_MOLECULES - 1) {
		throw coretools::TUserError("Number of molecules exceeds the maximum allowed (" +
		                            std::to_string(MAX_NUMBER_OF_MOLECULES) + ")");
	}
	// once the trees are checked we can initialize the MSMS data
	_dimensions_for_filters = {_number_of_filters, number_of_molecules};
	_msms_data.resize(_species_tree->get_number_of_leaves());

	_proba_to_pass_filter = filter_proba;
	_proba_contamination  = contamination_proba;
	// Register parameters with the DAG — same pattern as TLotus
	this->addPriorParameter({_proba_to_pass_filter, _proba_contamination});
}

double TMSMSData::calculateLLRatio(TypeParamMassSpecFilter *, size_t index) {
	const auto multidim       = _get_filter_molecule_pair_multidimensional_index_(index);
	const size_t filter_index = multidim[0];
	const size_t molecule_idx = multidim[1];

	const auto new_p  = LINEAR_SPACE_PROBA[_proba_to_pass_filter->value(index)];
	const auto old_p  = LINEAR_SPACE_PROBA[_proba_to_pass_filter->oldValue(index)];
	const double cont = (double)_proba_contamination->value();

	// Position 0 should be old, position 1 should be new
	std::array<coretools::TSumLogProbability, 2> log_lik{};
	const TStorageYMatrix &y_storage = _markov_field.get_Y_matrix();
	for (size_t species_idx = 0; species_idx < _species_tree->get_number_of_leaves();
	     ++species_idx) {
		const auto linear_index =
		    y_storage.get_linear_index_in_Y_space({species_idx, molecule_idx});
		// a missing cell reads as 0, so a direct point lookup covers both cases
		const bool y = y_storage.is_one(linear_index);
		// filter only enters the likelihood when the molecule is present
		if (!y) { continue; }
		for (const auto &run : get_ms_data_for_species(species_idx)) {
			if (run.size() == 0) { continue; }
			if (run.filter_index() != filter_index) { continue; }        // not this extraction
			const bool ms      = run.is_molecule_assigned(molecule_idx); // hook (2)
			const double p_new = TMassSpecRun::probability_of_assignment(y, ms, new_p, cont);
			log_lik[1].add(p_new);
			const double p_old = TMassSpecRun::probability_of_assignment(y, ms, old_p, cont);
			log_lik[0].add(p_old);
		}
	}
	return log_lik[1].getSum() - log_lik[0].getSum();
}

double TMSMSData::calculate_LL_ratio_for_assignment_move(size_t species_idx,
                                                         const TMassSpecRun &run,
                                                         const TAssignmentProposal &move) const {
	if (!move.is_valid()) { return 0.0; }

	// accumulate probabilities (TSumLogProbability::add takes a probability; getSum returns its
	// log).
	coretools::TSumLogProbability p_new{};
	coretools::TSumLogProbability p_old{};

	// --- feature term: P(feature | assigned molecule) = assigned binned likelihood through the
	// _log_lik_present table; the unknown molecule's binned value lives in the same assignment, so
	// it flows through the same table. ---
	auto feature_prob = [](const TFeatureLikelihood &assignment) {
		return PHRED_LIKE_PROBA[assignment.get_binned_likelihood()];
	};
	p_old.add(feature_prob(move.old_a));
	p_new.add(feature_prob(move.new_a));
	if (move.type == AssignmentMoveType::Swap) {
		p_old.add(feature_prob(move.old_b));
		p_new.add(feature_prob(move.new_b));
	}

	// --- filter / contamination term: only real molecules whose "is assigned" status flips matter.
	// Net the real molecules leaving (old assignments) against those arriving (new assignments);
	// any molecule present in both (a Swap) cancels and contributes nothing. ---
	std::array<uint32_t, 2> removed{};
	std::array<uint32_t, 2> added{};
	size_t n_removed = 0;
	size_t n_added   = 0;
	auto push_real   = [](std::array<uint32_t, 2> &out, size_t &n, const TFeatureLikelihood &a) {
		if (!a.is_unknown_molecule()) { out[n++] = a.get_molecule_index(); }
	};
	push_real(removed, n_removed, move.old_a);
	push_real(added, n_added, move.new_a);
	if (move.type == AssignmentMoveType::Swap) {
		push_real(removed, n_removed, move.old_b);
		push_real(added, n_added, move.new_b);
	}
	auto contains = [](const std::array<uint32_t, 2> &arr, size_t n, uint32_t m) {
		for (size_t k = 0; k < n; ++k) {
			if (arr[k] == m) { return true; }
		}
		return false;
	};

	const double cont                = (double)_proba_contamination->value();
	const TStorageYMatrix &y_storage = _markov_field.get_Y_matrix();
	auto y_present                   = [&](uint32_t molecule_idx) {
		return y_storage.is_one(
		    y_storage.get_linear_index_in_Y_space({species_idx, (size_t)molecule_idx}));
	};
	auto filter_prob = [&](uint32_t molecule_idx) {
		const size_t idx =
		    _get_linear_index_filter_molecule_pair({run.filter_index(), molecule_idx});
		return LINEAR_SPACE_PROBA[_proba_to_pass_filter->value(idx)];
	};

	// molecules that were assigned and now are not: ms 1 -> 0
	for (size_t k = 0; k < n_removed; ++k) {
		const uint32_t m = removed[k];
		if (contains(added, n_added, m)) { continue; }
		const bool y      = y_present(m);
		const double filt = filter_prob(m);
		p_old.add(TMassSpecRun::probability_of_assignment(y, true, filt, cont));
		p_new.add(TMassSpecRun::probability_of_assignment(y, false, filt, cont));
	}
	// molecules that were not assigned and now are: ms 0 -> 1
	for (size_t k = 0; k < n_added; ++k) {
		const uint32_t m = added[k];
		if (contains(removed, n_removed, m)) { continue; }
		const bool y      = y_present(m);
		const double filt = filter_prob(m);
		p_old.add(TMassSpecRun::probability_of_assignment(y, false, filt, cont));
		p_new.add(TMassSpecRun::probability_of_assignment(y, true, filt, cont));
	}

	return p_new.getSum() - p_old.getSum();
}

double TMSMSData::_calculate_log_likelihood_of_MSData(double contamination) const {
	const double cont             = contamination;
	const double log_1_minus_cont = std::log(1.0 - cont);
	const size_t n_molecules      = _molecules_tree->get_number_of_leaves();
	const TStorageYMatrix &Y      = _markov_field.get_Y_matrix();

	// Present (Y == 1) molecules grouped by species. Y is sparse, so this is O(nNonZero(Y)); the
	// streaming cursor walks the stored cells without materializing a vector.
	std::vector<std::vector<uint32_t>> present_per_species(_species_tree->get_number_of_leaves());
	for (auto cursor = Y.stored_cursor(); cursor.valid(); cursor.advance()) {
		if (!cursor.is_one()) { continue; }
		const auto md = Y.get_multi_dimensional_index(cursor.linear_index()); // {species, molecule}
		present_per_species[md[0]].push_back((uint32_t)md[1]);
	}

	coretools::TSumLogProbability sum_log{};
	double bulk_log = 0.0; // (1 - contamination) terms for absent & unassigned molecules

	for (size_t species_idx = 0; species_idx < _species_tree->get_number_of_leaves();
	     ++species_idx) {
		const auto runs = get_ms_data_for_species(species_idx);
		if (runs.empty()) { continue; }
		const auto &present = present_per_species[species_idx];

		for (const auto &run : runs) {
			if (run.size() == 0) { continue; }
			const size_t filter_index = run.filter_index();
			auto filter_prob          = [&](uint32_t m) {
				return LINEAR_SPACE_PROBA[_proba_to_pass_filter->value(
				    _get_linear_index_filter_molecule_pair({filter_index, m}))];
			};

			// feature term P(feature | assigned molecule); collect the real molecules currently
			// assigned in this run (ms == 1).
			std::vector<uint32_t> assigned;
			assigned.reserve(run.number_of_assignments());
			for (size_t f = 0; f < run.number_of_assignments(); ++f) {
				const auto &a = run.get_current_assignment(f);
				sum_log.add(PHRED_LIKE_PROBA[a.get_binned_likelihood()]);
				if (!a.is_unknown_molecule()) { assigned.push_back(a.get_molecule_index()); }
			}
			std::sort(assigned.begin(), assigned.end());

			// detection term. The default case (Y = 0, ms = 0) -> (1 - cont) is deferred to the
			// bulk term; only present (Y = 1) and assigned (ms = 1) molecules deviate from it.
			size_t n_special = 0;
			for (const uint32_t m : present) { // Y = 1
				const bool ms = std::binary_search(assigned.begin(), assigned.end(), m);
				sum_log.add(
				    TMassSpecRun::probability_of_assignment(true, ms, filter_prob(m), cont));
				++n_special;
			}
			for (const uint32_t m : assigned) { // ms = 1; skip those already scored as present
				if (Y.is_one(Y.get_linear_index_in_Y_space({species_idx, (size_t)m}))) { continue; }
				sum_log.add(
				    TMassSpecRun::probability_of_assignment(false, true, filter_prob(m), cont));
				++n_special;
			}
			bulk_log += (double)(n_molecules - n_special) * log_1_minus_cont;
		}
	}
	return sum_log.getSum() + bulk_log;
}

void TMSMSData::update_all_MS_assignments() {
	const size_t number_of_species = this->_species_tree->get_number_of_leaves();
	for (size_t i = 0; i < number_of_species; ++i) {
		if (this->number_of_ms_runs_for_species(i) == 0) continue;
		auto ms_runs = this->get_ms_data_for_species(i);
		for (auto &run : ms_runs) {
			const auto move = run.propose_move();
			if (!move.is_valid()) continue;
			const auto ll_ratio = this->calculate_LL_ratio_for_assignment_move(i, run, move);
			// full log odds ratio = log-likelihood ratio + log proposal (Hastings) ratio
			const bool accepted = coretools::TAcceptOddsRatio::accept(ll_ratio + move.log_hastings);
			if (accepted) run.apply_move(move);
		}
	}
}
