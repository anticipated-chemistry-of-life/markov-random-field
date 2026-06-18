#include "./msms_data.h"
#include "TMarkovField.h"
#include "Types.h"
#include "cli.h"
#include "constants.h"
#include "coretools/Main/TError.h"
#include "storages/y_storage/TStorageYMatrix.h"
#include <array>
#include <cstddef>

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

	// initialize _curLL
	_cur_LL_contamination = 0.0;
	_old_LL_contamination = _cur_LL_contamination;
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
	_old_LL_contamination = 0.0;
	_cur_LL_contamination = 0.0;
}

void TMSMSData::add_to_sumlog(coretools::TSumLogProbability &sum_log, uint8_t binned_value) const {
	sum_log.add(_log_lik_present[binned_value]);
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
	auto feature_prob = [this](const TFeatureLikelihood &assignment) {
		return _log_lik_present[assignment.get_binned_likelihood()];
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

// void TMSMSData::add_to_sumlog(std::array<coretools::TSumLogProbability, 2> &sum_log,
//                               uint8_t binned_value) const {
// 	// sum_log[0].add(_log_lik_absent[binned_value]);
// 	// sum_log[1].add(_log_lik_present[binned_value]);
// 	throw coretools::TDevError("We still don't know how to add this.");
// }
