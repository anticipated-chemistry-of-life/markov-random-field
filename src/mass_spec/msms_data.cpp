#include "./msms_data.h"
#include "TMarkovField.h"
#include "Types.h"
#include "cli.h"
#include "constants.h"
#include "coretools/Main/TError.h"
#include "storages/storage_backend.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

void TMSMSData::initialize(TDataModel *box) {
	// build molecule leaf names in leaf-index order
	const size_t n = _molecules_tree->get_number_of_leaves();
	std::vector<std::string> molecule_names;
	molecule_names.reserve(n);
	for (size_t i = 0; i < n; ++i) {
		molecule_names.push_back(
		    _molecules_tree->get_node_id(_molecules_tree->get_node_index_from_leaf_index(i)));
	}

	_proba_to_pass_filter->initStorage(
	    box, {n * _number_of_filters},
	    {std::make_shared<coretools::TNamesStrings>(molecule_names)});

	// scalar — default 1-slot names are correct
	_proba_contamination->initStorage(box, {1});
}

void TMSMSData::guess_initial_values() {
	for (size_t i = 0; i < _molecules_tree->get_number_of_leaves(); ++i) {
		_proba_to_pass_filter->set(i, ProgramOptions::INDEX_PROBA_TO_PASS_MS_FILTER);
	}

	_proba_contamination->set(ProgramOptions::PROBA_OF_MS_CONTAMINATION);
}

TMSMSData::TMSMSData(const std::vector<std::unique_ptr<TTree>> &trees,
                     const TMarkovField &markov_field, size_t number_of_filters,
                     TypeParamMassSpecFilter *filter_proba,
                     TypeParamContamination *contamination_proba)
    : _number_of_filters(number_of_filters), _markov_field(markov_field) {
	// Must run before initialize() sizes this parameter's storage: stattools's automatic
	// state-posterior tracking for a byte-domain parameter needs max() <= 254 (see constants.h).
	TypeFilterProbability::setMax(MAX_FILTER_PROBABILITY_INDEX);

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
	_number_of_species      = _species_tree->get_number_of_leaves();
	// once the trees are checked we can initialize the MSMS data
	_dimensions_for_filters = {_number_of_filters, number_of_molecules};
	_msms_data.resize(_number_of_species);
	_species_has_ms_data.resize(_number_of_species, false);

	// TODO: load ms data from file
	// load_ms_data_from_file()

	_proba_to_pass_filter = filter_proba;
	_proba_contamination  = contamination_proba;
	// Registering these two with the DAG, and the DAG edge to the Markov field's stattools
	// parameters, is now TDataModel's job: it owns this object's parameters the same way it owns
	// TLotus's and TSimpleErrorModel's (see TDataModel's constructor).
}

double TMSMSData::ll_ratio_after_filter_move(size_t index) const {
	const auto &multidim      = _get_filter_molecule_pair_multidimensional_index_(index);
	const size_t filter_index = multidim[0];
	const size_t molecule_idx = multidim[1];

	const auto new_index = _proba_to_pass_filter->value(index);
	const auto old_index = _proba_to_pass_filter->oldValue(index);
	const auto new_p_pos = LOG_LINEAR_SPACE_PROBA[new_index];
	const auto new_p_neg = LOG_LINEAR_SPACE_PROBA_REVERSED[new_index];
	const auto old_p_pos = LOG_LINEAR_SPACE_PROBA[old_index];
	const auto old_p_neg = LOG_LINEAR_SPACE_PROBA_REVERSED[old_index];

	// Position 0 should be old, position 1 should be new
	std::array<double, 2> log_lik{};

	// avoid heap allocation by using a local index array
	IndexArray index_in_Y_space    = {0, molecule_idx};
	const TFieldStorage &y_storage = _markov_field.get_Y_matrix();
	for (size_t species_idx = 0; species_idx < _number_of_species; ++species_idx) {
		// if the species has no MS data, skip it
		if (!_species_has_ms_data.at(species_idx)) { continue; }

		// allocate the species index
		index_in_Y_space[0] = species_idx;

		// a missing cell reads as 0, so a direct point lookup covers both cases
		const bool y = y_storage.is_one(index_in_Y_space);
		// filter only enters the likelihood when the molecule is present
		if (!y) { continue; }
		for (const auto &run : get_ms_data_for_species(species_idx)) {
			if (run.filter_index() != filter_index) { continue; }   // not this extraction
			const bool ms = run.is_detected(molecule_idx);
			if (ms) {
				log_lik[1] += new_p_pos;
				log_lik[0] += old_p_pos;
			} else {
				log_lik[1] += new_p_neg;
				log_lik[0] += old_p_neg;
			}
		}
	}
	return log_lik[1] - log_lik[0];
}

double TMSMSData::calculate_LL_ratio_for_assignment_move(size_t species_idx,
                                                         const TMassSpecRun &run,
                                                         const TAssignmentProposal &move) const {
	if (!move.is_valid()) {
		throw coretools::TDevError(
		    "Move is invalid and should have been handeled before calling this method");
	}

	// accumulate probabilities (TSumLogProbability::add takes a probability; getSum returns its
	// log).
	coretools::TSumLogProbability p_new{};
	coretools::TSumLogProbability p_old{};

	// --- feature term: P(feature | assigned molecule); the unknown molecule's binned value lives in
	// the same assignment, so it reads the same way. ---
	p_old.add(mass_spec_math::probability_of_feature(move.old_a));
	p_new.add(mass_spec_math::probability_of_feature(move.new_a));
	if (move.touches_two_features()) {
		p_old.add(mass_spec_math::probability_of_feature(move.old_b));
		p_new.add(mass_spec_math::probability_of_feature(move.new_b));
	}

	// --- filter / contamination term: only real molecules whose "is assigned" status flips matter.
	// Net the real molecules leaving (old assignments) against those arriving (new assignments);
	// any molecule present in both cancels and contributes nothing. Both swap types cancel
	// completely: a plain Swap because each molecule stays assigned, just to the other feature, and
	// a SwapWithUnknown because the molecule that changes hands leaves the other feature and
	// arrives at this one in the same move. ---
	std::array<uint32_t, 2> removed{};
	std::array<uint32_t, 2> added{};
	size_t n_removed = 0;
	size_t n_added   = 0;
	auto push_real   = [](std::array<uint32_t, 2> &out, size_t &n, const TFeatureLikelihood &a) {
		if (!a.is_unknown_molecule()) { out[n++] = a.get_molecule_index(); }
	};
	push_real(removed, n_removed, move.old_a);
	push_real(added, n_added, move.new_a);
	if (move.touches_two_features()) {
		push_real(removed, n_removed, move.old_b);
		push_real(added, n_added, move.new_b);
	}
	auto contains = [](const std::array<uint32_t, 2> &arr, size_t n, uint32_t m) {
		for (size_t k = 0; k < n; ++k) {
			if (arr[k] == m) { return true; }
		}
		return false;
	};

	const double cont              = (double)_proba_contamination->value();
	const TFieldStorage &y_storage = _markov_field.get_Y_matrix();
	auto y_present                 = [&](uint32_t molecule_idx) {
		return y_storage.is_one(
		    y_storage.get_linear_index_in_Y_space({species_idx, (size_t)molecule_idx}));
	};
	auto filter_prob = [&](uint32_t molecule_idx) {
		return _filter_probability(run.filter_index(), molecule_idx);
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

double TMSMSData::_log_likelihood(double contamination) const {
	return mass_spec_math::log_likelihood(
	    _msms_data, _markov_field.get_Y_matrix(),
	    [this](size_t filter, uint32_t molecule) { return _filter_probability(filter, molecule); },
	    contamination);
}

namespace {
/// Supplies `TMassSpecRun` with the part of the assignment likelihood ratio that lives in
/// `TMSMSData`: Y, the mass-spec filter probabilities and the contamination probability. One
/// instance is bound to a single (species, run) pair for the duration of one update.
class TMSMSDataScorer final : public TAssignmentScorer {
private:
	const TMSMSData &_data;
	size_t _species_idx;
	const TMassSpecRun &_run;

public:
	TMSMSDataScorer(const TMSMSData &data, size_t species_idx, const TMassSpecRun &run)
	    : _data(data), _species_idx(species_idx), _run(run) {}

	[[nodiscard]] double log_likelihood_ratio(const TAssignmentProposal &proposal) const override {
		return _data.calculate_LL_ratio_for_assignment_move(_species_idx, _run, proposal);
	}
};
} // namespace

void TMSMSData::update_all_MS_assignments() {
	_check_trees_exist();
	// Each run is swept by exactly one thread and everything the scorer reads (Y, the filter and
	// contamination parameters) is read-only for the duration. Note that the random generator is
	// thread-local, so results are not reproducible across different numbers of threads.
#pragma omp parallel default(none) num_threads(ProgramOptions::NUMBER_OF_THREADS)
	{
		for (size_t i = 0; i < _number_of_species; ++i) {
			if (!_species_has_ms_data.at(i)) continue;
			auto ms_runs = _msms_data.at(i);
#pragma omp for nowait
			for (auto &run : ms_runs) {
				const TMSMSDataScorer scorer(*this, i, run);
				run.update_all_assignments(scorer, ProgramOptions::MS_PROBA_MOVE_TO_UNKNOWN);
			}
		}
	}
}
