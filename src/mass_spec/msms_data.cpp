#include "./msms_data.h"
#include "cli.h"
#include "coretools/Main/TError.h"

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
	    this, {n}, {std::make_shared<coretools::TNamesStrings>(molecule_names)});

	// scalar — default 1-slot names are correct
	_contamination_probability->initStorage(this, {1});
}

void TMSMSData::guessInitialValues() {
	for (size_t i = 0; i < _molecules_tree->get_number_of_leaves(); ++i) {
		_proba_to_pass_filter->set(i,
		                           coretools::Probability(ProgramOptions::PROBA_TO_PASS_MS_FILTER));
	}

	_contamination_probability->set(ProgramOptions::PROBA_OF_MS_CONTAMINATION);

	// initialize _curLL
	_curLL = 0.0;
	_oldLL = _curLL;
}

TMSMSData::TMSMSData(
    const std::vector<std::unique_ptr<TTree>> &trees,
    const std::vector<std::unique_ptr<stattools::TParameter<SpecMarkovField, TLotus>>>
        &markov_field_stattools_param,
    TypeParamMassSpecFilter *filter_proba, TypeParamContamination *contamination_proba)
    : _markov_field_stattools_param(markov_field_stattools_param) {
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

	_msms_data.resize(_species_tree->get_number_of_leaves());

	_proba_to_pass_filter      = filter_proba;
	_contamination_probability = contamination_proba;
	// Register parameters with the DAG — same pattern as TLotus
	this->addPriorParameter({_proba_to_pass_filter, _contamination_probability});
	_oldLL = 0.0;
	_curLL = 0.0;
}

void TMSMSData::add_to_sumlog(coretools::TSumLogProbability &sum_log, uint8_t binned_value) const {
	sum_log.add(_log_lik_present[binned_value]);
}

// void TMSMSData::add_to_sumlog(std::array<coretools::TSumLogProbability, 2> &sum_log,
//                               uint8_t binned_value) const {
// 	// sum_log[0].add(_log_lik_absent[binned_value]);
// 	// sum_log[1].add(_log_lik_present[binned_value]);
// 	throw coretools::TDevError("We still don't know how to add this.");
// }
