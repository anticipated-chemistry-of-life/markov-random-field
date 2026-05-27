#include "./msms_data.h"
#include "coretools/Main/TError.h"

TMSMSData::TMSMSData(const std::vector<std::unique_ptr<TTree>> &trees) {
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
