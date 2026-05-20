#include "./msms_data.h"
#include "coretools/Main/TError.h"

TMSMSData::TMSMSData(const std::vector<std::unique_ptr<TTree>> &trees) {
	for (const auto &tree : trees) {
		if (tree->get_tree_name() == "species") {
			_species_tree = tree.get();
		} else if (tree->get_tree_name() == "molecules") {
			_molecules_tree = tree.get();
		}
	}
	if (!_species_tree || !_molecules_tree) {
		throw coretools::TDevError("TMSMSData : missing required tree(s) ");
	}

	_msms_data.resize(_species_tree->get_number_of_leaves());
}
