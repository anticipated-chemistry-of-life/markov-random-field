#include "./msms_data.h"

TMSMSData::TMSMSData(const std::unique_ptr<TTree> &species_tree,
                     const std::unique_ptr<TTree> &molecules_tree)
    : _species_tree(species_tree), _molecules_tree(molecules_tree) {
	_msms_data.resize(_species_tree->get_number_of_leaves());
}
