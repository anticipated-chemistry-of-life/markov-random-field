#include "msms_run.h"
#include "msms_data.h" // must be kepts else TMSData is incomplete

double
TMassSpecRun::calculate_probabiliy_from_y_to_assignment(bool y, bool ms,
                                                        const std::vector<size_t> &index_in_y_space,
                                                        size_t index_of_molecule_dimension) const {
	const size_t molecule_index = index_in_y_space[index_of_molecule_dimension];
	if (y & ms) { return _proba_to_pass_filter->value(molecule_index); }
	if (y & !ms) { return 1.0 - _proba_to_pass_filter->value(molecule_index); }
	if (!y & ms) { return _proba_of_contamination->value(); }
	return 1.0 - _proba_of_contamination->value();
}
