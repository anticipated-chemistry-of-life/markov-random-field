#include "TLotus.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/devtools.h"
#include <cstddef>
#include <vector>

void TLotus::load_from_file(const std::string &filename) {
	coretools::instances::logfile().listFlush("Reading links from file '", filename, "' ...");
	coretools::TInputFile file(filename, coretools::FileType::Header);
	const auto &header = file.header();
	std::vector<std::string> tree_names(this->_trees.size());
	for (size_t i = 0; i < this->_trees.size(); ++i) { tree_names[i] = this->_trees[i].get_tree_name(); }
	if (header[0] != "species") {
		UERROR("File '", filename, "' is expected to have a header with 'species' as first column, but has ", header[0],
		       " !");
	}

	if (header[1] != "molecules") {
		UERROR("File '", filename, "' is expected to have a header with 'molecules' as second column, but has ",
		       header[1], " !");
	}

	if (file.numCols() != 2) {
		UERROR("File '", filename, "' is expected to have 2 columns, but has ", file.numCols(), " !");
	}

	for (; !file.empty(); file.popFront()) {
		std::vector<size_t> coordinate(2); // since Lotus only has species and molecules, we know that the coordinate
		                                   // will always be of size 2

		for (size_t i = 0; i < 2; ++i) {
			std::string node = std::string(file.get(i));

			size_t tree_index_of_node      = this->_get_tree_index_of_node(node);
			// we get the index within the leaves of the tree
			coordinate[tree_index_of_node] = this->_trees[tree_index_of_node].get_index_within_leaves(node);
		}
		// we currently have the multidimensional index in the Y space. We need to convert it to a linear index in Y
		// space and store it in the data_Y vector.
		size_t linear_index_in_Y_space = this->_data_L.get_linear_index_in_container_space(coordinate);
		this->_data_L.insert_one(linear_index_in_Y_space);

		std::string species  = std::string(file.get(0));
		std::string molecule = std::string(file.get(1));

		// now we add the species and molecule to the counter. If the species is still not in the map, we add it and
		// add 1. If it is already in the map, we increment the counter by 1.
		if (auto it_species = _species_counter.find(species); it_species == _species_counter.end()) {
			_species_counter[species] = 1;
		} else {
			_species_counter[species]++;
		}

		// do the same for the molecules
		if (auto it_molecules = _molecules_counter.find(molecule); it_molecules == _molecules_counter.end()) {
			_molecules_counter[molecule] = 1;
		} else {
			_molecules_counter[molecule]++;
		}
	}
};
