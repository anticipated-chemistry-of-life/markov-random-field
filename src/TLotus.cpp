#include "TLotus.h"
#include "TStorageY.h"
#include "TStorageYVector.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Math/TSumLog.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

void TLotus::load_from_file(const std::string &filename) {
	coretools::instances::logfile().listFlush("Reading links from file '", filename, "' ...");
	coretools::TInputFile file(filename, coretools::FileType::Header);
	const auto &header = file.header();

	if (file.numCols() != 2) {
		UERROR("File '", filename, "' is expected to have 2 columns, but has ", file.numCols(), " !");
	}

	if (header[0] != "species") {
		UERROR("File '", filename, "' is expected to have a header with 'species' as first column, but has ", header[0],
		       " !");
	}

	if (header[1] != "molecules") {
		UERROR("File '", filename, "' is expected to have a header with 'molecules' as second column, but has ",
		       header[1], " !");
	}

	_species_counter.resize(_trees[0].get_number_of_leaves(), 0);
	_molecules_counter.resize(_trees[1].get_number_of_leaves(), 0);

	for (; !file.empty(); file.popFront()) {

		std::string species  = std::string(file.get(0));
		std::string molecule = std::string(file.get(1));

		if (!_trees[0].get_node(_trees[0].get_node_index(species)).isLeaf()) {
			UERROR("Node '", species, "' is not a leaf !");
		}

		if (!_trees[1].get_node(_trees[1].get_node_index(molecule)).isLeaf()) {
			UERROR("Node '", molecule, "' is not a leaf !");
		}

		size_t species_index  = this->_trees[0].get_index_within_leaves(species);
		size_t molecule_index = this->_trees[1].get_index_within_leaves(molecule);

		size_t linear_index_in_Y_space =
		    this->_L_sm.get_linear_index_in_container_space({species_index, molecule_index});
		this->_L_sm.insert_one(linear_index_in_Y_space);

		++_species_counter[species_index];
		++_molecules_counter[molecule_index];
	}
};

double TLotus::calculate_research_effort(size_t linear_index_in_L_space) const {
	auto coordinate = _L_sm.get_multi_dimensional_index(linear_index_in_L_space);
	auto Q_s        = static_cast<double>(_species_counter[coordinate[0]]);
	auto P_m        = static_cast<double>(_molecules_counter[coordinate[1]]);
	return (1 - exp(-0.1 * P_m)) * (1 - exp(-0.1 * Q_s));
};

/// Let x(m, s) denote the variable in X for NP m and species s, which, in case X contains additional dimensions, is
/// obtained by collapsing:
/// x(m, s) is the minimum between 1 and the sum of x(s, m, y) for all y in Y. Meaning that if there is at least one
/// y in Y such that x(s, m, y) is 1, then x(m, s) will be 1. Otherwise, it will be 0. If we have four dimensions,
/// then we will have to collapse the last two dimensions.
void TLotus::_initialize_x_sm(const TStorageYVector &Y) {

	// if we only have two dimensions, then the provided vector Y will be the same as x_ms
	if (_trees.size() == 2) {
		_x_sm = Y;
		return;
	}

	for (auto y = 0; Y.size(); ++y) {
		auto linear_index = Y[y].get_linear_index_in_container_space();
		// we now get the multidimensional index in the Y space
		auto coordinate   = Y.get_multi_dimensional_index(linear_index);

		if (Y[y].is_one()) {
			_x_sm.insert_one(_x_sm.get_linear_index_in_container_space({coordinate[0], coordinate[1]}));
		}
	}
};

double TLotus::calculate_log_likelihood_of_L() const {
	size_t index_in_x = 0;
	size_t index_in_L = 0;

	const auto length_x = _x_sm.size();
	const auto length_L = _L_sm.size();

	coretools::TSumLogProbability sum_log;

	for (size_t i = 0; i < std::max(length_x, length_L); ++i) {

		// if we reach the end of vector L but not the end of vector x, that means that all the other Ls are 0.
		// But x can still be 1 or 0.
		if (index_in_L == length_L) {

			sum_log.add(_calculate_probability_of_L_sm(_x_sm[index_in_x].is_one(), false,
			                                           _x_sm[index_in_x].get_linear_index_in_container_space()));
			++index_in_x;
		} else if (index_in_x == length_x) {
			// if we reach the end of vector x but not the end of vector L, that means that all the other x are 0.
			sum_log.add(_calculate_probability_of_L_sm(false, _L_sm[index_in_L].is_one(),
			                                           _L_sm[index_in_L].get_linear_index_in_container_space()));
			++index_in_L;
		} else if (_x_sm[index_in_x].get_linear_index_in_container_space() ==
		           _L_sm[index_in_L].get_linear_index_in_container_space()) {
			sum_log.add(_calculate_probability_of_L_sm(_x_sm[index_in_x].is_one(), _L_sm[index_in_L].is_one(),
			                                           _L_sm[index_in_L].get_linear_index_in_container_space()));
			++index_in_x;
			++index_in_L;
		} else if (_x_sm[index_in_x].get_linear_index_in_container_space() <
		           _L_sm[index_in_L].get_linear_index_in_container_space()) {
			sum_log.add(_calculate_probability_of_L_sm(_x_sm[index_in_x].is_one(), false,
			                                           _x_sm[index_in_x].get_linear_index_in_container_space()));
			++index_in_x;
		} else {
			sum_log.add(_calculate_probability_of_L_sm(false, _L_sm[index_in_L].is_one(),
			                                           _L_sm[index_in_L].get_linear_index_in_container_space()));
			++index_in_L;
		}
	}
	return sum_log.getSum();
}

double TLotus::_calculate_probability_of_L_sm(bool x_sm, bool L_sm, size_t linear_index_in_L_space) const {
	if (x_sm && L_sm) { return calculate_research_effort(linear_index_in_L_space); }
	if (x_sm && !L_sm) { return 1 - calculate_research_effort(linear_index_in_L_space); }
	if (!x_sm && L_sm) { return 0.0; }
	return 1.0;
}

/// TODO : make sure x_sm is updated when we update Y !
