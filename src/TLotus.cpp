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

TLotus::TLotus(const std::vector<TTree> &trees, TypeParamGamma *gamma) : _trees(trees), _gamma(gamma) {
	this->addPriorParameter(_gamma);

	_oldLL = 0.0;
	_curLL = 0.0;
}

TLotus::TLotus(const std::vector<TTree> &trees) : TLotus(trees, nullptr) {}

[[nodiscard]] std::string TLotus::name() const { return "lotus_likelihood"; }

void TLotus::initialize() {
	// initialize storage
	_gamma->initStorage(this, {_dimensions_to_keep.size()});
}

void TLotus::_get_dimensions_to_collapse(const std::vector<std::string> &header) {
	using namespace coretools::instances;
	// all dimensions that are not present in header will be collapsed
	if (header.size() > _trees.size()) {
		UERROR("Lotus can not have more dimensions than there are trees (", header.size(), " vs ", _trees.size(), ")");
	}

	for (const auto &tree_name : header) {
		bool found = false;
		for (size_t j = 0; j < _trees.size(); ++j) {
			if (_trees[j].get_tree_name() == tree_name) {
				// already found before -> duplicated!
				if (found) { UERROR("Duplicate column name '", tree_name, "' in lotus file!"); }
				// else: remember this dimension -> we will not collapse it
				_dimensions_to_keep.push_back(j);
				_len_per_dimension_lotus.push_back(_trees[j].get_number_of_leaves());
				found = true;
			}
		}
		if (!found) { UERROR("Could not find tree with name '", tree_name, "' in trees (required by lotus)."); }
	}

	if (_dimensions_to_keep.empty()) { UERROR("No dimensions in lotus file are kept!"); }
	if (_dimensions_to_keep.back() != _trees.size() - 1) {
		UERROR("Last dimension of trees and lotus must be identical (", _dimensions_to_keep.back(), " vs ",
		       _trees.size() - 1, ")!");
	}

	// find dimensions to collapse
	for (size_t i = 0; i < _trees.size(); ++i) {
		if (std::find(_dimensions_to_keep.begin(), _dimensions_to_keep.end(), i) == _dimensions_to_keep.end()) {
			// not found in keep -> collapse
			_dimensions_to_collapse.emplace_back(i);
		}
	}

	// report to logfile
	logfile().startIndent("Will keep the following dimensions for lotus: ", _dimensions_to_keep, ":");
	for (const auto i : _dimensions_to_keep) { logfile().list(_trees[i].get_tree_name()); }
	logfile().endIndent();
	if (_dimensions_to_collapse.empty()) {
		logfile().list("Will not collapse lotus.");
	} else {
		logfile().startIndent("Will collapse the following dimensions for lotus: ", _dimensions_to_collapse, ":");
		for (const auto i : _dimensions_to_collapse) { logfile().list(_trees[i].get_tree_name()); }
		logfile().endIndent();
	}

	// initialize the size of L
	_L.initialize(0, _len_per_dimension_lotus);
}

void TLotus::load_from_file(const std::string &filename) {
	coretools::instances::logfile().listFlush("Reading links from file '", filename, "' ...");
	coretools::TInputFile file(filename, coretools::FileType::Header);

	_get_dimensions_to_collapse(file.header());

	_occurrence_counters.resize(_dimensions_to_keep.size());
	for (size_t i = 0; i < _dimensions_to_keep.size(); ++i) {
		_occurrence_counters[i].resize(_trees[_dimensions_to_keep[i]].get_number_of_leaves(), 0);
	}

	std::vector<size_t> index_in_leaves(_dimensions_to_keep.size());
	for (; !file.empty(); file.popFront()) {
		// loop over all columns
		for (size_t i = 0; i < _dimensions_to_keep.size(); ++i) {
			std::string node_name = std::string(file.get(i));

			const size_t tree_index = _dimensions_to_keep[i];
			if (!_trees[tree_index].isLeaf(_trees[tree_index].get_node_index(node_name))) {
				UERROR("Node '", node_name, "' in tree '", _trees[tree_index].get_tree_name(), "' is not a leaf !");
			}

			const size_t ix    = _trees[tree_index].get_index_within_leaves(node_name);
			index_in_leaves[i] = ix;
			++_occurrence_counters[i][ix];
		}

		size_t linear_index_in_Y_space = _L.get_linear_index_in_container_space(index_in_leaves);
		_L.insert_one(linear_index_in_Y_space);
	}
}

double TLotus::calculate_research_effort(size_t linear_index_in_L_space) const {
	auto index_in_L_space = _L.get_multi_dimensional_index(linear_index_in_L_space);

	double prod = 1.0;
	for (size_t i = 0; i < _dimensions_to_keep.size(); ++i) {
		const auto counts = (double)_occurrence_counters[i][index_in_L_space[i]];
		prod *= 1.0 - exp(-_gamma->value(i) * counts);
	}
	return prod;
}

double TLotus::getSumLogPriorDensity(const Storage &) const { return calculate_log_likelihood_of_L(); };

bool TLotus::_x_is_one(const std::vector<size_t> &index_in_Lotus) const {
	// fill fixed indices of clique (e.g. species and molecule)
	std::vector<size_t> index_in_leaves(_trees.size());
	for (size_t i = 0; i < _dimensions_to_keep.size(); ++i) {
		index_in_leaves[_dimensions_to_keep[i]] = index_in_Lotus[i];
	}

	// now loop over all possible starting positions of cliques for the fixed indices
	// TODO: implement
	// for (size_t i = 0; i < _dimensions_to_collapse.size(); ++i) {
	//	_trees[_dimensions_to_collapse[i]].get_clique()
	//}
}

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

double TLotus::calculateLLRatio(TypeParamGamma *, size_t Index, const Storage &) {
	_oldLL = _curLL;                          // store current likelihood
	_curLL = calculate_log_likelihood_of_L(); // calculate likelihood of new gamma
	return _curLL - _oldLL;
}

void TLotus::updateTempVals(TypeParamGamma *, size_t Index, bool Accepted) {
	if (!Accepted) {
		_curLL = _oldLL; // reset
	}
}

void TLotus::set_x(bool state, size_t linear_index_in_L_space) {
	// TODO : make sure x_sm is updated when we update Y !
	// TODO: make sure _curLL is adjusted to match the new x
}

void TLotus::guessInitialValues() {
	// TODO: think of an initialization scheme
	_gamma->set(0.001);
}

void TLotus::_simulateUnderPrior(Storage *) {
	// TODO: implement simulator
}
