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

TLotus::TLotus(const std::vector<TTree> &trees, TypeParamGamma *gamma)
    : _trees(trees), _collapser(trees), _gamma(gamma), _tmp_state_along_last_dim(trees.back(), 1) {
	this->addPriorParameter(_gamma);

	_oldLL = 0.0;
	_curLL = 0.0;
}

TLotus::TLotus(const std::vector<TTree> &trees) : TLotus(trees, nullptr) {}

[[nodiscard]] std::string TLotus::name() const { return "lotus_likelihood"; }

void TLotus::initialize() {
	// initialize storage
	_gamma->initStorage(this, {_collapser.num_dim_to_keep()});
}

void TLotus::load_from_file(const std::string &filename) {
	coretools::instances::logfile().listFlush("Reading links from file '", filename, "' ...");
	coretools::TInputFile file(filename, coretools::FileType::Header);

	// initialize collapser: know which dimensions to keep and which to collapse
	// TODO: ensure somewhere that the order matches (at least if we don't collapse)
	const auto len_per_dimension_lotus = _collapser.initialize(file.header(), "LOTUS");

	// initialize the size of L
	_L.initialize(0, len_per_dimension_lotus);

	_occurrence_counters.resize(_collapser.num_dim_to_keep()); // for example, size is 2 if keep molecules and species
	for (size_t i = 0; i < _collapser.num_dim_to_keep(); ++i) {
		_occurrence_counters[i].resize(_trees[_collapser.dim_to_keep(i)].get_number_of_leaves(), 0);
	}

	std::vector<size_t> index_in_collapsed_space(_collapser.num_dim_to_keep());
	for (; !file.empty(); file.popFront()) {
		// loop over all columns
		for (size_t i = 0; i < _collapser.num_dim_to_keep(); ++i) {
			std::string node_name = std::string(file.get(i));

			const size_t tree_index = _collapser.dim_to_keep(i);
			if (!_trees[tree_index].isLeaf(_trees[tree_index].get_node_index(node_name))) {
				UERROR("Node '", node_name, "' in tree '", _trees[tree_index].get_tree_name(), "' is not a leaf !");
			}

			const size_t ix             = _trees[tree_index].get_index_within_leaves(node_name);
			index_in_collapsed_space[i] = ix;
			++_occurrence_counters[i][ix];
		}

		size_t linear_index_in_Y_space = _L.get_linear_index_in_container_space(index_in_collapsed_space);
		_L.insert_one(linear_index_in_Y_space);
	}
}

double TLotus::_calculate_research_effort(const std::vector<size_t> &index_in_collapsed_space) const {
	double prod = 1.0;
	for (size_t i = 0; i < _collapser.num_dim_to_keep(); ++i) {
		const size_t leaf_index = index_in_collapsed_space[i];
		const auto counts       = (double)_occurrence_counters[i][leaf_index];
		prod *= 1.0 - exp(-(double)_gamma->value(i) * counts);
	}
	return prod;
}

double TLotus::getSumLogPriorDensity(const Storage &) const { return calculate_log_likelihood_of_L(); };

void TLotus::fill_tmp_state_along_last_dim(const std::vector<size_t> &start_index_clique_along_last_dim, size_t K) {
	// collapse start_index_in_leaves (this is the index in Y)
	if (_collapser.do_collapse()) {
		_tmp_state_along_last_dim.fill_Y_along_last_dim(_collapser.collapse(start_index_clique_along_last_dim), K, _L);
	} else { // no need to collapse
		_tmp_state_along_last_dim.fill_Y_along_last_dim(start_index_clique_along_last_dim, K, _L);
	}
}

void TLotus::calculate_LL_update_Y(const std::vector<size_t> &index_in_leaves_space, bool new_state, bool old_state,
                                   std::array<coretools::TSumLogProbability, 2> &sum_log) {
	const bool x = _collapser.x_is_one(index_in_leaves_space, new_state, old_state);

	const size_t leaf_index_last_dim    = index_in_leaves_space.back();
	const auto index_in_collapsed_space = _collapser.collapse(index_in_leaves_space);
	_calculate_probability_of_L_given_x(x, _tmp_state_along_last_dim.get_Y(leaf_index_last_dim),
	                                    index_in_collapsed_space);
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
			const auto index_in_L_space = _x_sm[index_in_x].get_linear_index_in_container_space());
			sum_log.add(_calculate_probability_of_L_given_x(_x_sm[index_in_x].is_one(), false,
			                                           index_in_L_space);
			++index_in_x;
		} else if (index_in_x == length_x) {
			// if we reach the end of vector x but not the end of vector L, that means that all the other x are 0.
			const auto index_in_L_space = _L_sm[index_in_L].get_linear_index_in_container_space());
			sum_log.add(_calculate_probability_of_L_given_x(false, _L_sm[index_in_L].is_one(),
			                                                index_in_L_space);
			++index_in_L;
		} else if (_x_sm[index_in_x].get_linear_index_in_container_space() ==
		           _L_sm[index_in_L].get_linear_index_in_container_space()) {
			const auto index_in_L_space = _L_sm[index_in_L].get_linear_index_in_container_space());
			sum_log.add(_calculate_probability_of_L_given_x(_x_sm[index_in_x].is_one(), _L_sm[index_in_L].is_one(),
			                                                index_in_L_space);
			++index_in_x;
			++index_in_L;
		} else if (_x_sm[index_in_x].get_linear_index_in_container_space() <
		           _L_sm[index_in_L].get_linear_index_in_container_space()) {
			const auto index_in_L_space = _x_sm[index_in_x].get_linear_index_in_container_space());
			sum_log.add(_calculate_probability_of_L_given_x(_x_sm[index_in_x].is_one(), false,
			                                               index_in_L_space);
			++index_in_x;
		} else {
			const auto index_in_L_space = _L_sm[index_in_L].get_linear_index_in_container_space());
			sum_log.add(_calculate_probability_of_L_given_x(false, _L_sm[index_in_L].is_one(),
			                                                index_in_L_space);
			++index_in_L;
		}
	}
	return sum_log.getSum();
}

double TLotus::_calculate_probability_of_L_given_x(bool x, bool L,
                                                   const std::vector<size_t> &index_in_collapsed_space) const {
	if (x && L) { return _calculate_research_effort(index_in_collapsed_space); }
	if (x && !L) { return 1.0 - _calculate_research_effort(index_in_collapsed_space); }
	if (!x && L) { return 0.0; }
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

void TLotus::guessInitialValues() {
	// TODO: think of an initialization scheme
	_gamma->set(0.001);
}

void TLotus::_simulateUnderPrior(Storage *) {
	// TODO: implement simulator
}
