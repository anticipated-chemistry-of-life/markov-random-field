//
// Created by madleina on 03.03.25.
//
#include "lotus/TLotus.h"

#ifdef USE_LOTUS

#include "TSparseDataFile.h"
#include "Types.h"
#include "cli.h"
#include "constants.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Files/TOutputFile.h"
#include "coretools/Main/TError.h"
#include "coretools/Types/probability.h"
#include "lotus/paper_counts.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

TLotus::TLotus(const std::vector<std::unique_ptr<TTree>> &trees, TypeParamGamma *gamma,
               TypeParamErrorRate *error_rate)
    : _trees(trees), _gamma(gamma), _error_rate(error_rate) {}

void TLotus::initialize(TDataModel *box, bool simulate) {
	if (!simulate) { load_from_file(get_filename_lotus()); }

	// Note: when simulating, the collapser has not been initialized yet (that happens in
	// load_from_file, which is skipped), so gamma is sized to 0 here and re-initialized in
	// prepare_for_simulation once the kept dimensions are known.
	_gamma->initStorage(box, {NUMBER_OF_TREES},
	                    {std::make_shared<coretools::TNamesStrings>(kept_tree_names())});
	_error_rate->initStorage(box, {1});
}

void TLotus::load_from_file(const std::string &filename) {
	coretools::instances::logfile().startIndent("Reading links from file '", filename, "' ...");
	coretools::TInputFile file(filename, coretools::FileType::Header);

	// LOTUS may name fewer trees than exist; the remaining dimensions are collapsed away.
	sparse_data_file::validate_header_against_trees(file, _trees, filename);

	// initialize the size of L
	std::vector<size_t> len_per_dimension_lotus;
	for (size_t i = 0; i < _trees.size(); ++i) {
		len_per_dimension_lotus.push_back(_trees[i]->get_number_of_leaves());
	}
	_L.initialize(len_per_dimension_lotus);

	_gather_paper_counts();

	IndexArray index{};
	for (; !file.empty(); file.popFront()) {
		// loop over all columns
		for (size_t i = 0; i < NUMBER_OF_TREES; ++i) {
			index[i] = sparse_data_file::leaf_index_or_throw(*_trees[i], std::string(file.get(i)));
		}
		_L.insert_one(index);
	}
	coretools::instances::logfile().endIndent();
}

std::vector<std::string> TLotus::kept_tree_names() const {
	std::vector<std::string> tree_names;
	tree_names.reserve(NUMBER_OF_TREES);
	for (size_t i = 0; i < NUMBER_OF_TREES; ++i) {
		tree_names.push_back(_trees[i]->get_tree_name());
	}
	return tree_names;
}

std::vector<TNtfyNotifier::ParamStats> TLotus::gamma_stats() const {
	std::vector<TNtfyNotifier::ParamStats> stats;
	stats.reserve(NUMBER_OF_TREES);
	for (size_t i = 0; i < NUMBER_OF_TREES; ++i) {
		stats.push_back({_gamma->mean(i), _gamma->var(i), _gamma->sd(i)});
	}
	return stats;
}

TNtfyNotifier::ParamStats TLotus::error_rate_stats() const {
	return {_error_rate->mean(0), _error_rate->var(0), _error_rate->sd(0)};
}

double TLotus::calculate_log_likelihood_of_L(const TFieldStorage &Y) const {
	return _calculate_log_likelihood_of_L(Y);
}

/// This function will be used when we update Y.
void TLotus::calculate_LL_update_Y(const IndexArray &index_in_leaves_space,
                                   std::array<double, 2> &prob) const {
	// new Y = 0 -> x_is_one_for_Y_0 will always be false here (because of the previous
	// if-statement) new Y = 1 -> x will always be true
	for (size_t i = 0; i < 2; ++i) {
		prob[i] =
		    _reporting().probability(i, _L.is_one(index_in_leaves_space).is_one, index_in_leaves_space);
	}
}

double TLotus::ll_ratio_after_parameter_move(const TFieldStorage &Y) {
	// One function for both gamma and the error rate: the reporting model is built from both, so a
	// move on either replaces it wholesale. Rebuilding the factor table on an error-rate move is
	// strictly redundant -- the factors depend only on gamma -- but it is one exp() per leaf
	// against a likelihood sweep over the whole container space, and it buys a single code path
	// with no parameter-specific state.
	_oldLL           = _curLL;
	_reporting_model = _build_reporting_model();
	_curLL           = calculate_log_likelihood_of_L(Y);
	return _curLL - _oldLL;
}

void TLotus::revert_parameter_move() {
	// The candidate model is simply dropped: stattools has already restored the parameter, so the
	// next build reads the old value. Only the cached likelihood has to be put back by hand.
	_curLL           = _oldLL;
	_reporting_model = _build_reporting_model();
}

void TLotus::guess_initial_values(const TFieldStorage &Y) {
	for (size_t i = 0; i < NUMBER_OF_TREES; ++i) { _gamma->set(i, ProgramOptions::GAMMA); }
	_error_rate->set(ProgramOptions::EPSILON);

	_reporting_model = _build_reporting_model(); // parameters just set -> build before the LL

	// initialize _curLL
	_curLL = calculate_log_likelihood_of_L(Y);
	_oldLL = _curLL;
}

void TLotus::_gather_paper_counts() {
	// for example, size is 2 if keep molecules and species
	_paper_counts.resize(NUMBER_OF_TREES);
	for (size_t i = 0; i < NUMBER_OF_TREES; ++i) {
		const auto &tree = *_trees[i];
		_paper_counts[i] = read_paper_counts(tree.get_tree_name(), tree.phylogeny());
	}
}

lotus_math::TReportingModel TLotus::_build_reporting_model() const {
	std::vector<double> gammas(NUMBER_OF_TREES);
	for (size_t i = 0; i < gammas.size(); ++i) { gammas[i] = (double)_gamma->value(i); }
	return {gammas, (double)_error_rate->value(), _paper_counts};
}

double TLotus::_calculate_log_likelihood_of_L(const TFieldStorage &Y) const {
	const size_t total = Y.total_size_of_container_space();

	// Merge-join the two fields in ascending linear-index order without materializing their
	// entries. We only need to evaluate cells that are one in Y and/or L: for every other cell
	// both states are false, and _calculate_probability_of_L_given_x(false, false, i) is
	// position-independent, so those collapse into a single bulk term below.
	//
	// The ones and not the stored cells, because how the sum splits between the accumulator and
	// that bulk term has to be a property of the field's contents rather than of which backend is
	// holding them -- otherwise the same chain reaches answers that differ in the last bits under
	// the two, which a Metropolis ratio is quite capable of turning into two different chains.
	coretools::TSumLogProbability sum_log;

	for (size_t i = 0; i < total; ++i) {
		bool state_of_Y = Y.is_one(i).is_one;
		bool state_of_L = _L.is_one(i).is_one;

		sum_log.add(
		    _reporting().probability(state_of_Y, state_of_L, _L.get_multi_dimensional_index(i)));
	}
	return sum_log.getSum();
}

void TLotus::prepare_for_simulation(TDataModel *box) {
	// by default, we keep all the trees
	std::vector<std::string> tree_names_to_keep_default;
	tree_names_to_keep_default.reserve(_trees.size());
	for (const auto &tree : _trees) { tree_names_to_keep_default.push_back(tree->get_tree_name()); }

	// else get the tree names to keep from CLI
	std::vector<std::string> tree_names_to_keep =
	    coretools::instances::parameters().get("tree_names_to_keep", tree_names_to_keep_default);

	// initialize the size of L
	std::vector<size_t> len_per_dimension_lotus;
	for (size_t i = 0; i < _trees.size(); ++i) {
		len_per_dimension_lotus.push_back(_trees[i]->get_number_of_leaves());
	}
	_L.initialize(len_per_dimension_lotus);

	// initialize the error rate
	const auto error_rate =
	    coretools::instances::parameters().get<double>("error_rate", ProgramOptions::EPSILON);
	_error_rate->set(error_rate);

	// initialize the gamma parameters. Since we don't read Lotus from a file, the size of the
	// gamma is never initilized. That is why we need to do it here.
	_gamma->initStorage(box, {NUMBER_OF_TREES},
	                    {std::make_shared<coretools::TNamesStrings>(tree_names_to_keep)});
	const auto gamma = ProgramOptions::GAMMA;
	for (size_t i = 0; i < NUMBER_OF_TREES; ++i) { _gamma->set(i, gamma); }

	// 2025.06.16 after discussion of last week with Dan, we should be able to also simuate and
	// provide the number of papers prior to the simulation.
	_gather_paper_counts();

	_reporting_model = _build_reporting_model(); // counts + parameters ready -> ready to simulate L
}

void TLotus::simulate_L_from_Y(const TFieldStorage &Y) {
	for (size_t i = 0; i < _L.total_size_of_container_space(); ++i) {
		const auto multi_dim_index_in_L_space = _L.get_multi_dimensional_index(i);
		bool x                                = Y.is_one(i).is_one;
		const double proba = _reporting().probability(x, true, multi_dim_index_in_L_space);
		const coretools::Probability p(proba);
		if (coretools::instances::randomGenerator().pickOneOfTwo(p)) { _L.insert_one(i); }
	}
}

void TLotus::write_simulated_L(const std::string &prefix) const {
	const std::string file_name = prefix + "_simulated_lotus.tsv";

	// we get the tree name for the header of the file.
	const auto header = kept_tree_names();

	coretools::TOutputFile file(file_name, header, "\t");
	std::vector<std::string> line(NUMBER_OF_TREES);

	for (size_t i = 0; i < _L.total_size_of_container_space(); ++i) {
		if (!_L.is_one(i).is_one) { continue; }
		const auto multi_dim_index_in_L_space = _L.get_multi_dimensional_index(i);
		for (size_t j = 0; j < NUMBER_OF_TREES; ++j) {
			const size_t node_index_in_tree =
			    _trees[j]->get_node_index_from_leaf_index(multi_dim_index_in_L_space[j]);
			line[j] = _trees[j]->get_node_id(node_index_in_tree);
		}
		file.writeln(line);
	}
}

#endif // USE_LOTUS
