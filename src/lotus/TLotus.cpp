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

	// One gamma per tree, whether or not a file was read: L is indexed on every tree, so the
	// count no longer depends on anything load_from_file discovers.
	_gamma->initStorage(box, {_trees.size()},
	                    {std::make_shared<coretools::TNamesStrings>(tree_names())});
	_error_rate->initStorage(box, {1});
}

void TLotus::load_from_file(const std::string &filename) {
	coretools::instances::logfile().startIndent("Reading links from file '", filename, "' ...");
	coretools::TInputFile file(filename, coretools::FileType::Header);

	// L names every tree, in tree order, so column i is tree i.
	sparse_data_file::validate_header_against_trees(file, _trees, filename);

	_L.initialize(1, _leaf_counts());
	_gather_paper_counts();

	IndexArray index_in_L_space{};
	for (; !file.empty(); file.popFront()) {
		for (size_t i = 0; i < _trees.size(); ++i) {
			index_in_L_space[i] =
			    sparse_data_file::leaf_index_or_throw(*_trees[i], std::string(file.get(i)));
		}
		_L.insert_one(_L.get_linear_index_in_container_space(index_in_L_space));
	}
	coretools::instances::logfile().endIndent();
}

std::vector<std::string> TLotus::tree_names() const {
	std::vector<std::string> names;
	names.reserve(_trees.size());
	for (const auto &tree : _trees) { names.push_back(tree->get_tree_name()); }
	return names;
}

std::vector<size_t> TLotus::_leaf_counts() const {
	std::vector<size_t> counts;
	counts.reserve(_trees.size());
	for (const auto &tree : _trees) { counts.push_back(tree->get_number_of_leaves()); }
	return counts;
}

std::vector<TNtfyNotifier::ParamStats> TLotus::gamma_stats() const {
	std::vector<TNtfyNotifier::ParamStats> stats;
	stats.reserve(_trees.size());
	for (size_t i = 0; i < _trees.size(); ++i) {
		stats.push_back({_gamma->mean(i), _gamma->var(i), _gamma->sd(i)});
	}
	return stats;
}

TNtfyNotifier::ParamStats TLotus::error_rate_stats() const {
	return {_error_rate->mean(0), _error_rate->var(0), _error_rate->sd(0)};
}

/// This function will be used when we update Y.
void TLotus::calculate_LL_update_Y(const IndexArray &index_in_leaves_space, bool reports_the_cell,
                                   std::array<double, 2> &prob) const {
	// A field cell is its own L cell now, so the emission depends on this cell alone: no clique of
	// other cells can make L's latent state one while this one is zero.
	for (size_t i = 0; i < 2; ++i) {
		prob[i] = _reporting().probability(i, reports_the_cell, index_in_leaves_space);
	}
}

double TLotus::ll_ratio_after_parameter_move(const TFieldStorage &Y) {
	// One function for both gamma and the error rate: the reporting model is built from both, so a
	// move on either replaces it wholesale. Rebuilding the factor table on an error-rate move is
	// strictly redundant -- the factors depend only on gamma -- but it is one exp() per leaf
	// against a likelihood sum over the whole container space, and it buys a single code path
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
	for (size_t i = 0; i < _trees.size(); ++i) { _gamma->set(i, ProgramOptions::GAMMA); }
	_error_rate->set(ProgramOptions::EPSILON);

	_reporting_model = _build_reporting_model(); // parameters just set -> build before the LL

	// initialize _curLL
	_curLL = calculate_log_likelihood_of_L(Y);
	_oldLL = _curLL;
}

void TLotus::_gather_paper_counts() {
	_paper_counts.resize(_trees.size());
	for (size_t i = 0; i < _trees.size(); ++i) {
		const auto &tree = *_trees[i];
		_paper_counts[i] = read_paper_counts(tree.get_tree_name(), tree.phylogeny());
	}
}

lotus_math::TReportingModel TLotus::_build_reporting_model() const {
	std::vector<double> gammas(_trees.size());
	for (size_t i = 0; i < gammas.size(); ++i) { gammas[i] = (double)_gamma->value(i); }
	return {gammas, (double)_error_rate->value(), _paper_counts};
}

double TLotus::calculate_log_likelihood_of_L(const TFieldStorage &Y) const {
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
	auto y_cur = Y.ones_cursor();
	auto l_cur = _L.ones_cursor();

	size_t n_visited = 0; // distinct linear indices that are one in Y and/or L
	while (y_cur.valid() || l_cur.valid()) {
		const size_t yi = y_cur.valid() ? y_cur.linear_index() : total;
		const size_t li = l_cur.valid() ? l_cur.linear_index() : total;
		const size_t i  = std::min(yi, li);

		bool state_of_Y = false;
		bool state_of_L = false;
		if (yi == i) {
			state_of_Y = true; // the cursor yields only ones
			y_cur.advance();
		}
		if (li == i) {
			state_of_L = true;
			l_cur.advance();
		}

		// Only a present cell needs its position; for an absent one the answer is the same
		// everywhere, so the linear-to-multidimensional conversion is skipped entirely.
		if (state_of_Y) {
			sum_log.add(
			    _reporting().probability(true, state_of_L, _L.get_multi_dimensional_index(i)));
		} else {
			sum_log.add(_reporting().probability_absent(state_of_L));
		}
		++n_visited;
	}

	// Every remaining cell is (Y = 0, L = 0): a single constant, position-independent term.
	const double p_absent = _reporting().probability_absent(false);
	return sum_log.getSum() + static_cast<double>(total - n_visited) * std::log(p_absent);
}

void TLotus::prepare_for_simulation() {
	// L is indexed on every tree, so its shape follows from the trees alone.
	_L.initialize(1, _leaf_counts());

	// initialize the error rate
	const auto error_rate =
	    coretools::instances::parameters().get<double>("error_rate", ProgramOptions::EPSILON);
	_error_rate->set(error_rate);

	// gamma was already sized in initialize(), which no longer has to wait for a file to tell it
	// how many dimensions L has.
	const auto gamma = ProgramOptions::GAMMA;
	for (size_t i = 0; i < _trees.size(); ++i) { _gamma->set(i, gamma); }

	// 2025.06.16 after discussion of last week with Dan, we should be able to also simuate and
	// provide the number of papers prior to the simulation.
	_gather_paper_counts();

	_reporting_model = _build_reporting_model(); // counts + parameters ready -> ready to simulate L
}

void TLotus::simulate_L_from_Y(const TFieldStorage &Y) {
	for (size_t i = 0; i < _L.total_size_of_container_space(); ++i) {
		const auto multi_dim_index_in_L_space = _L.get_multi_dimensional_index(i);
		// L has the field's dimensions, so cell i of one is cell i of the other (a missing cell
		// reads as 0).
		const bool x       = Y.is_one(i);
		const double proba = _reporting().probability(x, true, multi_dim_index_in_L_space);
		const coretools::Probability p(proba);
		if (coretools::instances::randomGenerator().pickOneOfTwo(p)) { _L.insert_one(i); }
	}
}

void TLotus::write_simulated_L(const std::string &prefix) const {
	const std::string file_name = prefix + "_simulated_lotus.tsv";

	// we get the tree name for the header of the file.
	const auto header = tree_names();

	coretools::TOutputFile file(file_name, header, "\t");
	std::vector<std::string> line(_trees.size());
	for (const auto &[linear_index, storage] : _L.get_stored_entries()) {
		if (!storage.is_one()) { continue; }
		// for each draw we need to get the node name of the leaf in the correct tree and write it
		const auto multi_dim_index_in_L_space = _L.get_multi_dimensional_index(linear_index);
		for (size_t j = 0; j < _trees.size(); ++j) {
			const size_t node_index_in_tree =
			    _trees[j]->get_node_index_from_leaf_index(multi_dim_index_in_L_space[j]);
			line[j] = _trees[j]->get_node_id(node_index_in_tree);
		}
		file.writeln(line);
	}
}

#endif // USE_LOTUS
