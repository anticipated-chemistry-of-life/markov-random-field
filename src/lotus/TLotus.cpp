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
    : _trees(trees), _collapser(trees), _gamma(gamma), _error_rate(error_rate),
      _tmp_state_along_last_dim(*trees.back().get(), 1) {}

void TLotus::initialize(TDataModel *box, bool simulate) {
	if (!simulate) { load_from_file(get_filename_lotus()); }

	// Note: when simulating, the collapser has not been initialized yet (that happens in
	// load_from_file, which is skipped), so gamma is sized to 0 here and re-initialized in
	// prepare_for_simulation once the kept dimensions are known.
	_gamma->initStorage(box, {_collapser.num_dim_to_keep()},
	                    {std::make_shared<coretools::TNamesStrings>(kept_tree_names())});
	_error_rate->initStorage(box, {1});
}

void TLotus::load_from_file(const std::string &filename) {
	coretools::instances::logfile().startIndent("Reading links from file '", filename, "' ...");
	coretools::TInputFile file(filename, coretools::FileType::Header);

	// LOTUS may name fewer trees than exist; the remaining dimensions are collapsed away.
	sparse_data_file::validate_header_against_trees(file, _trees, filename);

	// initialize collapser: know which dimensions to keep and which to collapse
	const auto len_per_dimension_lotus = _collapser.initialize(file.header(), "LOTUS");

	// initialize the size of L
	_L.initialize(1, len_per_dimension_lotus);

	_gather_paper_counts();

	IndexArray index_in_collapsed_space{};
	for (; !file.empty(); file.popFront()) {
		// loop over all columns
		for (size_t i = 0; i < _collapser.num_dim_to_keep(); ++i) {
			const size_t tree_index     = _collapser.dim_to_keep(i);
			index_in_collapsed_space[i] = sparse_data_file::leaf_index_or_throw(
			    *_trees[tree_index], std::string(file.get(i)));
		}
		_L.insert_one(_L.get_linear_index_in_container_space(index_in_collapsed_space));
	}
	coretools::instances::logfile().endIndent();
}

std::vector<std::string> TLotus::kept_tree_names() const {
	std::vector<std::string> tree_names;
	tree_names.reserve(_collapser.num_dim_to_keep());
	for (size_t i = 0; i < _collapser.num_dim_to_keep(); ++i) {
		tree_names.push_back(_trees[_collapser.dim_to_keep(i)]->get_tree_name());
	}
	return tree_names;
}

std::vector<TNtfyNotifier::ParamStats> TLotus::gamma_stats() const {
	std::vector<TNtfyNotifier::ParamStats> stats;
	stats.reserve(_collapser.num_dim_to_keep());
	for (size_t i = 0; i < _collapser.num_dim_to_keep(); ++i) {
		stats.push_back({_gamma->mean(i), _gamma->var(i), _gamma->sd(i)});
	}
	return stats;
}

TNtfyNotifier::ParamStats TLotus::error_rate_stats() const {
	return {_error_rate->mean(0), _error_rate->var(0), _error_rate->sd(0)};
}

double TLotus::calculate_log_likelihood_of_L(const TStorageYMatrix &Y) const {
	if (_collapser.do_collapse()) {
		throw coretools::TDevError("calculate_log_likelihood_of_L: do_collapse() is true. This "
		                           "part has not been implemented yet.");
	}
	return _calculate_log_likelihood_of_L_no_collapsing(Y);
}

void TLotus::fill_tmp_state_along_last_dim(const IndexArray &start_index_clique_along_last_dim,
                                           size_t K) {
	// collapse start_index_in_leaves (this is the index in Y)
	if (_collapser.do_collapse()) {
		_tmp_state_along_last_dim.fill_Y_along_last_dim(
		    _collapser.collapse(start_index_clique_along_last_dim), K, _L);
	} else { // no need to collapse
		_tmp_state_along_last_dim.fill_Y_along_last_dim(start_index_clique_along_last_dim, K, _L);
	}
}

/// This function will be used when we update Y.
void TLotus::calculate_LL_update_Y(const IndexArray &index_in_leaves_space,
                                   size_t index_for_tmp_state, bool old_state,
                                   std::array<double, 2> &prob) const {
	// function gets the old_state and needs to calculate LL for new_state = 0 and 1
	// for state 1, we know that the new x will always be 1 (at least one is a one)
	bool x_is_one_for_Y_0 = false; // Y = 0 -> x = 0 if we don't collapse
	if (_collapser.do_collapse()) {
		x_is_one_for_Y_0 = _collapser.x_is_one(index_in_leaves_space, old_state);
	}

	if (x_is_one_for_Y_0) { // Y=0 also results in x=1 (others in clique are a one)
		// x is one for both states (due to collapsing) -> likelihood doesn't matter. `prob` is left
		// at its caller-provided neutral value (1.0), so both states get the same zero
		// contribution.
		return;
	}

	const auto index_in_collapsed_space = _collapser.collapse(index_in_leaves_space);
	// new Y = 0 -> x_is_one_for_Y_0 will always be false here (because of the previous
	// if-statement) new Y = 1 -> x will always be true
	for (size_t i = 0; i < 2; ++i) {
		prob[i] = _reporting().probability(i, _tmp_state_along_last_dim.get_Y(index_for_tmp_state),
		                                   index_in_collapsed_space);
	}
}

double TLotus::ll_ratio_after_parameter_move(const TStorageYMatrix &Y) {
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

void TLotus::guess_initial_values(const TStorageYMatrix &Y) {
	for (size_t i = 0; i < _collapser.num_dim_to_keep(); ++i) {
		_gamma->set(i, ProgramOptions::GAMMA);
	}
	_error_rate->set(ProgramOptions::EPSILON);

	_reporting_model = _build_reporting_model(); // parameters just set -> build before the LL

	// initialize _curLL
	_curLL = calculate_log_likelihood_of_L(Y);
	_oldLL = _curLL;
}

void TLotus::_gather_paper_counts() {
	// for example, size is 2 if keep molecules and species
	_paper_counts.resize(_collapser.num_dim_to_keep());
	for (size_t i = 0; i < _collapser.num_dim_to_keep(); ++i) {
		const auto &tree = *_trees[_collapser.dim_to_keep(i)];
		_paper_counts[i] = read_paper_counts(tree.get_tree_name(), tree.phylogeny());
	}
}

lotus_math::TReportingModel TLotus::_build_reporting_model() const {
	std::vector<double> gammas(_collapser.num_dim_to_keep());
	for (size_t i = 0; i < gammas.size(); ++i) { gammas[i] = (double)_gamma->value(i); }
	return {gammas, (double)_error_rate->value(), _paper_counts};
}

double TLotus::_calculate_log_likelihood_of_L_no_collapsing(const TStorageYMatrix &Y) const {
	const size_t total = Y.total_size_of_container_space();

	// Merge-join the two sparse matrices in ascending linear-index order without materializing
	// their entries. We only need to evaluate cells that are stored in Y and/or L: for every other
	// cell both states are false, and _calculate_probability_of_L_given_x(false, false, i) is
	// position-independent, so those collapse into a single bulk term below.
	coretools::TSumLogProbability sum_log;
	auto y_cur = Y.stored_cursor();
	auto l_cur = _L.stored_cursor();

	size_t n_visited = 0; // distinct linear indices stored in Y and/or L
	while (y_cur.valid() || l_cur.valid()) {
		const size_t yi = y_cur.valid() ? y_cur.linear_index() : total;
		const size_t li = l_cur.valid() ? l_cur.linear_index() : total;
		const size_t i  = std::min(yi, li);

		bool state_of_Y = false;
		bool state_of_L = false;
		if (yi == i) {
			state_of_Y = y_cur.is_one();
			y_cur.advance();
		}
		if (li == i) {
			state_of_L = l_cur.is_one();
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

void TLotus::prepare_for_simulation(TDataModel *box) {
	// by default, we keep all the trees
	std::vector<std::string> tree_names_to_keep_default;
	tree_names_to_keep_default.reserve(_trees.size());
	for (const auto &tree : _trees) { tree_names_to_keep_default.push_back(tree->get_tree_name()); }

	// else get the tree names to keep from CLI
	std::vector<std::string> tree_names_to_keep =
	    coretools::instances::parameters().get("tree_names_to_keep", tree_names_to_keep_default);

	const auto len_per_dimension_lotus = _collapser.initialize(tree_names_to_keep, "LOTUS");
	// initialize the size of L
	_L.initialize(1, len_per_dimension_lotus);

	// initialize the error rate
	const auto error_rate =
	    coretools::instances::parameters().get<double>("error_rate", ProgramOptions::EPSILON);
	_error_rate->set(error_rate);

	// initialize the gamma parameters. Since we don't read Lotus from a file, the size of the
	// gamma is never initilized. That is why we need to do it here.
	_gamma->initStorage(box, {_collapser.num_dim_to_keep()},
	                    {std::make_shared<coretools::TNamesStrings>(tree_names_to_keep)});
	const auto gamma = ProgramOptions::GAMMA;
	for (size_t i = 0; i < _collapser.num_dim_to_keep(); ++i) { _gamma->set(i, gamma); }

	// 2025.06.16 after discussion of last week with Dan, we should be able to also simuate and
	// provide the number of papers prior to the simulation.
	_gather_paper_counts();

	_reporting_model = _build_reporting_model(); // counts + parameters ready -> ready to simulate L
}

void TLotus::simulate_L_from_Y(const TStorageYMatrix &Y) {
	for (size_t i = 0; i < _L.total_size_of_container_space(); ++i) {
		const auto multi_dim_index_in_L_space = _L.get_multi_dimensional_index(i);
		bool x                                = false;
		if (_collapser.do_collapse()) {
			x = _collapser.x_is_one(multi_dim_index_in_L_space);
		} else {
			// When we don't collapse, the dimensions of Y equal the dimensions of Lotus, so we can
			// look up the state of Y directly (a missing cell reads as 0).
			x = Y.is_one(i);
		}
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
	std::vector<std::string> line(_collapser.num_dim_to_keep());
	for (const auto &[linear_index, storage] : _L.get_stored_entries()) {
		if (!storage.is_one()) { continue; }
		// for each draw we need to get the node name of the leaf in the correct tree and write it
		const auto multi_dim_index_in_L_space = _L.get_multi_dimensional_index(linear_index);
		for (size_t j = 0; j < _collapser.num_dim_to_keep(); ++j) {
			const size_t tree_index = _collapser.dim_to_keep(j);
			const size_t node_index_in_tree =
			    _trees[tree_index]->get_node_index_from_leaf_index(multi_dim_index_in_L_space[j]);
			line[j] = _trees[tree_index]->get_node_id(node_index_in_tree);
		}
		file.writeln(line);
	}
}

#endif // USE_LOTUS
