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

	_L.initialize_dimensions(sparse_data_file::leaf_shape(_trees));
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

double TLotus::ll_ratio_after_parameter_move() {
	// One function for both gamma and the error rate: the reporting model is built from both, so a
	// move on either replaces it wholesale. Neither moves the field, so the present cells the last
	// block update counted are the field's, and the likelihood is a sum over their kinds rather
	// than a pass over the cells (ADR-0015).
	_oldLL           = _curLL;
	_reporting_model = _build_reporting_model();
	_curLL           = _log_likelihood_of_present_cells();
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

	install_present_cells(count_present_cells(Y));
	_oldLL = _curLL;
}

void TLotus::_gather_paper_counts() {
	_paper_counts.resize(_trees.size());
	for (size_t i = 0; i < _trees.size(); ++i) {
		const auto &tree = *_trees[i];
		_paper_counts[i] = read_paper_counts(tree.get_tree_name(), tree.phylogeny());
	}
	_classes = lotus_math::TPaperCountClasses(_paper_counts);
}

lotus_math::TReportingModel TLotus::_build_reporting_model() const {
	return {_gammas(), (double)_error_rate->value(), _paper_counts};
}

std::vector<double> TLotus::_gammas() const {
	std::vector<double> gammas(_trees.size());
	for (size_t i = 0; i < gammas.size(); ++i) { gammas[i] = (double)_gamma->value(i); }
	return gammas;
}

double TLotus::_log_likelihood_of_present_cells() const {
	return lotus_math::log_likelihood_of_L(_present_cells, _classes, _gammas(),
	                                       (double)_error_rate->value());
}

lotus_math::TPresentCells TLotus::no_present_cells() const {
	return {_L.total_size_of_container_space(), _L.number_of_ones(), _classes};
}

void TLotus::install_present_cells(lotus_math::TPresentCells cells) {
	_present_cells = std::move(cells);
	_curLL         = _log_likelihood_of_present_cells();
}

lotus_math::TPresentCells TLotus::count_present_cells(const TFieldStorage &Y) const {
	static_assert(NUMBER_OF_TREES == 2, "a cell is a species leaf by a molecule leaf");
	lotus_math::TPresentCells cells = no_present_cells();
	const size_t n_columns          = Y.dimensions()[1];

	// Only the present cells are counted, so the walk is over the field's ones. Whether LOTUS
	// reports one is a forward cursor over its records, which are in the same order.
	auto record = _L.ones_cursor();
	for (auto present = Y.ones_cursor(); present.valid(); present.advance()) {
		const size_t cell = present.linear_index();
		while (record.valid() && record.linear_index() < cell) { record.advance(); }
		const bool reported = record.valid() && record.linear_index() == cell;
		const size_t slot   = cells.slot(species_class(cell / n_columns),
		                                 molecule_class(cell % n_columns));
		++cells.ones[static_cast<size_t>(reported)][slot];
	}
	return cells;
}

void TLotus::prepare_for_simulation() {
	// L is indexed on every tree, so its shape follows from the trees alone.
	_L.initialize_dimensions(sparse_data_file::leaf_shape(_trees));

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

void TLotus::simulate_from_Y(const TFieldStorage &Y) {
	for (size_t i = 0; i < _L.total_size_of_container_space(); ++i) {
		const auto multi_dim_index_in_L_space = _L.get_multi_dimensional_index(i);
		// L has the field's dimensions, so cell i of one is cell i of the other (a missing cell
		// reads as 0).
		const bool x                          = Y.is_one(i);
		const double proba = _reporting().probability(x, true, multi_dim_index_in_L_space);
		const coretools::Probability p(proba);
		if (coretools::instances::randomGenerator().pickOneOfTwo(p)) { _L.insert_one(i); }
	}
}

void TLotus::write_simulated(const std::string &prefix) const {
	const std::string file_name = prefix + "_simulated_lotus.tsv";

	// we get the tree name for the header of the file.
	const auto header = tree_names();

	coretools::TOutputFile file(file_name, header, "\t");
	std::vector<std::string> line(_trees.size());
	// The records that are one, in ascending linear-index order. A record is a one by definition,
	// so the walk needs no test of what it is handed.
	for (auto record = _L.ones_cursor(); record.valid(); record.advance()) {
		// for each draw we need to get the node name of the leaf in the correct tree and write it
		const auto multi_dim_index_in_L_space =
		    _L.get_multi_dimensional_index(record.linear_index());
		for (size_t j = 0; j < _trees.size(); ++j) {
			const size_t node_index_in_tree =
			    _trees[j]->get_node_index_from_leaf_index(multi_dim_index_in_L_space[j]);
			line[j] = _trees[j]->get_node_id(node_index_in_tree);
		}
		file.writeln(line);
	}
}

void TLotus::contribute_stats(TNotifierStats &stats) const {
	stats.dim_names   = tree_names();
	stats.gamma_stats = gamma_stats();
	stats.scalar_stats.push_back({"epsilon", error_rate_stats()});
}

#endif // USE_LOTUS
