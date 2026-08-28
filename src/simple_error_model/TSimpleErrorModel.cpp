//
// Created for the split of the data sources (LOTUS / simple error model / MS data).
//
#include "simple_error_model/TSimpleErrorModel.h"

#ifdef USE_SIMPLE_ERROR_MODEL

#include "TSparseDataFile.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Files/TOutputFile.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TLog.h"
#include <cstddef>
#include <string>
#include <vector>

TSimpleErrorModel::TSimpleErrorModel(const std::vector<std::unique_ptr<TTree>> &trees,
                                     TypeParamEpsilon *epsilon)
    : _trees(trees), _epsilon(epsilon), _tmp_state_along_last_dim(trees.back()->phylogeny(), 1) {}

void TSimpleErrorModel::initialize_storage() {
	// D lives in the same space as Y: one dimension per tree, sized by that tree's leaf count.
	std::vector<size_t> num_leaves_per_dim;
	num_leaves_per_dim.reserve(_trees.size());
	for (const auto &tree : _trees) { num_leaves_per_dim.push_back(tree->get_number_of_leaves()); }

	// n_iterations = 1: D is data, its per-cell counters are never used.
	_D.initialize(1, num_leaves_per_dim);
	_total_cells = _D.total_size_of_container_space();
}

void TSimpleErrorModel::load_from_file(const std::string &filename) {
	coretools::instances::logfile().startIndent("Reading simple error model data from file '",
	                                            filename, "' ...");
	coretools::TInputFile file(filename, coretools::FileType::Header);

	// Unlike LOTUS, D is never collapsed: it must name every tree, in tree order. The order is
	// checked by validate_header_against_trees, so pinning the count is enough to make column i
	// correspond to tree i.
	if (file.header().size() != _trees.size()) {
		throw coretools::TUserError("File '", filename, "' has ", file.header().size(),
		                            " columns but there are ", _trees.size(),
		                            " trees. The simple error model data is never collapsed, so it "
		                            "must have exactly one column per tree.");
	}
	sparse_data_file::validate_header_against_trees(file, _trees, filename);

	IndexArray index_in_D_space{};
	for (; !file.empty(); file.popFront()) {
		for (size_t i = 0; i < _trees.size(); ++i) {
			index_in_D_space[i] =
			    sparse_data_file::leaf_index_or_throw(*_trees[i], std::string(file.get(i)));
		}
		_D.insert_one(_D.get_linear_index_in_container_space(index_in_D_space));
	}

	coretools::instances::logfile().list("Read ", _D.number_of_ones(), " ones out of ",
	                                     _total_cells, " cells.");
	coretools::instances::logfile().endIndent();
}

void TSimpleErrorModel::guess_initial_values(const TFieldStorage &Y) {
	_n_disagree = simple_error_model::count_disagreements(Y, _D);
}

void TSimpleErrorModel::fill_tmp_state_along_last_dim(const IndexArray &start_index_in_leaves_space,
                                                      size_t K) {
	// D has the same dimensions as Y, so the index needs no translation (in contrast to LOTUS,
	// where the collapser maps Y space to L space first).
	_tmp_state_along_last_dim.fill_Y_along_last_dim(start_index_in_leaves_space, K, _D);
}

void TSimpleErrorModel::simulate_D_from_Y(const TFieldStorage &Y) {
	if (Y.dimensions() != _D.dimensions()) {
		throw coretools::TDevError("Cannot simulate the simple error model data: Y is ",
		                           Y.dimensions()[0], "x", Y.dimensions()[1], " but D is ",
		                           _D.dimensions()[0], "x", _D.dimensions()[1], ".");
	}

	const double eps = _eps();
	for (size_t i = 0; i < _total_cells; ++i) {
		// Only ones are stored; an absent cell reads as 0, which is exactly the encoding the
		// sparse output file uses.
		if (simple_error_model::draw_D_given_Y(Y.is_one(i), eps)) { _D.insert_one(i); }
	}

	_n_disagree = simple_error_model::count_disagreements(Y, _D);
	coretools::instances::logfile().list("Simulated simple error model data with epsilon = ", eps,
	                                     ": ", _n_disagree, " of ", _total_cells,
	                                     " cells differ from Y.");
}

void TSimpleErrorModel::write_simulated_D(const std::string &prefix) const {
	const std::string file_name = prefix + "_simulated_simple_data.tsv";
	coretools::TOutputFile file(file_name, sparse_data_file::header_from_trees(_trees), "\t");

	std::vector<std::string> line(_trees.size());
	for (const auto &[linear_index, storage] : _D.get_stored_entries()) {
		if (!storage.is_one()) { continue; }
		const auto index_in_D_space = _D.get_multi_dimensional_index(linear_index);
		for (size_t i = 0; i < _trees.size(); ++i) {
			const size_t node_index =
			    _trees[i]->get_node_index_from_leaf_index(index_in_D_space[i]);
			line[i] = _trees[i]->get_node_id(node_index);
		}
		file.writeln(line);
	}
}

#endif // USE_SIMPLE_ERROR_MODEL
