//
// Created by madleina on 10.01.25.
//

#include "TMarkovField.h"
#include "TClique.h"
#include "TDataModel.h"
#include "Types.h"
#include "cli.h"
#include "constants.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/progressTools.h"
#include "coretools/algorithms.h"
#include "omp.h"
#include "random/TCellUniforms.h"
#include "storages/storage_backend.h"
#include "tree/TTree.h"
#include "tree/io/write_Z.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

// The field update names the two trees. A thread takes a species leaf and walks the molecule
// leaves of its row. A third tree would be a third window and a third factor, so the update says
// so here rather than looping over a dimension count that cannot change.
static_assert(NUMBER_OF_TREES == 2,
              "The field update is written for one species tree and one molecule tree.");

TMarkovField::TMarkovField(size_t n_iterations, std::vector<std::unique_ptr<TTree>> &Trees,
                           std::string _prefix)
    : _trees(Trees), _prefix(std::move(_prefix)) {
	using namespace coretools::instances;

	// find molecule and species dimensions; construct mass spec data if both trees are present
	// _ms_data.emplace(_trees); // TODO: once we have data, we can remove this

	// read: fix Y or Z?
	_fix_Y = ProgramOptions::FIX_Y;
	if (_fix_Y) { logfile().list("Will fix Y during the MCMC."); }
	_fix_Z = ProgramOptions::FIX_Z;
	if (_fix_Z) { logfile().list("Will fix Z during the MCMC."); }

	// initialize Y: one dimension per tree, sized by that tree's leaf count. The loop is bounded by
	// the array and not by the tree count, because the array is NUMBER_OF_TREES long and nothing
	// stops --tree_others from making more trees than that.
	IndexArray num_leaves_per_dim{};
	for (size_t i = 0; i < num_leaves_per_dim.size(); ++i) {
		num_leaves_per_dim[i] = _trees[i]->get_number_of_leaves();
	}
	_Y.initialize(n_iterations, num_leaves_per_dim);
	if (parameters().exists("set_Y")) {
		std::string filename = parameters().get("set_Y", "acol_simulated_Y.txt");
		_read_Y_from_file(filename);
	}
}

void TDataUpdateAccumulator::commit(TDataModel &data_model) {
#ifdef USE_LOTUS
	double sum_new_LL = 0.0;
	for (auto &i : _lotus_LL) {
		// loop over all LL (stored per thread) and sum
		sum_new_LL += i.getSum();
	}
	data_model.get_lotus().update_cur_LL(sum_new_LL);
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
	size_t total_disagree = 0;
	for (const auto &i : _n_disagree) { total_disagree += i; }
	// The update visits every cell of Y exactly once, so this is the complete disagreement count.
	data_model.get_simple_error_model().set_n_disagree(total_disagree);
#endif
}

void TMarkovField::_read_Y_from_file(const std::string &filename) {
	coretools::TInputFile file(filename, coretools::FileType::Header);
	if (file.numCols() != 5) {
		throw coretools::TUserError("Simulated Y is expected to have 5 columns, but has ",
		                            file.numCols(), " !");
	}

	// read each line of the file
	for (; !file.empty(); file.popFront()) {
		auto linear_index_in_Y_space = file.get<uint64_t>(0);
		bool state                   = file.get<bool>(1);
		if (state) { _Y.insert_one(linear_index_in_Y_space); }
	}
}

//-----------------------------------
// The field update
//-----------------------------------

template<bool IsSimulation, bool initYFromData>
void TMarkovField::_update_field_row(size_t species_leaf, TDataModel &data_model,
                                     const TCellUniforms &uniforms,
                                     TDataUpdateAccumulator &accumulator,
                                     std::vector<size_t> &deferred_inserts) {
	const TTree &species_tree      = *_trees.front();
	const TTree &molecule_tree     = *_trees.back();
	const size_t n_molecule_leaves = molecule_tree.get_number_of_leaves();
	const auto thread              = static_cast<size_t>(omp_get_thread_num());
	const IndexArray row_start{species_leaf, 0};

	// The field's row: every cell this call draws, and the only cells it writes.
	auto field_row = _Y.open_window(row_start, n_molecule_leaves, /*stride=*/1);

	// The species parent's row. A leaf is never a root, so its parent is an internal node and the
	// species tree's node state holds it.
	//
	// The two node states are reached through the mutable trees because a window can write, and
	// opening one is a non-const operation. Neither window below is written: the two-state draw
	// writes the field alone.
	auto species_parent_row = _trees.front()->get_Z().open_window(
	    IndexArray{species_tree.parent_of(species_leaf), 0}, n_molecule_leaves, /*stride=*/1);

	// The molecule node state's row, over every molecule node. A cell's molecule parent is a
	// column of this one row, so one window covers every parent the row reads.
	auto molecule_row = _trees.back()->get_Z().open_window(
	    row_start, molecule_tree.get_number_of_nodes(), /*stride=*/1);

	// Each data source's own row of the same cells. Both have the field's dimensions, so the
	// field's index is already theirs.
#ifdef USE_LOTUS
	auto lotus_row = data_model.get_lotus().open_row(row_start, n_molecule_leaves);
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
	auto simple_data_row =
	    data_model.get_simple_error_model().open_row(row_start, n_molecule_leaves);
#endif

	// The branch from the species leaf to its parent is the same for every cell of the row. So is
	// the molecule tree's clique, which a species leaf names.
	const auto species_branch      = species_tree.get_binned_branch_length(species_leaf);
	const TClique &molecule_clique = molecule_tree.get_clique(row_start);

	for (size_t molecule_leaf = 0; molecule_leaf < n_molecule_leaves; ++molecule_leaf) {
		const IndexArray cell{species_leaf, molecule_leaf};

		// prepare log probabilities for the two possible states
		std::array<coretools::TSumLogProbability, 2> sum_log;

		// calculate probabilities in Markov random field: one factor per tree, each the two-state
		// process on the branch from that tree's parent
		if constexpr (!initYFromData) {
			// The species tree's clique is named by the molecule leaf, so it changes along the row.
			species_tree.get_clique(cell).calculate_log_prob_parent_to_node(
			    species_branch, species_parent_row.is_one(molecule_leaf), sum_log);
			molecule_clique.calculate_log_prob_parent_to_node(
			    molecule_tree.get_binned_branch_length(molecule_leaf),
			    molecule_row.is_one(molecule_tree.parent_of(molecule_leaf)), sum_log);
		}
		// getSum() is not const, so the copy the joint density reads from is not either.
		std::array<coretools::TSumLogProbability, 2> sum_log_field = sum_log;

		// Declared outside the IsSimulation branch so the simulation instantiation does not warn
		// about an unused variable. 1.0 is the neutral value, adding log(1) = 0.
#ifdef USE_LOTUS
		std::array<double, 2> prob_lotus{1.0, 1.0};
#endif
		if constexpr (!IsSimulation) {
			// calculate log likelihood (lotus)
#ifdef USE_LOTUS
			data_model.get_lotus().calculate_LL_update_Y(cell, lotus_row.is_one(molecule_leaf),
			                                             prob_lotus);
			for (size_t i = 0; i < 2; ++i) { sum_log[i].add(prob_lotus[i]); }
#endif
			// calculate log likelihood (simple error model)
#ifdef USE_SIMPLE_ERROR_MODEL
			std::array<double, 2> prob_simple{};
			data_model.get_simple_error_model().probabilities_for_Y_update(
			    simple_data_row.is_one(molecule_leaf), prob_simple);
			for (size_t i = 0; i < 2; ++i) { sum_log[i].add(prob_simple[i]); }
#endif
			// calculate log likelihood mass spec data
			if (_ms_data.has_value()) { _ms_data->add_log_likelihood(cell, sum_log); }
		}

		// sample state. The cell names the uniform it draws, so the thread that happens to reach
		// this cell does not decide what it gets.
		const bool new_state = sample(sum_log, uniforms.at(field_row.linear_index(molecule_leaf)));

		// Write the cell back through the window it was read from. A window that holds the cell
		// writes it in place, and one that does not buffers the insert for the caller.
		if (field_row.is_one(molecule_leaf) != new_state) {
			field_row.set_state(molecule_leaf, new_state);
		}

		TYUpdateResult result;
#ifdef USE_LOTUS
		result.prob_lotus_new_state = prob_lotus[new_state];
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
		if constexpr (!IsSimulation) {
			// The observed cell contradicts the state the field was just given.
			result.simple_model_disagrees = simple_data_row.is_one(molecule_leaf) != new_state;
		}
#endif
		if constexpr (!IsSimulation) { accumulator.add(thread, result); }

		_complete_log_density[thread] += sum_log_field[static_cast<size_t>(new_state)].getSum();
	}

	// The window ends here, inside the parallel region, so it hands its inserts out rather than
	// making them: an insert writes one row and one column of a sparse field, and the threads
	// share every column. See ADR-0006.
	deferred_inserts = field_row.take_buffered_inserts();
}

template<bool IsSimulation, bool initYFromData>
void TMarkovField::_update_all_Y(TDataModel &data_model, size_t iteration) {
	_reset_log_joint_density();

	if (iteration == 0 && ProgramOptions::WRITE_Y_TRACE && !_Y_trace_file.isOpen() && !_fix_Y) {
		std::vector<size_t> Y_trace_header;
		Y_trace_header.reserve(_Y.total_size_of_container_space());
		for (size_t i = 0; i < _Y.total_size_of_container_space(); ++i) {
			Y_trace_header.push_back(i);
		}
		if constexpr (IsSimulation) {
			_Y_trace_file.open(_prefix + "_simulated_Y_trace.txt", Y_trace_header, "\t");
		} else {
			_Y_trace_file.open(_prefix + "_Y_trace.txt", Y_trace_header, "\t");
		}
	}

	if (_fix_Y) {
		// keep the two ifs separate because if Y is not empty, then we just return
		if (_Y.empty()) {
			throw coretools::TUserError("Y is currently empty and fixed. Was Y read from a file ? "
			                            "(--set_Y)");
		}
		return;
	}

	// The stream the field's cells draw from this iteration, built before the parallel region
	// (see run_seed).
	const TCellUniforms field_uniforms(run_seed(), TCellStream::field, iteration);

	TDataUpdateAccumulator accumulator(ProgramOptions::NUMBER_OF_THREADS);

	// One list per species leaf, filled by the window that leaf's row is written through. No two
	// leaves share a list, so the lists need neither a thread index nor a lock, and nothing has to
	// be drained inside the region.
	// Not const: OpenMP made a const variable predetermined shared before version 4.0 and does not
	// now, so a `default(none)` clause that names one is right under some compilers and wrong under
	// others.
	size_t n_species_leaves = _trees.front()->get_number_of_leaves();
	std::vector<std::vector<size_t>> deferred_inserts(n_species_leaves);

	// A thread takes a species leaf. Two rows of the field share no cell, and a field cell's
	// Markov blanket holds no other field cell, so the rows are conditionally independent given
	// the node states. A cell's uniform is hashed from its position (ADR-0007), so which thread
	// reaches a cell, and in which order, moves no number.
	//
	// The species tree's leaf count is therefore the widest this team can run. Splitting the other
	// way would give a thread a column, and a column of the sparse field is what no window may
	// insert into (ADR-0006).
#pragma omp parallel for num_threads(ProgramOptions::NUMBER_OF_THREADS)                            \
    schedule(static) default(none)                                                                 \
    shared(accumulator, deferred_inserts, data_model, field_uniforms, n_species_leaves)
	for (size_t species_leaf = 0; species_leaf < n_species_leaves; ++species_leaf) {
		_update_field_row<IsSimulation, initYFromData>(species_leaf, data_model, field_uniforms,
		                                               accumulator, deferred_inserts[species_leaf]);
	}

	// The inserts every row deferred, in one batch. The dense field is handed empty lists and does
	// nothing with them.
	_Y.insert_in_Y(deferred_inserts);
	// at the very end: sum the per-thread accumulators and store them in the data sources
	if constexpr (!IsSimulation) { accumulator.commit(data_model); }
	if (ProgramOptions::WRITE_Y_TRACE && (iteration % _Y.get_thinning_factor() == 0) && !_fix_Y) {
		_Y_trace_file.writeln(_Y.get_full_Y_binary_vector());
	}
}

void TMarkovField::update(TDataModel &data_model, size_t iteration) {
	if (ProgramOptions::WRITE_JOINT_LOG_PROB_DENSITY && iteration == 0 &&
	    !_joint_density_file.isOpen()) {
		_joint_density_file.open(_prefix + "_simulated_joint_density.txt",
		                         {
		                             "joint_density",
		                         },
		                         "\t");
	}

	if (iteration == 0 && !_z_initialized_from_children) {
		_update_all_Y<false, true>(data_model, iteration);
		for (auto &tree : _trees) { tree->initialize_Z_from_children(_Y); }
		_z_initialized_from_children = true;
	} else {
		_update_all_Y<false, false>(data_model, iteration);
	}
	if (_fix_Z) {
		_update_all_Z<false, true>(iteration);
	} else {
		_update_all_Z<false, false>(iteration);
	}
	if (_ms_data.has_value()) _ms_data->update_all_MS_assignments();
	_Y.add_to_counter(iteration);
	// calculate joint density
	if (ProgramOptions::WRITE_JOINT_LOG_PROB_DENSITY && iteration % _Y.get_thinning_factor() == 0) {
		auto sum_log_field = _calculate_complete_joint_density();
		_joint_density_file.writeln(sum_log_field);
	}
}

void TMarkovField::simulate(TDataModel &data_model) {
	// For simulation we always draw from the prior. (top-down)
	// 1. Draw branch len -> draw mus
	// 2. For every tree draw the root from those mus and the we BFS sample all the internal nodes
	// based on the state of the parent. At the end also draw the state of the leaves. For each tree
	// we know which leaves and which internal node are 0 or 1.
	// 2. draw Z since its unique per tree
	// 3. Draw Y (this function)
	// 4. Draw the data of every compiled-in source from that Y (done by the caller,
	//    TDataModel::_simulateUnderPrior, so that all sources see the same Y).
	// 5. From that simulated data can we go back and infer the prior params
	//
	for (size_t tree_index = 0; tree_index < _trees.size(); ++tree_index) {
		auto &tree = _trees[tree_index];
		tree->simulate_Z(tree_index);
	}

	_simulate_Y();

	// for iteration in 1->max_iteration, (max_iteration should be passed from CLI)
	// we use tree.update_Z(); and then
	// update Y where likelihood of data is always one so it doesn't matter.
	size_t max_iteration = get_num_iterations_simulation();

	// create the Markov field density file
	if (ProgramOptions::WRITE_JOINT_LOG_PROB_DENSITY) {
		_joint_density_file.open(_prefix + "_simulated_joint_density.txt",
		                         {
		                             "joint_density",
		                         },
		                         "\t");
	}

	std::string report =
	    "Running an MCMC chain of " + coretools::str::toString(max_iteration) + " iterations";
	coretools::TProgressReporter prog(max_iteration, report);
	for (size_t iteration = 0; iteration < max_iteration; ++iteration) {

		_update_all_Y<true, false>(data_model, iteration);

		if (_fix_Z) {
			_update_all_Z<true, true>(iteration);
		} else {
			_update_all_Z<true, false>(iteration);
		}
		_Y.add_to_counter(iteration);

		// calculate joint density
		if (iteration % _Y.get_thinning_factor() == 0 &&
		    ProgramOptions::WRITE_JOINT_LOG_PROB_DENSITY) {
			auto sum_log_field = _calculate_complete_joint_density();
			_joint_density_file.writeln(sum_log_field);
		}
		prog.next();
	}
	prog.done();
	if (ProgramOptions::WRITE_Y) { _write_Y_to_file<true>(_prefix + "_simulated_Y.txt"); }
	if (ProgramOptions::WRITE_Z) {
		for (const auto &tree : _trees) {
			write_Z_to_file(_prefix + "_simulated_Z_" + tree->get_tree_name() + ".txt", *tree,
			                _trees, /*write_full_Z =*/true);
			if (ProgramOptions::WRITE_BRANCH_LENGTHS) { write_branch_length_grid(*tree); }
		}
	}
}

void TMarkovField::_simulate_Y() {
	// to sample Y we need to know the state of the parent for each leaf that is represented in a Y
	// entry. We are going to iterate over all possible Y and sample givent the product of
	// probabilities of the child given the parent. set number of leaves per dimension (set the last
	// dimension to one)
	if (ProgramOptions::SIMULATION_NO_Y_INITIALIZATION) { return; }
	// A stream of its own, so this draw and the chain's first update are two draws (ADR-0007).
	const TCellUniforms uniforms(run_seed(), TCellStream::field_at_start, 0);
	for (size_t linear_index_in_leaves_space = 0;
	     linear_index_in_leaves_space < _Y.total_size_of_container_space();
	     ++linear_index_in_leaves_space) {
		auto multidim_index_in_Y = _Y.get_multi_dimensional_index(linear_index_in_leaves_space);
		std::array<coretools::TSumLogProbability, 2> sum_log;
		for (size_t dim = 0; dim < _trees.size(); ++dim) {
			// get relevant clique
			const auto &clique = _trees[dim]->get_clique(multidim_index_in_Y);
			// translate index in leaves to the index in tree
			const size_t index_in_tree =
			    _trees[dim]->get_node_index_from_leaf_index(multidim_index_in_Y[dim]);
			// calculate P(parent | node = 0) and P(parent | node = 1)
			// Note: leaves can never be roots -> they always have a parent (no need to bother with
			// stationary)
			// The parent reads as 0, which is what this simulator has always drawn from: it built
			// a current-state cache per cell and never filled it. Reading the node state instead
			// would change every simulated field, so #42 repairs it with the rest of the
			// simulator.
			constexpr bool state_of_parent = false;
			clique.calculate_log_prob_parent_to_node(
			    _trees[dim]->get_binned_branch_length(index_in_tree), state_of_parent, sum_log);
		}
		bool y_state = sample(sum_log, uniforms.at(linear_index_in_leaves_space));
		if (y_state) { _Y.insert_one(linear_index_in_leaves_space); }
	}
}

void TMarkovField::burninHasFinished() {
	_Y.reset_counts();
	_Y.remove_zeros();
}

void TMarkovField::oneBurninHasFinished() { _Y.remove_zeros(); }

void TMarkovField::MCMCHasFinished() {
	// write function to write the posterior state of Y to file
	_write_Y_to_file<false>(_prefix + "_Y_posterior.txt");
}

const TFieldStorage &TMarkovField::get_Y_matrix() const { return _Y; }

double TMarkovField::_calculate_complete_joint_density() {

	// we can initialize the sum_log_field for the joint probability of the Markov random field

	// Easy case: Y
	double sum_log_field = coretools::containerSum(_complete_log_density);

	// now we loop over all Z to get the joint probability
	for (const auto &tree : _trees) { sum_log_field += tree->get_complete_joint_density(); }

	return sum_log_field;
};
