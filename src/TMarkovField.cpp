//
// Created by madleina on 10.01.25.
//

#include "TMarkovField.h"
#include "TClique.h"
#include "TDataModel.h"
#include "Types.h"
#include "cli.h"
#include "constants.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/progressTools.h"
#include "coretools/algorithms.h"
#include "field/TBlockUpdate.h"
#include "field/TBlockModel.h"
#include "field/link_backend.h"
#include "random/TCellUniforms.h"
#include "storages/storage_backend.h"
#include "tree/TTree.h"
#include "tree/io/write_Z.h"
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

TMarkovField::TMarkovField(size_t n_iterations, std::vector<std::unique_ptr<TTree>> &Trees,
                           TypeParamErrorProbability *omega, std::string _prefix)
    : _trees(Trees), _prefix(std::move(_prefix)), _omega(omega) {
	using namespace coretools::instances;

	// find molecule and species dimensions; construct mass spec data if both trees are present
	// _ms_data.emplace(_trees); // TODO: once we have data, we can remove this

	// read: fix Y or Z?
	_fix_Y = ProgramOptions::FIX_Y;
	if (_fix_Y) {
		// The block draws the field and both tree fields together, so holding one holds all three.
		logfile().list("Will fix Y, and with it both tree fields, during the MCMC.");
	}
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
// The error probability
//-----------------------------------

void TMarkovField::set_error_probability_support() {
	TypeErrorProbability::setMin(std::numeric_limits<double>::min());
	TypeErrorProbability::setMax(std::nextafter(0.5, 0.0));
}

field_math::TErrorProbability TMarkovField::_error_probability() const {
	// The type keeps the value inside (0, 0.5), and TErrorProbability checks it again. The second
	// check costs nothing here: this is called once per update, not once per cell.
	return field_math::TErrorProbability(static_cast<double>(_omega->value()));
}

double TMarkovField::_link_log_likelihood() const {
	return TLinkPolicy::log_likelihood(_link_counters, _error_probability());
}

double TMarkovField::link_log_likelihood_ratio() const {
	const field_math::TErrorProbability old_omega(static_cast<double>(_omega->oldValue()));
	return TLinkPolicy::log_likelihood_ratio(_link_counters, old_omega, _error_probability());
}

//-----------------------------------
// The block update
//-----------------------------------

void TMarkovField::_open_Y_trace_file(bool is_simulation) {
	std::vector<size_t> Y_trace_header;
	Y_trace_header.reserve(_Y.total_size_of_container_space());
	for (size_t i = 0; i < _Y.total_size_of_container_space(); ++i) {
		Y_trace_header.push_back(i);
	}
	const std::string suffix = is_simulation ? "_simulated_Y_trace.txt" : "_Y_trace.txt";
	_Y_trace_file.open(_prefix + suffix, Y_trace_header, "\t");
}

void TMarkovField::_trace_link_counters(size_t iteration, bool is_simulation) {
	if (!_link_counters_file.isOpen()) {
		std::vector<std::string> header;
		header.reserve(2 * field_math::TLinkCounters::n_buckets);
		for (size_t bucket = 0; bucket < field_math::TLinkCounters::n_buckets; ++bucket) {
			header.push_back("n_bucket" + std::to_string(bucket) + "_field0");
			header.push_back("n_bucket" + std::to_string(bucket) + "_field1");
		}
		const std::string suffix =
		    is_simulation ? "_simulated_link_counters_trace.txt" : "_link_counters_trace.txt";
		_link_counters_file.open(_prefix + suffix, header, "\t");
	}

	if (iteration % _Y.get_thinning_factor() != 0) { return; }

	std::vector<size_t> line;
	line.reserve(2 * field_math::TLinkCounters::n_buckets);
	for (size_t bucket = 0; bucket < field_math::TLinkCounters::n_buckets; ++bucket) {
		line.push_back(_link_counters.count(bucket, false));
		line.push_back(_link_counters.count(bucket, true));
	}
	_link_counters_file.writeln(line);

	// The diagnostic reads what this file wrote, and nothing else.
	_traced_link_counters.merge(_link_counters);
}

void TMarkovField::_report_link_diagnostic() const {
	using namespace coretools::instances;

	const auto diagnostic = TLinkPolicy::diagnose(_traced_link_counters);
	logfile().startIndent("The AND link, checked against the counters (ADR-0005):");
	if (!diagnostic.is_complete()) {
		logfile().list("A bucket held no cell, so neither constraint says anything.");
		logfile().endIndent();
		return;
	}
	logfile().list("P_0 = ", diagnostic.prob[0], ", P_1 = ", diagnostic.prob[1],
	               ", P_2 = ", diagnostic.prob[2], ".");
	logfile().list("P_1^2 - P_0 * P_2 = ", diagnostic.and_identity_residual,
	               ". It is 0 when the link corrupts the two tree fields independently and ANDs "
	               "them.");
	logfile().list("sqrt(P_0) + sqrt(P_2) - 1 = ", diagnostic.shared_error_probability_residual,
	               ". It is 0 when both trees share one error probability.");
	logfile().list("Both read ", _traced_link_counters.total(),
	               " counted cells, pooled over the chain. They carry that chain's noise. A "
	               "residual that survives a longer chain means the link is wrong. That is a "
	               "finding, and it fails nothing.");
	logfile().endIndent();
}

/// A held field holds both tree fields with it, because the block draws all three together.
///
/// The tree fields take the field's own states. Under the AND link that is the configuration the
/// field is most likely to have come from, and it is what the internal nodes are built from
/// (TTree::initialize_Z_from_children). Before the block, that walk read the field itself.
void TMarkovField::_hold_tree_fields_at_the_field() {
	const TTree &species_tree      = *_trees.front();
	const TTree &molecule_tree     = *_trees.back();
	auto &species_field            = _trees.front()->get_Z();
	auto &molecule_field           = _trees.back()->get_Z();
	const size_t n_molecule_leaves = molecule_tree.get_number_of_leaves();

	_link_counters = field_math::TLinkCounters();
	for (size_t species_leaf = 0; species_leaf < species_tree.get_number_of_leaves();
	     ++species_leaf) {
		for (size_t molecule_leaf = 0; molecule_leaf < n_molecule_leaves; ++molecule_leaf) {
			const IndexArray cell{species_leaf, molecule_leaf};
			const bool y = _Y.is_one(_Y.get_linear_index_in_container_space(cell));
			// Only a cell that disagrees is written. A sparse node state does not hold every cell,
			// and inserting one that already reads as zero would grow it for nothing.
			for (auto *tree_field : {&species_field, &molecule_field}) {
				const size_t index = tree_field->get_linear_index_in_container_space(cell);
				if (tree_field->is_one(index) == y) { continue; }
				if (y) {
					tree_field->insert_one(index);
				} else {
					tree_field->insert_zero(index);
				}
			}
			_link_counters.add(TLinkPolicy::bucket(y, y), y);
		}
	}
}

template<bool IsSimulation, bool InitYFromData>
void TMarkovField::_update_block(TDataModel &data_model, size_t iteration) {
	if (iteration == 0 && ProgramOptions::WRITE_Y_TRACE && !_Y_trace_file.isOpen() && !_fix_Y) {
		_open_Y_trace_file(IsSimulation);
	}

	if (_fix_Y) {
		// keep the two ifs separate because if Y is not empty, then we just return
		if (_Y.empty()) {
			throw coretools::TUserError("Y is currently empty and fixed. Was Y read from a file ? "
			                            "(--set_Y)");
		}
		// The block draws the two tree fields with the field, so holding one holds all three. The
		// leaf layer never moves after this, and one tally stands for the whole chain.
		if (iteration == 0) { _hold_tree_fields_at_the_field(); }
		// The error probability still moves against that tally, so the trace still carries it.
		_trace_link_counters(iteration, IsSimulation);
		return;
	}

	// The stream the leaf pairs draw from this iteration, built before the parallel region
	// (see run_seed).
	const TCellUniforms field_uniforms(run_seed(), TCellStream::field, iteration);

	TDataUpdateAccumulator accumulator(ProgramOptions::NUMBER_OF_THREADS);
	TBlockModel<IsSimulation, InitYFromData> model(_trees, data_model, accumulator);
	std::vector<block_update::TThreadTally> tallies(ProgramOptions::NUMBER_OF_THREADS);

	const field_math::TErrorProbability omega = _error_probability();
	block_update::run<TLinkPolicy>(_Y, _trees.front()->get_Z(), _trees.back()->get_Z(),
	                               _trees.front()->phylogeny(), _trees.back()->phylogeny(), omega,
	                               model, field_uniforms, tallies);

	// The update covered every leaf pair, so the tally it built is the configuration itself.
	_link_counters     = field_math::TLinkCounters();
	_block_log_density = 0.0;
	for (const auto &tally : tallies) {
		_link_counters.merge(tally.counters);
		_block_log_density += tally.log_density;
	}

	_trace_link_counters(iteration, IsSimulation);

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

	// The block update is the whole of the leaf layer's turn: it draws the field and both tree
	// fields together. Then each tree walks its own internal nodes, and then the parameters move.
	if (iteration == 0 && !_z_initialized_from_children) {
		_update_block<false, true>(data_model, iteration);
		for (auto &tree : _trees) { tree->initialize_Z_from_children(); }
		_z_initialized_from_children = true;
	} else {
		_update_block<false, false>(data_model, iteration);
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
	for (auto &tree : _trees) { tree->simulate_Z(); }

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

		_update_block<true, false>(data_model, iteration);

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
	// The diagnostic reports the chain, not the burn-in that preceded it.
	_traced_link_counters = field_math::TLinkCounters();
}

void TMarkovField::oneBurninHasFinished() { _Y.remove_zeros(); }

void TMarkovField::MCMCHasFinished() {
	// write function to write the posterior state of Y to file
	_write_Y_to_file<false>(_prefix + "_Y_posterior.txt");
	_report_link_diagnostic();
}

const TFieldStorage &TMarkovField::get_Y_matrix() const { return _Y; }

double TMarkovField::_calculate_complete_joint_density() {
	// The leaf layer: the two tree factors of every leaf pair, at the states the block drew. The
	// fixed-field path leaves this at zero, because no block ran to compute it.
	double sum_log_field = _block_log_density;

	// The link, for the whole field at once. Its likelihood is a function of the six counters and
	// the error probability alone, so this is six integers rather than a walk over the cells
	// (ADR-0005).
	sum_log_field += _link_log_likelihood();

	// now we loop over all Z to get the joint probability
	for (const auto &tree : _trees) { sum_log_field += tree->get_complete_joint_density(); }

	// This is not yet the ADR-0005 factorisation. A tree's own term scores every node against its
	// parent *and* against each child, so each edge is counted twice. That predates the block, and
	// #41 rebuilds this trace.
	return sum_log_field;
};
