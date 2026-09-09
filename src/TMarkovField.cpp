//
// Created by madleina on 10.01.25.
//

#include "TMarkovField.h"
#include "TDataModel.h"
#include "Types.h"
#include "cli.h"
#include "constants.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/progressTools.h"
#include "coretools/algorithms.h"
#include "field/TBlockModel.h"
#include "field/TBlockUpdate.h"
#include "field/leaf_layer_start.h"
#include "field/link_backend.h"
#include "field/simulate_field.h"
#include "random/TCellUniforms.h"
#include "storages/storage_backend.h"
#include "tree/TTree.h"
#include "tree/io/node_state_columns.h"
#include "tree/io/write_Z.h"
#include "tree/io/write_tree_field.h"
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

TMarkovField::TMarkovField(size_t n_iterations, std::vector<std::unique_ptr<TTree>> &Trees,
                           TypeParamErrorProbability *omega, std::string _prefix, bool simulate)
    : _trees(Trees), _prefix(std::move(_prefix)), _simulate(simulate), _omega(omega) {
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

	// What the posterior field is worth, said before the chain runs. A field cell holds its
	// posterior counter in 15 bits, so a chain longer than 32767 iterations is thinned to fit, and
	// the factor is the resolution of the posterior field, of both tree field posteriors, and of
	// every trace the run writes. Both backends hold the same cell and thin the same chain
	// identically. A dense run written before that was true counted one iteration in
	// ceil(n / 65535), so its posterior field is at twice the resolution of this one.
	logfile().list("The posterior field counts one iteration in ", _Y.get_thinning_factor(),
	               ", which is what a 15-bit counter holds over ", n_iterations, " iterations.");

	// The block update is written for one species tree and one molecule tree, and a cell index
	// holds two coordinates. A third tree would run past the end of both. The loop above is
	// bounded by the array so that it reaches this line, which says so.
	if (_trees.size() != NUMBER_OF_TREES) {
		throw coretools::TUserError("The model is written for ", NUMBER_OF_TREES,
		                            " trees, but the run has ", _trees.size(), ".");
	}

	// One posterior per tree, over the leaf-pair space, thinned the way the field thins its own.
	// Both files then count the same iterations and their denominators agree.
	//
	// Only when the run asks. A counter per leaf pair per tree is what a run that chose the sparse
	// field chose not to pay for the field itself (ADR-0006), so an empty list is what the counting
	// and the writing below read as "not this run".
	if (ProgramOptions::WRITE_TREE_FIELD_POSTERIORS) {
		_tree_field_posteriors.resize(_trees.size());
		for (auto &posterior : _tree_field_posteriors) {
			posterior.initialize(_Y.total_size_of_container_space(), _Y.get_thinning_factor());
		}
	}

	if (parameters().exists("set_Y")) {
		std::string filename = parameters().get("set_Y", "acol_simulated_Y.txt");
		_read_Y_from_file(filename);
		_field_came_from_a_file = true;
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

void TMarkovField::_open_Y_trace_file() {
	std::vector<size_t> Y_trace_header;
	Y_trace_header.reserve(_Y.total_size_of_container_space());
	for (size_t i = 0; i < _Y.total_size_of_container_space(); ++i) { Y_trace_header.push_back(i); }
	_Y_trace_file.open(_prefix + "_Y_trace.txt", Y_trace_header, "\t");
}

void TMarkovField::_trace_link_counters(size_t iteration) {
	if (!_link_counters_file.isOpen()) {
		std::vector<std::string> header;
		header.reserve(2 * field_math::TLinkCounters::n_buckets);
		for (size_t bucket = 0; bucket < field_math::TLinkCounters::n_buckets; ++bucket) {
			header.push_back("n_bucket" + std::to_string(bucket) + "_field0");
			header.push_back("n_bucket" + std::to_string(bucket) + "_field1");
		}
		const std::string suffix =
		    _simulate ? "_simulated_link_counters_trace.txt" : "_link_counters_trace.txt";
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
	               " counted cells, pooled over every tally the trace holds. They carry that "
	               "run's noise. A residual that survives a longer chain means the link is wrong. "
	               "That is a finding, and it fails nothing.");
	logfile().endIndent();
}

/// A held field holds both tree fields with it, because the block draws all three together.
///
/// The tree fields take the field's own states. Under the AND link that is the configuration the
/// field is most likely to have come from, and it is what the internal nodes are built from
/// (TTree::initialize_Z_from_children). Before the block, that walk read the field itself.
void TMarkovField::_hold_tree_fields_at_the_field() {
	_link_counters = leaf_layer_start::hold_tree_fields_at_the_field<TLinkPolicy>(
	    _Y, _trees.front()->get_Z(), _trees.back()->get_Z());
}

void TMarkovField::_throw_if_the_fixed_field_is_empty() const {
	if (!_fix_Y || !_Y.empty()) { return; }
	throw coretools::TUserError("Y is currently empty and fixed. Was Y read from a file ? "
	                            "(--set_Y)");
}

/// The chain start, which CONTEXT.md names and ADR-0005 argues.
///
/// The field takes the LOTUS records. Both tree fields take the field. Each tree then initialises
/// every internal node from its children. With no record anywhere the start is all zeros.
///
/// A field the run gave stands as it is. --set_Y is the start the run asked for, and --fix_Y holds
/// one for the whole chain. --fix_Y without a file leaves the field empty, which this reports as
/// the user error it is.
///
/// The six counters this leaves are degenerate. Bucket 1 holds nothing, so the AND diagnostic says
/// nothing about them. That is accepted. The first block update recounts every leaf pair and
/// replaces the tally. A fixed field runs no block update, so its degenerate tally stands for the
/// whole chain -- which is the truth of a field both tree fields match exactly.
void TMarkovField::_start_the_chain([[maybe_unused]] const TDataModel &data_model) {
	using namespace coretools::instances;

	// Before the log line below reports a start, and before every clique is walked for nothing.
	_throw_if_the_fixed_field_is_empty();

	// A fixed field is a field from a file, so the flag alone decides.
	if (!_field_came_from_a_file) {
#ifdef USE_LOTUS
		leaf_layer_start::start_the_field_at(data_model.get_lotus().get_L(), _Y);
#endif
	}
	_hold_tree_fields_at_the_field();
	logfile().list("The chain starts with the field and both tree fields at one in ",
	               _Y.number_of_ones(), " of ", _Y.total_size_of_container_space(), " cells.");

	for (auto &tree : _trees) { tree->initialize_Z_from_children(); }
}

void TMarkovField::_update_block(TDataModel &data_model, size_t iteration) {
	if (iteration == 0 && ProgramOptions::WRITE_Y_TRACE && !_Y_trace_file.isOpen() && !_fix_Y) {
		_open_Y_trace_file();
	}

	if (_fix_Y) {
		_throw_if_the_fixed_field_is_empty();
		// The block draws the two tree fields with the field, so holding one holds all three. The
		// leaf layer never moves after this, and one tally stands for the whole chain; the chain's
		// start is where it was built (_start_the_chain).
		//
		// The error probability still moves against that tally, so the trace still carries it.
		_trace_link_counters(iteration);
		return;
	}

	// The stream the leaf pairs draw from this iteration, built before the parallel region
	// (see run_seed).
	const TCellUniforms field_uniforms(run_seed(), TCellStream::field, iteration);

	TDataUpdateAccumulator accumulator(ProgramOptions::NUMBER_OF_THREADS);
	TBlockModel model(_trees, data_model, accumulator);
	std::vector<block_update::TThreadTally> tallies(ProgramOptions::NUMBER_OF_THREADS);

	const field_math::TErrorProbability omega = _error_probability();
	block_update::run<TLinkPolicy>(_Y, _trees.front()->get_Z(), _trees.back()->get_Z(),
	                               _trees.front()->phylogeny(), _trees.back()->phylogeny(), omega,
	                               model, field_uniforms, tallies);

	// The update covered every leaf pair, so the tally it built is the configuration itself.
	_link_counters = field_math::TLinkCounters();
	for (const auto &tally : tallies) { _link_counters.merge(tally.counters); }

	_trace_link_counters(iteration);

	// at the very end: sum the per-thread accumulators and store them in the data sources
	accumulator.commit(data_model);
	if (ProgramOptions::WRITE_Y_TRACE && (iteration % _Y.get_thinning_factor() == 0) && !_fix_Y) {
		_Y_trace_file.writeln(_Y.get_full_Y_binary_vector());
	}
}

void TMarkovField::update(TDataModel &data_model, size_t iteration) {
	// The chain is started before its first update, so that update reads a node state both trees
	// have something to say about.
	if (iteration == 0 && !_chain_started) {
		_start_the_chain(data_model);
		_chain_started = true;
	}

	// The block update is the whole of the leaf layer's turn: it draws the field and both tree
	// fields together. Then each tree walks its own internal nodes, and then the parameters move.
	_update_block(data_model, iteration);
	if (_fix_Z) {
		_update_all_Z<true>(iteration);
	} else {
		_update_all_Z<false>(iteration);
	}
	if (_ms_data.has_value()) _ms_data->update_all_MS_assignments();
	_Y.add_to_counter(iteration);
	_count_the_tree_fields(iteration);

	// Last of all, because it is the density of the configuration this iteration leaves behind.
	_trace_joint_density(iteration, data_model);
}

void TMarkovField::simulate(TDataModel &data_model) {
	using namespace coretools::instances;

	// Simulation is a forward draw, and nothing more (ADR-0005):
	//
	//   1. each tree draws its whole node state top-down, leaves included, so each tree ends up
	//      holding its own tree field;
	//   2. the link turns the two tree fields into the field;
	//   3. the caller derives every compiled-in data source from that one field, so that all of
	//      them see the same one (TDataModel::_simulateUnderPrior).
	//
	// No chain runs here. Every factor is a proper conditional density, so the draw is exact and
	// there is nothing left for a chain to repair. ADR-0005 says why the model this replaces
	// needed one.

	// A forward draw leaves one configuration, and the files below hold all of it. Saying so is
	// cheaper than letting a run ask for a per-iteration trace and find none.
	if (ProgramOptions::WRITE_Y_TRACE || ProgramOptions::WRITE_Z_TRACE) {
		logfile().list("A simulated configuration is one draw, so there is no trace of it. The "
		               "field and both node states are written in full instead.");
	}

	for (auto &tree : _trees) { tree->simulate_Z(); }

	_simulate_Y();

	// One draw, counted once, so the field file's fraction column reads the state beside it.
	_Y.add_to_counter(0);
	_trace_link_counters(0);
	_trace_joint_density(0, data_model);
	// The link's two parameter-free constraints, read off the field just drawn. A simulated field
	// satisfies both up to its own sampling noise, so a residual here is a defect in this simulator
	// rather than a finding about data (ADR-0005).
	_report_link_diagnostic();

	if (ProgramOptions::WRITE_Y) { _write_Y_to_file<true>(_prefix + "_simulated_Y.txt"); }
	if (ProgramOptions::WRITE_Z) {
		for (const auto &tree : _trees) {
			write_Z_to_file(_prefix + "_simulated_Z_" + tree->get_tree_name() + ".txt",
			                tree->get_Z(), node_state_columns(_trees), /*write_full_Z =*/true);
			if (ProgramOptions::WRITE_BRANCH_LENGTHS) { write_branch_length_grid(*tree); }
		}
	}
}

void TMarkovField::_simulate_Y() {
	// A stream of its own, so this draw and the chain's first update are two draws (ADR-0007).
	const TCellUniforms uniforms(run_seed(), TCellStream::field_at_start, 0);

	// The field is a noisy AND of the two tree fields, cell by cell. The draw hands back the six
	// counters over it.
	_link_counters = simulate_field::draw_from_the_tree_fields<TLinkPolicy>(
	    _Y, _trees.front()->get_Z(), _trees.back()->get_Z(), _error_probability(), uniforms);
}

void TMarkovField::burninHasFinished() {
	_Y.reset_counts();
	_Y.remove_zeros();
	// Every posterior reports the chain, not the burn-in that preceded it.
	for (auto &posterior : _tree_field_posteriors) { posterior.reset_counts(); }
	// The diagnostic reports the chain, not the burn-in that preceded it.
	_traced_link_counters = field_math::TLinkCounters();
}

void TMarkovField::oneBurninHasFinished() { _Y.remove_zeros(); }

void TMarkovField::MCMCHasFinished() {
	// write function to write the posterior state of Y to file
	_write_Y_to_file<false>(_prefix + "_Y_posterior.txt");
	// Each tree field's posterior stands beside it, in a file of its own. The field's own says
	// nothing about how the two trees split the rate between them (ADR-0005, derivation 3).
	_write_tree_field_posteriors();
	_report_link_diagnostic();
}

const TFieldStorage &TMarkovField::get_Y_matrix() const { return _Y; }

//-----------------------------------
// The joint density, and the tree fields' posteriors
//-----------------------------------

joint_density::TJointDensity TMarkovField::_calculate_joint_density(const TDataModel &data_model) {
	// One column per tree, and the constructor has already said there are NUMBER_OF_TREES of them.
	joint_density::TJointDensity density;

	// Each tree's own factor, over its whole node state, leaves included. One term per node, so
	// each branch is counted once and each factor is a density (tree/node_state_density.h).
	for (size_t tree_idx = 0; tree_idx < _trees.size(); ++tree_idx) {
		density.node_states[tree_idx] = _trees[tree_idx]->log_node_state_density();
	}

	// The link, for the whole field at once. Its likelihood is a function of the six counters and
	// the error probability alone, so this is six integers rather than a walk over the cells
	// (ADR-0005).
	density.link = _link_log_likelihood();

	// A simulated chain draws from the prior. Every data term is neutral in the block update, and
	// no source has scored the field, so there is nothing here to add.
	if (!_simulate) { density.data = data_model.data_log_likelihood(); }

	return density;
}

void TMarkovField::_trace_joint_density(size_t iteration, const TDataModel &data_model) {
	if (!ProgramOptions::WRITE_JOINT_LOG_PROB_DENSITY) { return; }
	if (iteration % _Y.get_thinning_factor() != 0) { return; }

	if (!_joint_density_file.isOpen()) {
		std::vector<std::string> tree_names;
		tree_names.reserve(_trees.size());
		for (const auto &tree : _trees) { tree_names.push_back(tree->get_tree_name()); }
		const std::string suffix =
		    _simulate ? "_simulated_joint_density.txt" : "_joint_density.txt";
		_joint_density_file.open(_prefix + suffix, joint_density::trace_header(tree_names), "\t");
	}

	_joint_density_file.writeln(joint_density::trace_row(_calculate_joint_density(data_model)));
}

void TMarkovField::_count_the_tree_fields(size_t iteration) {
	// Empty unless the run asked for the files, and then there is nothing to count.
	for (size_t tree_idx = 0; tree_idx < _tree_field_posteriors.size(); ++tree_idx) {
		_tree_field_posteriors[tree_idx].add_to_counter(iteration, _Y, _trees[tree_idx]->get_Z());
	}
}

void TMarkovField::_write_tree_field_posteriors() const {
	if (_tree_field_posteriors.empty()) { return; }
	const auto columns = node_state_columns(_trees);
	for (size_t tree_idx = 0; tree_idx < _tree_field_posteriors.size(); ++tree_idx) {
		write_tree_field_posterior(
		    _prefix + "_" + _trees[tree_idx]->get_tree_name() + "_tree_field_posterior.txt", _Y,
		    _trees[tree_idx]->get_Z(), _tree_field_posteriors[tree_idx], columns);
	}
}
