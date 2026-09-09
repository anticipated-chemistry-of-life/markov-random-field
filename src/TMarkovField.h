//
// Created by madleina on 10.01.25.
//

#ifndef ACOL_TMARKOVFIELD_H
#define ACOL_TMARKOVFIELD_H

#include "Types.h"
#include "cli.h"
#include "constants.h"
#include "coretools/Files/TOutputFile.h"
#include "coretools/Main/TError.h"
#include "field/TFieldMath.h"
#include "field/joint_density.h"
#include "field/link_backend.h"
#include "field/tree_field_link.h"
#include "field/tree_field_posterior.h"
#include "mass_spec/msms_data.h"
#include "random/TCellUniforms.h"
#include "storages/storage_backend.h"
#include "tree/TTree.h"
#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

class TDataModel; // forward declaration

//-----------------------------------
// TMarkovField
//-----------------------------------

class TMarkovField {
public:
	/// The error probability, as stattools moves it. It hangs off TDataModel, which owns this
	/// class and forwards the MCMC callbacks, the way the simple error model's rate does.
	using TypeParamErrorProbability = stattools::TParameter<SpecErrorProbability, TDataModel>;

private:
	// trees and Y
	std::vector<std::unique_ptr<TTree>> &_trees;
	TFieldStorage _Y;
	std::string _prefix;

	/// Whether this run simulates. It decides what every file this class writes is called, and
	/// nothing else: a simulated configuration and an inferred chain go to different names so that
	/// one directory can hold both.
	bool _simulate = false;

	// fix values?
	bool _fix_Y = false;
	bool _fix_Z = false;

	// Mass spectrometry data, still dormant. Nothing builds it. The field update does not read it
	// either: it takes the LOTUS and the simple-error term, and adapting a third source is that
	// source's own work.
	std::optional<TMSMSData> _ms_data;

	// The error probability standing between the two tree fields and the field. stattools owns the
	// value and moves it; this is where the field reads it.
	TypeParamErrorProbability *_omega = nullptr;

	// The link's sufficient statistic over the whole field, n(bucket, field state). The field
	// update tallies it as it goes and commits it here, which is what makes the error
	// probability's likelihood O(1) in the number of cells (ADR-0005).
	field_math::TLinkCounters _link_counters;

	// Every counter tally the trace file has written, added together. The AND diagnostic is
	// reported from these, so it pools the cells of the whole chain instead of reading one
	// iteration. Burn-in clears them.
	field_math::TLinkCounters _traced_link_counters;

	/// Whether the chain has been started. The start runs on the first update and not in the
	/// constructor. Initialising the internal nodes reads each clique's transition grid, and the
	/// parameters build those.
	bool _chain_started = false;

	/// Whether --set_Y gave the field its states. The chain leaves such a field as it is, and does
	/// not start it at the LOTUS records.
	bool _field_came_from_a_file = false;

	// The posterior of each tree field, one per tree, over the leaf-pair space the tree fields
	// share with the field. The field's own posterior lives inside its storage; a node state
	// carries no counter, so these stand beside them (field/tree_field_posterior.h).
	std::vector<TTreeFieldPosterior> _tree_field_posteriors;

	// output files
	coretools::TOutputFile _Y_trace_file;
	std::vector<coretools::TOutputFile> _Z_trace_files;
	coretools::TOutputFile _joint_density_file;
	coretools::TOutputFile _link_counters_file;

	/// One field update: every cell of the field's container space, split over the threads.
	/// Defined in TMarkovField.cpp, where the model it hands the pass is complete.
	///
	/// It is the last state update of an iteration, so the six counters it leaves describe the
	/// configuration the error probability then proposes against. A held field draws nothing and is
	/// retallied instead, because the two tree fields moved underneath it.
	///
	/// An inferred chain is the only caller. A simulated one draws its whole configuration forward
	/// and runs no update at all (`simulate`).
	void _update_the_field(TDataModel &data_model, size_t iteration);

	/// Opens the field's trace file on the first iteration of a chain.
	void _open_Y_trace_file();

	/// The error probability the chain holds now.
	[[nodiscard]] field_math::TErrorProbability _error_probability() const;

	/// The link's log-likelihood, from the six counters and the current error probability.
	[[nodiscard]] double _link_log_likelihood() const;

	/// Writes the six counters of one tally, and adds them to the total the diagnostic reads.
	/// Opens the file on first use.
	///
	/// The trace is not behind a flag. Six integers an iteration is what the error probability's
	/// whole likelihood rests on, and the AND diagnostic reads nothing else.
	void _trace_link_counters(size_t iteration);

	/// Reports the two parameter-free constraints to the log file. A violation means the link is
	/// wrong, which is a finding rather than a defect, so this throws nothing and fails nothing.
	void _report_link_diagnostic() const;

	/// The chain start: the field at the LOTUS records, both tree fields at the field, and every
	/// internal node at what its children make most likely.
	void _start_the_chain(const TDataModel &data_model);

	/// A fixed field has to come from somewhere. Throws when the run fixed the field and gave it
	/// no states, which is a user error rather than a chain to run.
	void _throw_if_the_fixed_field_is_empty() const;

	/// Puts both tree fields at the field, and tallies the six counters over them. This is the
	/// second half of the chain start, and its only caller.
	void _hold_tree_fields_at_the_field();

	void _simulate_Y();

	/// The joint density of the configuration the iteration leaves behind: both trees' node
	/// states, the link, and the data (ADR-0005). The data term is what TDataModel sums, which is
	/// the LOTUS records and the simple error model.
	///
	/// A simulated chain draws from the prior and scores no data, so its data term is zero.
	[[nodiscard]] joint_density::TJointDensity
	_calculate_joint_density(const TDataModel &data_model);

	/// Writes one row of the joint density trace, opening the file on first use. Does nothing
	/// unless the run asked for the trace: the density is a pass over both node states, which is
	/// the cost `--write_joint_log_prob_density` exists to let a run skip.
	void _trace_joint_density(size_t iteration, const TDataModel &data_model);

	/// Counts every tree field cell that is a one now, on the iterations the field counts.
	void _count_the_tree_fields(size_t iteration);

	/// Writes the posterior of every tree field, one file per tree, beside the field's own.
	void _write_tree_field_posteriors() const;

	void _read_Y_from_file(const std::string &filename);

	template<bool FixZ> void _update_all_Z(size_t iteration) {
		if (iteration == 0 && ProgramOptions::WRITE_Z_TRACE && _Z_trace_files.empty() && !_fix_Z) {
			for (const auto &tree : _trees) {
				std::vector<size_t> Z_trace_header;
				Z_trace_header.reserve(tree->get_Z().total_size_of_container_space());
				for (size_t i = 0; i < tree->get_Z().total_size_of_container_space(); ++i) {
					Z_trace_header.push_back(i);
				}
				_Z_trace_files.emplace_back(_prefix + "_" + tree->get_tree_name() + "_Z_trace.txt",
				                            Z_trace_header, "\t");
			}
		}

		// The link each tree's leaves draw through. It carries the field, the *other* tree's node
		// state and the error probability -- none of which a tree may name (ADR-0005) -- and it is
		// built here, outside the parallel region the tree opens over its cliques.
		//
		// The trees are drawn one after another, so the second reads the leaf states the first was
		// just given.
		const field_math::TErrorProbability omega = _error_probability();
		for (size_t tree_idx = 0; tree_idx < _trees.size(); ++tree_idx) {
			// The other tree is the one this is not. The constructor has already said there are
			// exactly NUMBER_OF_TREES of them, and the link is written for two.
			const TNodeStateStorage &other_tree_field =
			    _trees[NUMBER_OF_TREES - 1 - tree_idx]->get_Z();
			const tree_field_link::TLeafLinks<TLinkPolicy, TFieldStorage, TNodeStateStorage> links(
			    _Y, other_tree_field, omega);
			_trees[tree_idx]->update_Z_and_nus_and_alphas_and_branch_lengths<FixZ>(iteration,
			                                                                      links);
		}
		if (_fix_Z) { return; }
		if (iteration % _Y.get_thinning_factor() == 0 && ProgramOptions::WRITE_Z_TRACE) {
			for (size_t tree_idx = 0; tree_idx < _trees.size(); ++tree_idx) {
				const auto &tree = _trees[tree_idx];
				_Z_trace_files[tree_idx].writeln(tree->get_Z().get_full_Z_binary_vector());
			}
		}
	}

	template<bool WriteFullY> void _write_Y_to_file(const std::string &filename) const {
		std::vector<std::string> header;
		header.emplace_back("position");
		header.emplace_back("Y_state");
		for (const auto &tree : _trees) { header.push_back(tree->get_tree_name()); }
		header.emplace_back("fraction_of_one");

		coretools::TOutputFile file(filename, header, "\t");
		if constexpr (WriteFullY) {
			this->_write_full_Y(file);
		} else {
			this->_write_only_values_in_Y_vector(file);
		}
	}

	void _write_full_Y(coretools::TOutputFile &file) const {
		std::array<size_t, 2> line{};
		for (size_t i = 0; i < _Y.total_size_of_container_space(); ++i) {
			// a missing cell reads as state 0, so a direct point lookup covers both cases
			line                 = {i, _Y.is_one(i)};
			auto leaf_index_of_Y = _Y.get_multi_dimensional_index(i);
			std::vector<std::string> node_names;
			for (size_t idx = 0; idx < leaf_index_of_Y.size(); ++idx) {
				size_t node_idx = _trees[idx]->get_node_index_from_leaf_index(leaf_index_of_Y[idx]);
				node_names.push_back(_trees[idx]->get_node_id(node_idx));
			};
			const double fraction = _Y.get_fraction_of_ones(i);
			file.writeln(line, node_names, fraction);
		}
	}

	void _write_only_values_in_Y_vector(coretools::TOutputFile &file) const {
		std::array<size_t, 2> line{};
		// iterate only the stored (non-default) cells, in ascending linear-index order
		for (const auto &[linear_index_in_y, storage] : _Y.get_stored_entries()) {
			// A cell that is not a one now and was never counted a one carries no posterior. Which
			// cells are *stored* is a property of the backend -- the sparse field holds the cells
			// it was given, ones and zeros alike, and the dense one holds the whole container
			// space -- so leaving those cells out is what makes this file say the same thing under
			// either backend, and it drops no information: every column of such a row is the
			// default.
			if (!storage.is_one() && storage.get_counter() == 0) { continue; }
			line                 = {linear_index_in_y, storage.is_one()};
			auto leaf_index_of_Y = _Y.get_multi_dimensional_index(linear_index_in_y);
			std::vector<std::string> node_names;
			for (size_t idx = 0; idx < leaf_index_of_Y.size(); ++idx) {
				size_t node_idx = _trees[idx]->get_node_index_from_leaf_index(leaf_index_of_Y[idx]);
				node_names.push_back(_trees[idx]->get_node_id(node_idx));
			};
			const double fraction = _Y.get_fraction_of_ones(linear_index_in_y);
			file.writeln(line, node_names, fraction);
		}
	}

public:
	TMarkovField(size_t n_iterations, std::vector<std::unique_ptr<TTree>> &Trees,
	             TypeParamErrorProbability *omega, std::string _prefix, bool simulate);
	~TMarkovField() = default;

	/// Puts the error probability's support at the open interval (0, 0.5).
	///
	/// The bound is a statement about the model, not a range on an argument (ADR-0005): at 0 the
	/// link is the deterministic AND and a draw takes log(0), and at 0.5 and above the
	/// tree fields are anti-correlated with the field. The type carries it, so the Metropolis
	/// proposal mirrors at both ends and never leaves the interval.
	///
	/// Must run before stattools sizes the parameter, because the value it creates then is checked
	/// against these bounds.
	static void set_error_probability_support();

	/// The log-likelihood ratio of the proposed error probability against the one it replaces.
	///
	/// O(1) in the number of cells: the counters do not move with the error probability, so this
	/// reads six integers and no cell. stattools has already proposed when this is called, so the
	/// parameter holds the proposal and remembers the old value.
	[[nodiscard]] double link_log_likelihood_ratio() const;

	// updates
	void update(TDataModel &data_model, size_t iteration);

	// simulation
	void simulate(TDataModel &data_model);

	// get Y
	[[nodiscard]] const TFieldStorage &get_Y_matrix() const;


	// functions to perform stuff on Y after burnin / MCMC finished
	void burninHasFinished();
	void MCMCHasFinished();
	void oneBurninHasFinished();
};

#endif // ACOL_TMARKOVFIELD_H
