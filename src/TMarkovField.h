//
// Created by madleina on 10.01.25.
//

#ifndef ACOL_TMARKOVFIELD_H
#define ACOL_TMARKOVFIELD_H

#include "TClique.h"

#include "Types.h"
#include "cli.h"
#include "constants.h"
#include "coretools/Files/TOutputFile.h"
#include "coretools/Main/TError.h"
#include "field/TFieldMath.h"
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
private:
	// trees and Y
	std::vector<std::unique_ptr<TTree>> &_trees;
	TFieldStorage _Y;
	std::string _prefix;

	// fix values?
	bool _fix_Y = false;
	bool _fix_Z = false;

	// Mass spectrometry data, still dormant. Nothing builds it. The block update does not read it
	// either: the eight-state block takes the LOTUS and the simple-error term, and adapting a third
	// source is that source's own work.
	std::optional<TMSMSData> _ms_data;

	// The error probability standing between the two tree fields and the field. Held fixed for
	// now; the Metropolis move that estimates it is #37.
	field_math::TErrorProbability _omega{ProgramOptions::ERROR_PROBABILITY};

	// The link's sufficient statistic over the whole field, n(bucket, field state). The block
	// update tallies it as it goes and commits it here, which is what makes the error
	// probability's likelihood O(1) in the number of cells (ADR-0005).
	field_math::TLinkCounters _link_counters;

	// The two tree factors of every leaf pair, at the states the last block update drew. The
	// link's own term is not in here: it comes from the counters above, for the whole field at
	// once.
	double _block_log_density = 0.0;

	/// Was Z initialized from children ?
	bool _z_initialized_from_children = false;

	// output files
	coretools::TOutputFile _Y_trace_file;
	std::vector<coretools::TOutputFile> _Z_trace_files;
	coretools::TOutputFile _joint_density_file;

	/// One block update: the field and both tree fields at every leaf pair, one species leaf per
	/// thread. Defined in TMarkovField.cpp, where the model it hands the traversal is complete.
	template<bool IsSimulation, bool InitYFromData>
	void _update_block(TDataModel &data_model, size_t iteration);

	/// Opens the field's trace file on the first iteration of a chain.
	void _open_Y_trace_file(bool is_simulation);

	/// Puts both tree fields at the field, and tallies the six counters over them. Only the fixed
	/// field path needs it: a block update writes all three and leaves the tally behind as it goes.
	void _hold_tree_fields_at_the_field();

	void _simulate_Y();
	double _calculate_complete_joint_density();

	void _read_Y_from_file(const std::string &filename);

	template<bool IsSimulation, bool FixZ> void _update_all_Z(size_t iteration) {
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

		for (auto &_tree : _trees) {
			_tree->update_Z_and_nus_and_alphas_and_branch_lengths<IsSimulation, FixZ>(iteration);
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
	             std::string _prefix);
	~TMarkovField() = default;

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

	static size_t get_num_iterations_simulation() { return ProgramOptions::NUM_ITERATIONS; }
};

#endif // ACOL_TMARKOVFIELD_H
