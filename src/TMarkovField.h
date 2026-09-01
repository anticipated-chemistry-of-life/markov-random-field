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
#include "mass_spec/msms_data.h"
#include "random/TCellUniforms.h"
#include "storages/storage_backend.h"
#include "tree/TTree.h"
#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

//-----------------------------------
// Field update bookkeeping
//-----------------------------------

class TDataModel; // forward declaration

/// Per-cell outcome of one Y update. Each data source keeps its own likelihood bookkeeping (they
/// are independent terms), so the results are handed back separately instead of merged.
struct TYUpdateResult {
#ifdef USE_LOTUS
	/// P(L_cell | x = new_state). Neutral value 1.0 (log 0).
	double prob_lotus_new_state = 1.0;
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
	/// Whether the observed D cell contradicts the state Y was just set to.
	bool simple_model_disagrees = false;
#endif
};

/// Per-thread accumulators for one full field update, committed to the data sources at the end.
///
/// The accumulators are bundled into a single object on purpose: `#ifdef` cannot appear inside a
/// `#pragma omp` line, and the update's `default(none) shared(...)` clause has to name every
/// variable it touches. One object keeps that clause identical in every build configuration.
class TDataUpdateAccumulator {
private:
#ifdef USE_LOTUS
	std::vector<coretools::TSumLogProbability> _lotus_LL;
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
	std::vector<size_t> _n_disagree;
#endif

public:
	/// Sizing happens in the body rather than in a member-initializer list, so that adding or
	/// removing a source does not require rebalancing the commas of a #ifdef'd init list.
	explicit TDataUpdateAccumulator(size_t n_threads) {
#ifdef USE_LOTUS
		_lotus_LL.resize(n_threads);
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
		_n_disagree.assign(n_threads, 0);
#endif
	}

	/// Hot path: called once per updated Y cell, from inside the parallel region. Only ever touches
	/// the slot of the calling thread.
	void add(size_t thread, const TYUpdateResult &result) {
#ifdef USE_LOTUS
		_lotus_LL[thread].add(result.prob_lotus_new_state);
#endif
#ifdef USE_SIMPLE_ERROR_MODEL
		_n_disagree[thread] += static_cast<size_t>(result.simple_model_disagrees);
#endif
	}

	/// Sums the per-thread slots and installs the results in the data sources. Called once, after
	/// the parallel region.
	void commit(TDataModel &data_model);
};

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

	// mass spectrometry data and the dimension indices of molecules/species within _trees
	std::optional<TMSMSData> _ms_data;

	// complete joint density of the markov random field
	std::vector<double> _complete_log_density;

	/// Was Z initialized from children ?
	bool _z_initialized_from_children = false;

	// output files
	coretools::TOutputFile _Y_trace_file;
	std::vector<coretools::TOutputFile> _Z_trace_files;
	coretools::TOutputFile _joint_density_file;

	// The field update: one species leaf per thread, read and written through windows.
	//
	// Both templates are defined in TMarkovField.cpp. Every caller of them is in that file, and
	// the update reads the data sources, which this header only forward-declares.

	/// Draws every cell of one species leaf's row of the field.
	///
	/// The row is a thread's whole share of the update. A field cell's Markov blanket holds no
	/// other field cell, so the rows are conditionally independent given the node states. The call
	/// hands out the cells its window could not write in place, through `deferred_inserts`.
	template<bool IsSimulation, bool initYFromData>
	void _update_field_row(size_t species_leaf, TDataModel &data_model,
	                       const TCellUniforms &uniforms, TDataUpdateAccumulator &accumulator,
	                       std::vector<size_t> &deferred_inserts);

	/// One field update: every species leaf, one thread each.
	template<bool IsSimulation, bool initYFromData>
	void _update_all_Y(TDataModel &data_model, size_t iteration);

	void _simulate_Y();
	double _calculate_complete_joint_density();
	void _reset_log_joint_density() {
		_complete_log_density.clear();
		_complete_log_density.resize(ProgramOptions::NUMBER_OF_THREADS);
	}

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
			_tree->update_Z_and_nus_and_alphas_and_branch_lengths<IsSimulation, FixZ>(_Y,
			                                                                          iteration);
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
