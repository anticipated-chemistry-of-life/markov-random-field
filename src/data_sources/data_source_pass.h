//
// A data source's part in one block update: its pass (ADR-0017).
//
// The block update draws every cell of the field once, and at each cell it multiplies in what
// every compiled-in data source says about that cell (the cell data factor, GLOSSARY.md). A data
// source takes part through a pass: a value it hands out when a block update begins
// (`begin_block_update()`) and takes back when it ends (`end_block_update(pass)`), to install the
// statistic the pass tallied. The pass is built from plain data -- what the source observes the
// field with, and its current parameter values -- so a test can build one without a tree or
// stattools.
//
// This is not `DataSource` (data_sources/data_source.h). That concept is a source's setup,
// likelihood, simulation and notifier stats, and neither it nor ACOL_DATA_SOURCES knows about
// passes. Mass spectrometry takes part in the block update without satisfying `DataSource`.
//
// The fold over every compiled-in pass is field/TBlockData.h, and the list of them is
// data_sources/block_sources.h.
//

#pragma once

#include "storages/storage_backend.h"
#include <array>
#include <concepts>
#include <cstddef>
#include <type_traits>
#include <utility>

/// The type a pass's `begin_row` hands back: the row's own state, which only that row's thread
/// touches.
template<typename P>
using PassRow =
    std::remove_cvref_t<decltype(std::declval<const P &>().begin_row(size_t{}, size_t{}))>;

/// A data source's part in one block update, through six hooks, in the order the block update
/// calls them:
///
/// - `begin_pass()`: once, serial, before the parallel region. Readies whatever cannot be readied
///   concurrently.
/// - `begin_row(species_leaf, first_field_cell)`: once per species leaf, inside the parallel
///   region. `first_field_cell` is the field's linear index of the row's cell at molecule leaf 0.
///   Returns the row's own state.
/// - `cell_factor(row, molecule_leaf)`: once per cell. Returns `{P(obs | Y = 0), P(obs | Y = 1)}`
///   for the cell, and may advance row-local state such as a cursor.
/// - `record(row, molecule_leaf, drawn_y)`: once per cell, after its `cell_factor`. Gets the field
///   state the cell was given, and nothing else: not the factors, not the tree field states.
/// - `end_row(row)`: once per row, inside the parallel region. Folds the row's tally into the
///   pass-wide one. Rows end concurrently, so this hook writes the pass from many threads at once.
///   The pass must do that safely, and so that the tally is the same at any thread count.
/// - `end_pass(Y)`: once, serial, after the parallel region. Finishes the pass-wide tally. Gets
///   the field the update left, read-only, so that a debug build can recount the statistic.
///
/// Within a row, cells arrive in ascending molecule-leaf order, each once, with the cell factor
/// before the record. That holds across tiles because a row's state survives from one tile to the
/// next, and a debug build of the fold asserts it (field/TBlockData.h). There is no per-column
/// hook.
///
/// `begin_row`, `cell_factor` and `record` are const: many rows run at once, and the only state
/// they may write is the row's own.
template<typename P>
concept DataSourcePass =
    std::movable<P> &&
    requires(P &pass, const P &read_only, PassRow<P> &row, size_t species_leaf,
	         size_t first_field_cell, size_t molecule_leaf, bool drawn_y, const TFieldStorage &Y) {
	    { pass.begin_pass() } -> std::same_as<void>;
	    { read_only.begin_row(species_leaf, first_field_cell) } -> std::same_as<PassRow<P>>;
	    { read_only.cell_factor(row, molecule_leaf) } -> std::same_as<std::array<double, 2>>;
	    { read_only.record(row, molecule_leaf, drawn_y) } -> std::same_as<void>;
	    { pass.end_row(row) } -> std::same_as<void>;
	    { pass.end_pass(Y) } -> std::same_as<void>;
    };
