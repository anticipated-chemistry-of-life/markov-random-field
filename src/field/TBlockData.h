//
// What the data sources say to one block update: one pass per data source, folded into one cell
// data factor per cell (ADR-0017).
//
// The block model (field/TBlockModel.h) asks this about the data sources, and asks the trees about
// the rest. This holds a pass for each data source (data_sources/data_source_pass.h), and begins
// and ends them. It fans each row's hooks out to every pass, and multiplies their cell factors into
// the cell data factor (GLOSSARY.md). It is a compile-time fold, because the cell factor runs once
// per leaf pair (ADR-0017).
//
// The fold multiplies in list order, into a running product that starts from 1.0. ADR-0017 says
// why that keeps the chain byte-identical. `TCompiledBlockData`, at the end of this file, is the
// only place the block update reads the compiled-in list (data_sources/block_sources.h).
//

#pragma once

#include "coretools/Main/TError.h"
#include "data_sources/block_sources.h"
#include "data_sources/data_source_pass.h"
#include "storages/storage_concepts.h"
#include <array>
#include <cstddef>
#include <tuple>
#include <type_traits>
#include <utility>

template<DataSourcePass... Passes> class TBlockData {
private:
	std::tuple<Passes...> _passes;

	/// Calls `hook(pass, row)` on every pass and its own row state, in list order.
	template<typename Self, typename Rows, typename Hook>
	static void _each(Self &passes, Rows &rows, Hook &&hook) {
		[&]<size_t... I>(std::index_sequence<I...>) {
			(hook(std::get<I>(passes), std::get<I>(rows)), ...);
		}(std::index_sequence_for<Passes...>{});
	}

public:
	/// One row's state for every pass, and, in a debug build, where the row is in its walk.
	struct TRow {
		std::tuple<PassRow<Passes>...> passes;
#ifndef NDEBUG
		/// The lowest molecule leaf the row may still be asked about.
		size_t next_molecule_leaf = 0;
		/// Whether a cell has had its factor asked for and not yet been recorded.
		bool awaiting_record      = false;
		/// The molecule leaf of the cell whose factor was last asked for.
		size_t asked              = 0;
#endif
	};

	explicit TBlockData(Passes... passes) : _passes(std::move(passes)...) {}

	/// The fold over the passes `begin(std::type_identity<Pass>{})` hands out, one per pass type of
	/// the list, in list order. How a data source hands its pass out is the caller's business.
	template<typename Begin> [[nodiscard]] static TBlockData begin_with(Begin &&begin) {
		return TBlockData(begin(std::type_identity<Passes>{})...);
	}

	/// Once, serial, before the parallel region.
	void begin_pass() {
		std::apply([](Passes &...pass) { (pass.begin_pass(), ...); }, _passes);
	}

	/// Once per species leaf, inside the parallel region.
	[[nodiscard]] TRow begin_row(size_t species_leaf, size_t first_field_cell) const {
		return std::apply(
		    [&](const Passes &...pass) {
			    return TRow{.passes = {pass.begin_row(species_leaf, first_field_cell)...}};
		    },
		    _passes);
	}

	/// Multiplies every pass's answer for one cell into `product`, in list order, and returns it.
	/// `product` is 1.0 unless a data source not yet migrated to a pass was multiplied in ahead of
	/// the list (field/TBlockModel.h).
	///
	/// A debug build asserts the ordering guarantee every pass relies on: within a row, ascending
	/// molecule leaves, each once, the cell factor before the record.
	[[nodiscard]] std::array<double, 2>
	cell_factor(TRow &row, size_t molecule_leaf,
	            std::array<double, 2> product = {1.0, 1.0}) const {
#ifndef NDEBUG
		DEBUG_ASSERT(!row.awaiting_record);
		DEBUG_ASSERT(molecule_leaf >= row.next_molecule_leaf);
		row.awaiting_record = true;
		row.asked           = molecule_leaf;
#endif
		_each(_passes, row.passes, [&](const auto &pass, auto &pass_row) {
			const std::array<double, 2> factor = pass.cell_factor(pass_row, molecule_leaf);
			product[0] *= factor[0];
			product[1] *= factor[1];
		});
		return product;
	}

	/// Tells every pass the field state the cell was given.
	void record(TRow &row, size_t molecule_leaf, bool drawn_y) const {
#ifndef NDEBUG
		DEBUG_ASSERT(row.awaiting_record && row.asked == molecule_leaf);
		row.awaiting_record    = false;
		row.next_molecule_leaf = molecule_leaf + 1;
#endif
		_each(_passes, row.passes, [&](const auto &pass, auto &pass_row) {
			pass.record(pass_row, molecule_leaf, drawn_y);
		});
	}

	/// Once per row, inside the parallel region: every pass folds the row's tally into its own.
	void end_row(TRow &row) {
#ifndef NDEBUG
		DEBUG_ASSERT(!row.awaiting_record);
#endif
		_each(_passes, row.passes, [](auto &pass, auto &pass_row) { pass.end_row(pass_row); });
	}

	/// Once, serial, after the parallel region, with the field the update left.
	template<FieldStorage Field> void end_pass(const Field &Y) {
		std::apply([&](Passes &...pass) { (pass.end_pass(Y), ...); }, _passes);
	}

	/// One pass, to read what it tallied.
	template<typename Pass> [[nodiscard]] const Pass &pass() const {
		return std::get<Pass>(_passes);
	}

	/// Every pass, to hand each back to its data source once the update has ended.
	[[nodiscard]] std::tuple<Passes...> take_passes() && { return std::move(_passes); }
};

template<typename List> struct TBlockDataOver;
template<DataSourcePass... Passes> struct TBlockDataOver<TPassList<Passes...>> {
	using type = TBlockData<Passes...>;
};

/// The fold over every pass this build compiles in.
using TCompiledBlockData = typename TBlockDataOver<TBlockSources>::type;
