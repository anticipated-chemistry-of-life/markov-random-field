//
// The likelihood of one LOTUS cell.
//
// A LOTUS record is an occurrence someone reported in the literature, so its probability is not a
// flat error rate: a present pair is reported only if anyone looked, and how hard anyone looked is
// the research effort spent on that pair (CONTEXT.md, "Research effort", "LOTUS record").
//
// Research effort is a product over the trees of `1 - exp(-gamma_i * log(count_i + 1))`,
// where `count_i` is the paper count of the leaf that cell occupies in dimension `i` and `gamma_i`
// is that dimension's detection rate. Note there is one gamma *per tree*, not one
// overall. Given the latent state `x` of the cell, the emission is then
//
//     P(L = 1 | x = 1) = research effort          P(L = 1 | x = 0) = error rate
//     P(L = 0 | x = 1) = 1 - research effort      P(L = 0 | x = 0) = 1 - error rate
//
// so an absent pair is reported at a flat rate while a present one is reported at a rate that
// depends on where it sits. That asymmetry is the whole point of the source: a missing record for a
// well-studied pair is evidence of absence, and for an unstudied pair it is evidence of nothing.
//
// The transform `log(count + 1)` lives here rather than in the reader, because CONTEXT.md makes
// the raw publication count the input to research effort. Applying it while reading the file gave
// a reader that returned something other than paper counts, which the validation harness records
// as a trap: simulating LOTUS data from raw counts is indistinguishable from an inference bug.
//
// A grid is immutable. It is a value derived from (gammas, error rate, paper counts), so a
// Metropolis proposal on gamma builds a second one and the caller keeps whichever it accepts; there
// is no refresh hook and nothing to revert.
//
// The likelihood of the whole of L does not need the field cell by cell. A cell's term depends on
// its field state, on whether LOTUS reports it, and -- for a present cell -- on the paper counts of
// its two leaves, nothing else. So the field is summarised by how many present cells there are of
// each (species count, molecule count, reported) kind, and the likelihood is a sum over the kinds:
// a few hundred thousand terms in place of a pass over every cell. The field moves only in the
// block update, which keeps that summary as it draws; gamma and the error rate move against it.
// See `TPresentCells` and `log_likelihood_of_L` below, and ADR-0015.
//
// Everything here is free-standing: it is deliberately NOT guarded by USE_LOTUS, so it compiles
// (and is unit-tested) in every build configuration.
//

#pragma once

#include "constants.h"
#include "storages/storage_concepts.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace lotus_math {

/// The paper count as research effort consumes it: `log(count + 1)`.
///
/// A leaf with no papers gives 0, which drives its factor to `1 - exp(0) = 0`, so a pair nobody
/// has studied is never reported. That is the intended behaviour, not a degenerate case.
inline double occurrence_count(size_t paper_count) {
	return std::log(static_cast<double>(paper_count) + 1.0);
}

/// The probability that one tree's side of a present pair is studied: `1 - exp(-gamma *
/// log(count + 1))`. Research effort is the product of this over the trees.
inline double detection_factor(double gamma, size_t paper_count) {
	return 1.0 - std::exp(-gamma * occurrence_count(paper_count));
}

/// The leaves of every tree, grouped by paper count.
///
/// Research effort reads a leaf only through its paper count, so two leaves with one count are the
/// same leaf as far as LOTUS is concerned. A tree has tens of thousands of leaves and a few hundred
/// distinct counts, which is what makes `TPresentCells` small. A class is the rank of its count
/// among the tree's distinct counts.
class TPaperCountClasses {
private:
	std::vector<std::vector<uint32_t>> _class_of_leaf; ///< [dimension][leaf]
	std::vector<std::vector<size_t>> _count_of_class;  ///< [dimension][class], ascending

public:
	TPaperCountClasses() = default;
	explicit TPaperCountClasses(const std::vector<std::vector<size_t>> &paper_counts) {
		_class_of_leaf.resize(paper_counts.size());
		_count_of_class.resize(paper_counts.size());
		for (size_t dimension = 0; dimension < paper_counts.size(); ++dimension) {
			const auto &counts = paper_counts[dimension];
			auto &distinct     = _count_of_class[dimension];
			distinct.assign(counts.begin(), counts.end());
			std::sort(distinct.begin(), distinct.end());
			distinct.erase(std::unique(distinct.begin(), distinct.end()), distinct.end());

			auto &class_of_leaf = _class_of_leaf[dimension];
			class_of_leaf.resize(counts.size());
			for (size_t leaf = 0; leaf < counts.size(); ++leaf) {
				const auto at = std::lower_bound(distinct.begin(), distinct.end(), counts[leaf]);
				class_of_leaf[leaf] = static_cast<uint32_t>(at - distinct.begin());
			}
		}
	}

	[[nodiscard]] size_t n_dimensions() const { return _count_of_class.size(); }
	[[nodiscard]] size_t n_classes(size_t dimension) const {
		return _count_of_class[dimension].size();
	}
	[[nodiscard]] uint32_t class_of(size_t dimension, size_t leaf) const {
		return _class_of_leaf[dimension][leaf];
	}
	[[nodiscard]] size_t count_of(size_t dimension, size_t paper_class) const {
		return _count_of_class[dimension][paper_class];
	}
};

/// What the likelihood of L needs to know about the field: how many cells are present in it, by
/// the paper-count classes of their two leaves and by whether LOTUS reports them.
///
/// An absent cell's term does not depend on where it is, so absent cells are not counted: how many
/// there are, reported or not, follows from the totals. `ones[L]` is indexed by
/// `species_class * n_molecule_classes + molecule_class`.
///
/// Every entry is an integer, so counts added up in any order are the same counts: the block update
/// adds its rows in from many threads at once, and the likelihood read from them does not depend
/// on the schedule, the thread count or the backend.
struct TPresentCells {
	/// Every cell of the field.
	size_t n_cells            = 0;
	/// The cells LOTUS reports, present or not.
	size_t n_records          = 0;
	size_t n_species_classes  = 0;
	size_t n_molecule_classes = 0;
	/// The present cells, `ones[reported][species_class * n_molecule_classes + molecule_class]`.
	std::array<std::vector<uint64_t>, 2> ones;

	TPresentCells() = default;
	TPresentCells(size_t n_cells, size_t n_records, const TPaperCountClasses &classes)
	    : n_cells(n_cells), n_records(n_records), n_species_classes(classes.n_classes(0)),
	      n_molecule_classes(classes.n_classes(1)) {
		static_assert(NUMBER_OF_TREES == 2, "a kind is a species class by a molecule class");
		for (auto &table : ones) { table.assign(n_species_classes * n_molecule_classes, 0); }
	}

	[[nodiscard]] size_t slot(uint32_t species_class, uint32_t molecule_class) const {
		return static_cast<size_t>(species_class) * n_molecule_classes + molecule_class;
	}
	bool operator==(const TPresentCells &) const = default;
};

/// Counts a field's present cells from scratch. `L` holds the LOTUS records and has the field's
/// dimensions. A cell of the field is therefore the same cell of `L`.
///
/// The count walks the field's ones and nothing else. A forward cursor over the records says
/// whether LOTUS reports each one. The records come in the same ascending order. The cost is
/// O(ones(Y) + ones(L)), not O(total cells). The block update counts the present cells as it
/// draws. This count is for a field that arrived some other way, and for checking that count.
template<FieldStorage Field, BinaryStorage Records>
[[nodiscard]] TPresentCells count_present_cells(const Field &Y, const Records &L,
                                                const TPaperCountClasses &classes) {
	static_assert(NUMBER_OF_TREES == 2, "a cell is a species leaf by a molecule leaf");
	TPresentCells cells(L.total_size_of_container_space(), L.number_of_ones(), classes);
	const size_t n_columns = Y.dimensions()[1];

	auto record = L.ones_cursor();
	for (auto present = Y.ones_cursor(); present.valid(); present.advance()) {
		const size_t cell = present.linear_index();
		while (record.valid() && record.linear_index() < cell) { record.advance(); }
		const bool reported = record.valid() && record.linear_index() == cell;
		const size_t slot   = cells.slot(classes.class_of(0, cell / n_columns),
		                                 classes.class_of(1, cell % n_columns));
		++cells.ones[static_cast<size_t>(reported)][slot];
	}
	return cells;
}

/// `n * log(p)`, read as 0 when there is nothing to count: a kind of cell that the field does not
/// hold contributes nothing, even where the model gives it no probability at all.
inline double weighted_log(uint64_t n, double p) {
	return n == 0 ? 0.0 : static_cast<double>(n) * std::log(p);
}

/// The log-likelihood of the whole of L, given the field's present cells:
///
///     sum over present cells of log P(L | x = 1, effort of the cell)
///   + (absent reported cells) log(error rate) + (absent unreported cells) log(1 - error rate)
///
/// The present cells are summed by kind, in class order, so the answer is a function of the counts
/// and the parameters alone.
[[nodiscard]] inline double log_likelihood_of_L(const TPresentCells &cells,
                                                const TPaperCountClasses &classes,
                                                const std::vector<double> &gammas,
                                                double error_rate) {
	if (gammas.size() != 2 || classes.n_dimensions() != 2) {
		throw std::invalid_argument("The LOTUS likelihood needs a gamma and classes per tree.");
	}
	std::array<std::vector<double>, 2> factor;
	for (size_t dimension = 0; dimension < 2; ++dimension) {
		factor[dimension].resize(classes.n_classes(dimension));
		for (size_t c = 0; c < classes.n_classes(dimension); ++c) {
			factor[dimension][c] =
			    detection_factor(gammas[dimension], classes.count_of(dimension, c));
		}
	}

	uint64_t present_reported   = 0;
	uint64_t present_unreported = 0;
	double present              = 0.0;
	for (uint32_t s = 0; s < cells.n_species_classes; ++s) {
		// One row of kinds at a time, so that its terms are added among themselves before they
		// meet the running total, which is the larger by far.
		double row = 0.0;
		for (uint32_t m = 0; m < cells.n_molecule_classes; ++m) {
			const size_t slot         = cells.slot(s, m);
			const uint64_t reported   = cells.ones[1][slot];
			const uint64_t unreported = cells.ones[0][slot];
			const double effort       = factor[0][s] * factor[1][m];
			row += weighted_log(reported, effort) + weighted_log(unreported, 1.0 - effort);
			present_reported += reported;
			present_unreported += unreported;
		}
		present += row;
	}
	// More present cells than the field holds, of either kind, are the counts of another field.
	if (cells.n_records > cells.n_cells || present_reported > cells.n_records ||
	    present_unreported > cells.n_cells - cells.n_records) {
		throw std::invalid_argument(
		    "The LOTUS present-cell counts hold more cells than the field they summarise.");
	}
	const uint64_t absent_reported   = cells.n_records - present_reported;
	const uint64_t absent_unreported = cells.n_cells - cells.n_records - present_unreported;
	return present + weighted_log(absent_reported, error_rate) +
	       weighted_log(absent_unreported, 1.0 - error_rate);
}

/// One clique of the observation model: the probability of a LOTUS record, per cell, given the
/// latent state of that cell.
class TReportingModel {
private:
	/// `_factor[dimension][leaf]` caches `1 - exp(-gamma_i * occurrence_count(count))`, so the
	/// per-cell path in the field update is a product of lookups with no exp() and no parameter
	/// access.
	std::vector<std::vector<double>> _factor;
	double _error_rate = 0.0;

public:
	/// @param gammas       One detection rate per tree.
	/// @param error_rate   The rate at which an absent pair is reported anyway.
	/// @param paper_counts Raw publication counts, `paper_counts[dimension][leaf]`.
	TReportingModel(const std::vector<double> &gammas, double error_rate,
	                const std::vector<std::vector<size_t>> &paper_counts)
	    : _error_rate(error_rate) {
		if (gammas.size() != paper_counts.size()) {
			throw std::invalid_argument(
			    "TReportingModel needs one gamma per dimension of paper counts.");
		}

		_factor.resize(gammas.size());
		for (size_t dimension = 0; dimension < gammas.size(); ++dimension) {
			const auto &counts = paper_counts[dimension];
			auto &factor       = _factor[dimension];
			factor.resize(counts.size());
			for (size_t leaf = 0; leaf < counts.size(); ++leaf) {
				factor[leaf] = detection_factor(gammas[dimension], counts[leaf]);
			}
		}
	}

	[[nodiscard]] size_t n_dimensions() const { return _factor.size(); }
	[[nodiscard]] double error_rate() const { return _error_rate; }

	/// The probability that a *present* pair at `cell` is reported. `cell` is indexed in leaf
	/// space: one leaf index per tree.
	[[nodiscard]] double research_effort(const IndexArray &cell) const {
		double product = 1.0;
		for (size_t dimension = 0; dimension < _factor.size(); ++dimension) {
			product *= _factor[dimension][cell[dimension]];
		}
		return product;
	}

	/// `P(L | x)` for the cell at `cell`.
	[[nodiscard]] double probability(bool x, bool L, const IndexArray &cell) const {
		if (!x) { return probability_absent(L); }
		const double effort = research_effort(cell);
		return L ? effort : 1.0 - effort;
	}

	/// `P(L | x = 0)`, which does not depend on where the cell is.
	///
	/// Separate from `probability` so a caller walking a sparse container can answer the absent
	/// case without building a cell index for it. Most cells are absent, so that conversion is the
	/// bulk of the work it would otherwise do.
	[[nodiscard]] double probability_absent(bool L) const {
		return L ? _error_rate : 1.0 - _error_rate;
	}
};

} // namespace lotus_math
