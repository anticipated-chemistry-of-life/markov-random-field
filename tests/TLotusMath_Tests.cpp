// Properties of the LOTUS reporting model.
//
// The first four mirror the research-effort checks in model_validation/tests/test_independent.py,
// which until now could only be asked of the Python reference. The rest cover ground the reference
// cannot reach at all: it simulates with a single scalar gamma, while the C++ infers one per
// tree.

#include "lotus/TLotusMath.h"
#include "storages/TSparse.h"
#include "storages/y_storage/TStorageYDense.h"
#include "storages/y_storage/TStorageYSparse.h"
#include "gtest/gtest.h"

#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using lotus_math::TReportingModel;

IndexArray cell(size_t first, size_t second) { return IndexArray{first, second}; }

/// Two dimensions of three leaves each, with paper counts that differ within and across dimensions.
TReportingModel two_dimensions(double gamma_a, double gamma_b, double error_rate = 0.01) {
	return {{gamma_a, gamma_b}, error_rate, {{0, 1, 7}, {0, 3, 40}}};
}

} // namespace

// --------------------------------------------------------------------------
// Research effort
// --------------------------------------------------------------------------

TEST(LotusMath, research_effort_uses_log_paper_counts_not_raw_ones) {
	// The trap the validation README calls out: driving research effort from raw counts is
	// indistinguishable from an inference bug. log(count + 1) is the transform, and it now lives
	// here rather than in whatever reads the counts off disk.
	// A small gamma keeps both routes away from saturation, where they would agree at ~1 whatever
	// the transform: 40 papers give 0.17 through the log and 0.86 through the raw count.
	const double gamma = 0.05;
	const TReportingModel model({gamma}, 0.01, {{40}});

	const double from_log = 1.0 - std::exp(-gamma * std::log(41.0));
	const double from_raw = 1.0 - std::exp(-gamma * 40.0);

	EXPECT_DOUBLE_EQ(model.research_effort(IndexArray{0}), from_log);
	EXPECT_GT(std::abs(from_log - from_raw), 0.5) << "the two would be indistinguishable here";
}

TEST(LotusMath, a_leaf_with_no_papers_is_never_reported) {
	// occurrence_count(0) = log(1) = 0, so the factor is 1 - exp(0) = 0 exactly, and the product
	// across dimensions collapses. A pair nobody studied contributes no evidence of presence.
	const auto model = two_dimensions(1.5, 1.5);

	EXPECT_DOUBLE_EQ(model.research_effort(cell(0, 2)), 0.0);
	EXPECT_DOUBLE_EQ(model.research_effort(cell(2, 0)), 0.0);
	EXPECT_DOUBLE_EQ(model.research_effort(cell(0, 0)), 0.0);
	EXPECT_GT(model.research_effort(cell(2, 2)), 0.0);
}

TEST(LotusMath, research_effort_is_a_probability) {
	for (const double gamma : {0.01, 0.5, 1.1, 5.0, 50.0}) {
		const auto model = two_dimensions(gamma, gamma);
		for (size_t first = 0; first < 3; ++first) {
			for (size_t second = 0; second < 3; ++second) {
				const double effort = model.research_effort(cell(first, second));
				EXPECT_GE(effort, 0.0) << "gamma " << gamma;
				EXPECT_LE(effort, 1.0) << "gamma " << gamma;
			}
		}
	}
}

TEST(LotusMath, research_effort_rises_with_papers_and_with_gamma) {
	const auto model = two_dimensions(1.0, 1.0);
	// leaf order within a dimension is (0, 1, 7) and (0, 3, 40) papers
	EXPECT_LT(model.research_effort(cell(1, 1)), model.research_effort(cell(2, 1)));
	EXPECT_LT(model.research_effort(cell(1, 1)), model.research_effort(cell(1, 2)));

	const auto keener = two_dimensions(3.0, 3.0);
	EXPECT_LT(model.research_effort(cell(2, 2)), keener.research_effort(cell(2, 2)));
}

TEST(LotusMath, effort_saturates_as_gamma_grows) {
	const auto model = two_dimensions(500.0, 500.0);
	EXPECT_NEAR(model.research_effort(cell(2, 2)), 1.0, 1e-9);
}

// --------------------------------------------------------------------------
// One gamma per tree
// --------------------------------------------------------------------------

TEST(LotusMath, each_dimension_uses_its_own_gamma) {
	// The Python reference simulates with a single scalar gamma for both trees, so this asymmetry
	// is invisible to it. Here the two dimensions are given different rates and identical counts,
	// which must produce different factors.
	const TReportingModel model({0.2, 2.0}, 0.01, {{9}, {9}});

	const double slow = 1.0 - std::exp(-0.2 * std::log(10.0));
	const double fast = 1.0 - std::exp(-2.0 * std::log(10.0));
	ASSERT_NE(slow, fast);

	EXPECT_DOUBLE_EQ(model.research_effort(IndexArray{0, 0}), slow * fast);
}

TEST(LotusMath, effort_is_the_product_across_dimensions) {
	const auto model = two_dimensions(0.9, 1.4);
	const TReportingModel only_first({0.9}, 0.01, {{0, 1, 7}});
	const TReportingModel only_second({1.4}, 0.01, {{0, 3, 40}});

	EXPECT_DOUBLE_EQ(model.research_effort(cell(2, 1)),
	                 only_first.research_effort(IndexArray{2}) *
	                     only_second.research_effort(IndexArray{1}));
}

TEST(LotusMath, a_gamma_per_dimension_is_required) {
	EXPECT_THROW((TReportingModel({1.0}, 0.01, {{1}, {2}})), std::invalid_argument);
	EXPECT_THROW((TReportingModel({1.0, 2.0}, 0.01, {{1}})), std::invalid_argument);
}

// --------------------------------------------------------------------------
// The emission
// --------------------------------------------------------------------------

TEST(LotusMath, a_present_pair_is_reported_at_its_research_effort) {
	const auto model    = two_dimensions(1.1, 1.1);
	const auto position = cell(2, 2);
	const double effort = model.research_effort(position);

	EXPECT_DOUBLE_EQ(model.probability(true, true, position), effort);
	EXPECT_DOUBLE_EQ(model.probability(true, false, position), 1.0 - effort);
}

TEST(LotusMath, an_absent_pair_is_reported_at_the_flat_error_rate) {
	const double error_rate = 0.03;
	const auto model        = two_dimensions(1.1, 1.1, error_rate);

	EXPECT_DOUBLE_EQ(model.probability(false, true, cell(2, 2)), error_rate);
	EXPECT_DOUBLE_EQ(model.probability(false, false, cell(2, 2)), 1.0 - error_rate);
}

TEST(LotusMath, the_absent_case_does_not_depend_on_position) {
	// This is what lets the likelihood update answer absent cells without converting a linear
	// index, and what makes the bulk term for never-stored cells a single constant.
	const auto model = two_dimensions(1.1, 1.1, 0.03);
	for (size_t first = 0; first < 3; ++first) {
		for (size_t second = 0; second < 3; ++second) {
			EXPECT_DOUBLE_EQ(model.probability(false, true, cell(first, second)),
			                 model.probability_absent(true));
			EXPECT_DOUBLE_EQ(model.probability(false, false, cell(first, second)),
			                 model.probability_absent(false));
		}
	}
}

TEST(LotusMath, both_states_give_a_distribution_over_L) {
	const auto model = two_dimensions(1.1, 0.4, 0.03);
	for (size_t first = 0; first < 3; ++first) {
		for (size_t second = 0; second < 3; ++second) {
			for (const bool x : {false, true}) {
				const double reported = model.probability(x, true, cell(first, second));
				const double missing  = model.probability(x, false, cell(first, second));
				EXPECT_GE(reported, 0.0);
				EXPECT_GE(missing, 0.0);
				EXPECT_NEAR(reported + missing, 1.0, 1e-12)
				    << "x " << x << " at (" << first << ", " << second << ")";
			}
		}
	}
}

TEST(LotusMath, an_unstudied_present_pair_looks_absent_but_is_not_evidence) {
	// The asymmetry the source exists for: with no papers, a present pair is never reported, so a
	// missing record says nothing. A well-studied pair reports almost always, so a missing record
	// there is real evidence.
	const auto model = two_dimensions(2.0, 2.0);

	EXPECT_DOUBLE_EQ(model.probability(true, false, cell(0, 0)), 1.0);
	EXPECT_LT(model.probability(true, false, cell(2, 2)), 0.6);
}

// --------------------------------------------------------------------------
// The likelihood of the whole of L, read off the field's present cells
// --------------------------------------------------------------------------

namespace {

/// A small field and LOTUS on it, every cell drawn at random, with paper counts that repeat within
/// each tree -- so that several leaves share a class and the counts have to add them up. A leaf
/// with no papers is never reported when present, because the model gives that no probability.
struct TSmallField {
	std::vector<std::vector<size_t>> paper_counts{{0, 3, 3, 7, 0, 12, 3}, {5, 0, 5, 1, 40, 1}};
	std::vector<std::vector<bool>> present;
	std::vector<std::vector<bool>> reported;

	explicit TSmallField(uint64_t seed) {
		std::mt19937_64 rng(seed);
		std::bernoulli_distribution coin(0.4);
		for (size_t s = 0; s < paper_counts[0].size(); ++s) {
			present.emplace_back();
			reported.emplace_back();
			for (size_t m = 0; m < paper_counts[1].size(); ++m) {
				const bool x           = coin(rng);
				const bool unstudied   = paper_counts[0][s] == 0 || paper_counts[1][m] == 0;
				const bool reported_if = coin(rng);
				present.back().push_back(x);
				reported.back().push_back(reported_if && !(x && unstudied));
			}
		}
	}

	/// The counts, written out from what they are: the present cells by kind.
	[[nodiscard]] lotus_math::TPresentCells
	present_cells(const lotus_math::TPaperCountClasses &classes) const {
		size_t n_records = 0;
		for (const auto &row : reported) {
			for (const bool r : row) { n_records += static_cast<size_t>(r); }
		}
		lotus_math::TPresentCells cells(present.size() * present[0].size(), n_records, classes);
		for (size_t s = 0; s < present.size(); ++s) {
			for (size_t m = 0; m < present[s].size(); ++m) {
				if (!present[s][m]) { continue; }
				++cells.ones[static_cast<size_t>(reported[s][m])]
				            [cells.slot(classes.class_of(0, s), classes.class_of(1, m))];
			}
		}
		return cells;
	}

	/// The likelihood the counts stand in for: cell by cell, under the reporting model.
	[[nodiscard]] double cell_by_cell(const std::vector<double> &gammas, double error_rate) const {
		const TReportingModel model(gammas, error_rate, paper_counts);
		double sum = 0.0;
		for (size_t s = 0; s < present.size(); ++s) {
			for (size_t m = 0; m < present[s].size(); ++m) {
				sum += std::log(model.probability(present[s][m], reported[s][m], cell(s, m)));
			}
		}
		return sum;
	}
};

} // namespace

TEST(LotusMath, the_present_cells_give_the_likelihood_of_every_cell) {
	for (const uint64_t seed : {1U, 2U, 3U}) {
		const TSmallField field(seed);
		const lotus_math::TPaperCountClasses classes(field.paper_counts);
		const auto cells = field.present_cells(classes);
		for (const auto &[gammas, error_rate] : std::vector<std::pair<std::vector<double>, double>>{
		         {{0.3, 1.2}, 0.01}, {{2.0, 0.05}, 0.2}, {{0.9, 0.9}, 0.5}}) {
			const double expected = field.cell_by_cell(gammas, error_rate);
			ASSERT_TRUE(std::isfinite(expected));
			EXPECT_NEAR(lotus_math::log_likelihood_of_L(cells, classes, gammas, error_rate),
			            expected, 1e-9 * std::abs(expected))
			    << "seed " << seed << ", gamma " << gammas[0] << " " << gammas[1];
		}
	}
}

TEST(LotusMath, a_reported_present_cell_nobody_studied_has_no_probability) {
	// Research effort is 0 for a leaf with no papers, so the field that puts a reported cell there
	// is impossible, and says so rather than reading the empty kind as no evidence.
	const lotus_math::TPaperCountClasses classes({{0, 4}, {2}});
	lotus_math::TPresentCells cells(2, 1, classes);
	++cells.ones[1][cells.slot(classes.class_of(0, 0), classes.class_of(1, 0))];
	EXPECT_EQ(lotus_math::log_likelihood_of_L(cells, classes, {1.0, 1.0}, 0.1),
	          -std::numeric_limits<double>::infinity());
}

TEST(LotusMath, counts_of_more_cells_than_the_field_holds_are_refused) {
	const lotus_math::TPaperCountClasses classes({{1, 4}, {2}});
	lotus_math::TPresentCells cells(2, 1, classes);
	cells.ones[0][0] = 2;
	EXPECT_THROW((void)lotus_math::log_likelihood_of_L(cells, classes, {1.0, 1.0}, 0.1),
	             std::invalid_argument);
}

// --------------------------------------------------------------------------
// The field's present cells, counted from scratch
// --------------------------------------------------------------------------

namespace {

constexpr size_t N_ITERATIONS = 100;

/// The count, cell by cell. The tally reads every cell of the container space by point lookup,
/// present or not. It therefore shares no walk with the count it checks.
template<typename Field>
lotus_math::TPresentCells tally_every_cell(const Field &Y, const TSparseBinary &L,
                                           const lotus_math::TPaperCountClasses &classes) {
	const size_t n_cells   = Y.total_size_of_container_space();
	const size_t n_columns = Y.dimensions()[1];
	size_t n_records       = 0;
	for (size_t i = 0; i < n_cells; ++i) { n_records += static_cast<size_t>(L.is_one(i)); }

	lotus_math::TPresentCells cells(n_cells, n_records, classes);
	for (size_t i = 0; i < n_cells; ++i) {
		if (!Y.is_one(i)) { continue; }
		const uint32_t species  = classes.class_of(0, i / n_columns);
		const uint32_t molecule = classes.class_of(1, i % n_columns);
		++cells.ones[static_cast<size_t>(L.is_one(i))][cells.slot(species, molecule)];
	}
	return cells;
}

/// A field of random shape and fill, with paper counts drawn from a few values so that several
/// leaves share a class. The field also stores zeros, which are cells the count must not take
/// for present.
template<typename Field> struct TRandomField {
	IndexArray shape;
	std::vector<std::vector<size_t>> paper_counts;
	Field Y;

	explicit TRandomField(std::mt19937_64 &rng, double fill) {
		std::uniform_int_distribution<size_t> n_leaves(1, 40);
		shape = IndexArray{n_leaves(rng), n_leaves(rng)};

		std::uniform_int_distribution<size_t> papers(0, 6);
		paper_counts.resize(2);
		for (size_t dimension = 0; dimension < 2; ++dimension) {
			for (size_t leaf = 0; leaf < shape[dimension]; ++leaf) {
				paper_counts[dimension].push_back(papers(rng));
			}
		}

		Y = Field(N_ITERATIONS, shape);
		std::bernoulli_distribution present(fill);
		std::bernoulli_distribution stored_zero(0.2);
		for (size_t i = 0; i < Y.total_size_of_container_space(); ++i) {
			if (present(rng)) {
				Y.insert_one(i);
			} else if (stored_zero(rng)) {
				Y.insert_zero(i);
			}
		}
	}

	/// Records at random, each cell reported with probability `rate`.
	[[nodiscard]] TSparseBinary random_records(std::mt19937_64 &rng, double rate) const {
		TSparseBinary L(shape);
		std::bernoulli_distribution reported(rate);
		for (size_t i = 0; i < L.total_size_of_container_space(); ++i) {
			if (reported(rng)) { L.insert_one(i); }
		}
		return L;
	}

	/// Records on exactly the present cells, or on exactly the absent ones.
	[[nodiscard]] TSparseBinary records_where_present_is(bool state) const {
		TSparseBinary L(shape);
		for (size_t i = 0; i < L.total_size_of_container_space(); ++i) {
			if (Y.is_one(i) == state) { L.insert_one(i); }
		}
		return L;
	}
};

uint64_t total(const std::vector<uint64_t> &table) {
	uint64_t sum = 0;
	for (const uint64_t n : table) { sum += n; }
	return sum;
}

class FieldNames {
public:
	template<typename Field> static std::string GetName(int /*index*/) {
		if constexpr (std::is_same_v<Field, TStorageYSparse>) {
			return "sparse";
		} else {
			return "dense";
		}
	}
};

} // namespace

template<typename Field> class PresentCellCount : public ::testing::Test {};
using Fields = ::testing::Types<TStorageYSparse, TStorageYDense>;
TYPED_TEST_SUITE(PresentCellCount, Fields, FieldNames);

TYPED_TEST(PresentCellCount, matches_a_tally_of_every_cell_on_random_fields) {
	std::mt19937_64 rng(20261022);
	std::uniform_real_distribution<double> rate(0.05, 0.7);
	uint64_t reported   = 0;
	uint64_t unreported = 0;
	for (size_t trial = 0; trial < 50; ++trial) {
		const TRandomField<TypeParam> field(rng, rate(rng));
		const lotus_math::TPaperCountClasses classes(field.paper_counts);
		const TSparseBinary L = field.random_records(rng, rate(rng));

		const auto counted = lotus_math::count_present_cells(field.Y, L, classes);
		ASSERT_EQ(counted, tally_every_cell(field.Y, L, classes)) << "trial " << trial;
		reported += total(counted.ones[1]);
		unreported += total(counted.ones[0]);
	}
	// The trials count both halves. A swap of the two fails this check.
	EXPECT_GT(reported, 0U);
	EXPECT_GT(unreported, 0U);
}

TYPED_TEST(PresentCellCount, every_present_cell_reported_fills_the_reported_half_alone) {
	std::mt19937_64 rng(1);
	for (size_t trial = 0; trial < 10; ++trial) {
		const TRandomField<TypeParam> field(rng, 0.4);
		const lotus_math::TPaperCountClasses classes(field.paper_counts);
		const TSparseBinary L = field.records_where_present_is(true);

		const auto counted = lotus_math::count_present_cells(field.Y, L, classes);
		ASSERT_EQ(counted, tally_every_cell(field.Y, L, classes)) << "trial " << trial;
		EXPECT_EQ(total(counted.ones[0]), 0U) << "trial " << trial;
		EXPECT_EQ(total(counted.ones[1]), L.number_of_ones()) << "trial " << trial;
	}
}

TYPED_TEST(PresentCellCount, no_present_cell_reported_fills_the_unreported_half_alone) {
	std::mt19937_64 rng(2);
	for (size_t trial = 0; trial < 10; ++trial) {
		const TRandomField<TypeParam> field(rng, 0.4);
		const lotus_math::TPaperCountClasses classes(field.paper_counts);
		const TSparseBinary L = field.records_where_present_is(false);

		const auto counted = lotus_math::count_present_cells(field.Y, L, classes);
		ASSERT_EQ(counted, tally_every_cell(field.Y, L, classes)) << "trial " << trial;
		EXPECT_EQ(total(counted.ones[1]), 0U) << "trial " << trial;
		EXPECT_EQ(total(counted.ones[0]),
		          field.Y.total_size_of_container_space() - L.number_of_ones())
		    << "trial " << trial;
	}
}

TYPED_TEST(PresentCellCount, a_field_with_no_present_cell_counts_none) {
	std::mt19937_64 rng(3);
	for (size_t trial = 0; trial < 10; ++trial) {
		// Stored zeros only: the field has entries, and none of them is present.
		const TRandomField<TypeParam> field(rng, 0.0);
		const lotus_math::TPaperCountClasses classes(field.paper_counts);
		const TSparseBinary L = field.random_records(rng, 0.3);

		const auto counted = lotus_math::count_present_cells(field.Y, L, classes);
		const lotus_math::TPresentCells none(field.Y.total_size_of_container_space(),
		                                     L.number_of_ones(), classes);
		EXPECT_EQ(counted, none) << "trial " << trial;
		EXPECT_EQ(counted, tally_every_cell(field.Y, L, classes)) << "trial " << trial;
	}
}
