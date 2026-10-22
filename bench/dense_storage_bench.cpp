//
// Memory and pass time of the dense storages against armadillo matrices holding the same cells.
//
// Four backings of the field and the two node states run the same chain:
//
//   own          TStorageYDense / TStorageZDense, driven through the production code paths:
//                block_update::run, clique_traversal::run over TCliqueView, node_state_walk and
//                node_state_density. Row-major: a field row (one species leaf) is contiguous, a
//                species clique (one molecule leaf) is a column a whole row width apart per node.
//   own_clique   the same storages and code paths, with each node state laid out as TTree lays
//                out its own: the tree's own dimension contiguous, so every clique is one run of
//                adjacent cells (ADR-0013).
//   arma         arma::Mat<u16> / arma::Mat<u8> addressed at the same (row, column) as `own`.
//                Armadillo is column-major, so this flips which runs are contiguous: a species
//                clique becomes contiguous, a molecule clique and a field row become strided.
//   arma_clique  arma::Mat as well, but every node state stored with its own tree's cliques as
//                contiguous columns (Z_s as nodes x molecule leaves, Z_m as nodes x species
//                leaves), and the field stored with a species leaf's row contiguous. Both walks are
//                contiguous; the block update reads the species tree field across a stride.
//
// Every cell draws the uniform hashed from its row-major (row, column) position (ADR-0007), so all
// four run one chain: the final states, the link counters and the clique densities are compared
// and the run fails if they differ. What differs is only where the cells sit in memory, which is
// what the timings measure.
//
// The model is reduced to what reaches memory: one transition grid per tree (production holds one
// per clique, in a contiguous table), random fixed bins, and data factors of 1. The block update's
// loop, the eight-state draw and the walks are the production ones.
//

#include "cli.h"
#include "constants.h"
#include "coretools/Types/probability.h"
#include "field/TBlockUpdate.h"
#include "field/TFieldMath.h"
#include "random/TCellUniforms.h"
#include "storages/y_storage/TStorageY.h"
#include "storages/y_storage/TStorageYDense.h"
#include "storages/z_storage/TStorageZDense.h"
#include "tree/TPhylogeny.h"
#include "tree/branch/TBinGrid.h"
#include "tree/branch/TTransitionGrid.h"
#include "tree/clique/TCliqueSpace.h"
#include "tree/clique/TCliqueView.h"
#include "tree/clique/clique_traversal.h"
#include "tree/node_state_density.h"
#include "tree/node_state_walk.h"

#include <armadillo>
#include <omp.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Policy = field_math::TAndLinkPolicy;

constexpr uint64_t SEED  = 20261020;
constexpr double ALPHA   = 0.05;
constexpr double NU      = 1.0;
constexpr double OMEGA   = 0.01;
constexpr double MIB     = 1024.0 * 1024.0;
constexpr size_t SPECIES = 0;
constexpr size_t MOLECULE = 1;

// The shapes of the LOTUS run in inference_lotus/ (leaves, internal nodes including roots, roots).
// `--scale` multiplies the first two.
constexpr std::array<size_t, 3> LOTUS_SPECIES  = {32347, 10832, 3};
constexpr std::array<size_t, 3> LOTUS_MOLECULE = {141618, 752, 7};

// -------------------------------------------------------------------------------------------------
// Options
// -------------------------------------------------------------------------------------------------

struct TOptions {
	std::string species_tree;
	std::string molecule_tree;
	double scale        = 0.1;
	size_t iterations   = 3;
	size_t warmup       = 1;
	size_t threads      = 0;
	size_t n_bins       = 10;
	std::vector<std::string> variants = {"own", "own_clique", "arma", "arma_clique"};
};

[[noreturn]] void usage(const std::string &problem) {
	std::cerr << problem << "\n\n"
	          << "acol_dense_bench [--species FILE --molecules FILE | --scale X] [--iterations N]\n"
	          << "                 [--warmup N] [--threads N] [--n_bins N]\n"
	          << "                 [--variants own,own_clique,arma,arma_clique]\n\n"
	          << "  --species/--molecules  real tree files (child, parent, branch_length), as acol\n"
	          << "                         infer reads them; both or neither\n"
	          << "  --scale                synthetic trees shaped like inference_lotus/, leaf and\n"
	          << "                         internal counts multiplied by X (default 0.1; 1 = full)\n"
	          << "  --iterations           timed iterations (default 3)\n"
	          << "  --warmup               untimed iterations before them (default 1)\n"
	          << "  --threads              OpenMP threads (default: all)\n";
	std::exit(2);
}

TOptions parse(int argc, char **argv) {
	TOptions options;
	for (int i = 1; i < argc; ++i) {
		const std::string flag = argv[i];
		if (flag == "--help" || flag == "-h") { usage("Benchmarks the dense storages."); }
		if (i + 1 >= argc) { usage("Missing a value after " + flag); }
		const std::string value = argv[++i];
		if (flag == "--species") {
			options.species_tree = value;
		} else if (flag == "--molecules") {
			options.molecule_tree = value;
		} else if (flag == "--scale") {
			options.scale = std::stod(value);
		} else if (flag == "--iterations") {
			options.iterations = std::stoul(value);
		} else if (flag == "--warmup") {
			options.warmup = std::stoul(value);
		} else if (flag == "--threads") {
			options.threads = std::stoul(value);
		} else if (flag == "--n_bins") {
			options.n_bins = std::stoul(value);
		} else if (flag == "--variants") {
			options.variants.clear();
			std::stringstream list(value);
			for (std::string name; std::getline(list, name, ',');) {
				if (name != "own" && name != "own_clique" && name != "arma" &&
				    name != "arma_clique") {
					usage("Unknown variant " + name);
				}
				options.variants.push_back(name);
			}
		} else {
			usage("Unknown flag " + flag);
		}
	}
	if (options.species_tree.empty() != options.molecule_tree.empty()) {
		usage("Give both --species and --molecules, or neither.");
	}
	if (options.iterations == 0) { usage("--iterations must be at least 1."); }
	return options;
}

// -------------------------------------------------------------------------------------------------
// Trees and the reduced model
// -------------------------------------------------------------------------------------------------

/// A forest shaped like the LOTUS trees: `n_internal` internal nodes of which the first `n_roots`
/// are roots, each other internal node hung below an earlier one, and the leaves hung below the
/// internal nodes. Every internal node gets one leaf first, so none is childless. The order the
/// leaves are listed in does not reach the benchmark: canonical order numbers them depth-first,
/// siblings adjacent, whatever the file says (ADR-0014).
TPhylogeny synthetic_tree(size_t n_leaves, size_t n_internal, size_t n_roots, uint64_t seed) {
	n_roots    = std::max<size_t>(1, n_roots);
	n_internal = std::max(n_internal, n_roots);
	if (n_leaves < n_internal) {
		throw std::invalid_argument("A synthetic tree needs at least one leaf per internal node.");
	}
	std::mt19937_64 rng(seed);
	std::vector<TEdge> edges;
	edges.reserve(n_leaves + n_internal);
	auto internal = [](size_t i) { return "i" + std::to_string(i); };
	for (size_t i = n_roots; i < n_internal; ++i) {
		std::uniform_int_distribution<size_t> pick(0, i - 1);
		edges.push_back({internal(i), internal(pick(rng)), 1.0});
	}
	std::uniform_int_distribution<size_t> any_internal(0, n_internal - 1);
	for (size_t leaf = 0; leaf < n_leaves; ++leaf) {
		const size_t parent = leaf < n_internal ? leaf : any_internal(rng);
		edges.push_back({"l" + std::to_string(leaf), internal(parent), 1.0});
	}
	return build_phylogeny(edges);
}

/// One tree's process: a single transition grid, and a fixed random bin per node.
struct TProcess {
	TTransitionGrid grid;
	std::vector<size_t> bins;

	TProcess(const TPhylogeny &tree, const TBinGrid &bin_grid, uint64_t seed)
	    : grid(ALPHA, NU, bin_grid), bins(tree.n_nodes()) {
		std::mt19937_64 rng(seed);
		std::uniform_int_distribution<size_t> pick(0, bin_grid.n_bins() - 1);
		for (auto &bin : bins) { bin = pick(rng); }
	}

	[[nodiscard]] auto bin_of() const {
		return [this](size_t node) { return bins[node]; };
	}
};

struct TModelInputs {
	const TPhylogeny &species;
	const TPhylogeny &molecule;
	const TProcess &species_process;
	const TProcess &molecule_process;
};

/// The block update's model, reduced to the two transition grids. Every data factor is 1.
class TBenchBlockModel {
	const TModelInputs *_in;

public:
	explicit TBenchBlockModel(const TModelInputs &in) : _in(&in) {}

	struct TRowState {
		size_t species_bin = 0;
	};

	[[nodiscard]] TRowState begin_row(size_t species_leaf, size_t /*first_field_cell*/) const {
		return {_in->species_process.bins[species_leaf]};
	}

	struct TColumnState {
		size_t molecule_bin = 0;
	};

	[[nodiscard]] TColumnState begin_column(size_t molecule_leaf) const {
		return {_in->molecule_process.bins[molecule_leaf]};
	}

	[[nodiscard]] block_update::TLeafPairFactors factors(const TRowState &row,
	                                                     const TColumnState &column,
	                                                     size_t /*molecule_leaf*/,
	                                                     bool species_parent,
	                                                     bool molecule_parent) const {
		const coretools::Probability one = coretools::P(1.0);
		return {.prob_z_s_is_one = coretools::P(_in->species_process.grid.probability(
		            row.species_bin, species_parent, true)),
		        .prob_z_m_is_one = coretools::P(_in->molecule_process.grid.probability(
		            column.molecule_bin, molecule_parent, true)),
		        .lotus           = {one, one},
		        .simple_error    = {one, one},
		        .mass_spec       = {one, one}};
	}

	void record(TRowState & /*row*/, size_t /*molecule_leaf*/,
	            const block_update::TLeafPairFactors & /*factors*/,
	            const field_math::TBlockStates & /*drawn*/) const {}
	void end_row(TRowState & /*row*/) const {}
	void prepare_for_traversal() const {}
};

static_assert(block_update::BlockModel<TBenchBlockModel>);

// -------------------------------------------------------------------------------------------------
// Measurement
// -------------------------------------------------------------------------------------------------

/// The resident set of this process, from /proc/self/statm.
double resident_mib() {
	std::ifstream statm("/proc/self/statm");
	size_t total_pages    = 0;
	size_t resident_pages = 0;
	statm >> total_pages >> resident_pages;
	return static_cast<double>(resident_pages) * static_cast<double>(sysconf(_SC_PAGESIZE)) / MIB;
}

using Clock = std::chrono::steady_clock;

template<typename Fn> double time_ms(Fn &&fn) {
	const auto start = Clock::now();
	fn();
	return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

/// An order-independent fingerprint of a container: how many cells hold something, and a sum of
/// the hashes of (row-major position, value). Two backings agree on it when they hold the same
/// cells, whichever order their memory is walked in.
struct TFingerprint {
	uint64_t n_nonzero = 0;
	uint64_t hash      = 0;

	static uint64_t mix(uint64_t z) {
		z ^= z >> 30U;
		z *= 0xBF58476D1CE4E5B9ULL;
		z ^= z >> 27U;
		z *= 0x94D049BB133111EBULL;
		z ^= z >> 31U;
		return z;
	}

	bool operator==(const TFingerprint &) const = default;
};

std::string to_string(const TFingerprint &f) {
	std::ostringstream out;
	out << f.n_nonzero << " / " << std::hex << f.hash;
	return out.str();
}

/// What one variant measured, and the state it ended in.
struct TResult {
	std::string name;
	double nominal_mib   = 0.0; ///< the bytes the three containers allocate
	double resident_mib  = 0.0; ///< the growth of the resident set once they are built and run
	std::map<std::string, std::vector<double>> pass_ms{};
	TFingerprint y{}, z_species{}, z_molecule{};
	std::array<std::array<size_t, 2>, 3> link_counters{};
	double species_density  = 0.0;
	double molecule_density = 0.0;
};

const std::vector<std::string> PASSES = {"block update",   "species walk",     "species density",
                                         "molecule walk",  "molecule density", "posterior counter"};

std::array<std::array<size_t, 2>, 3> counts_of(const field_math::TLinkCounters &counters) {
	std::array<std::array<size_t, 2>, 3> out{};
	for (size_t bucket = 0; bucket < 3; ++bucket) {
		for (size_t y = 0; y < 2; ++y) { out[bucket][y] = counters.count(bucket, y != 0); }
	}
	return out;
}

/// The shared driver: one iteration is every pass once, in the order a chain runs them. The
/// variant supplies each pass, and what is left is the bookkeeping.
struct TPasses {
	std::function<field_math::TLinkCounters(size_t)> block;
	std::function<void(size_t)> species_walk;
	std::function<double()> species_density;
	std::function<void(size_t)> molecule_walk;
	std::function<double()> molecule_density;
	std::function<void(size_t)> counter;
};

void run_chain(const TOptions &options, TPasses &passes, TResult &result) {
	for (size_t iteration = 0; iteration < options.warmup + options.iterations; ++iteration) {
		const bool timed = iteration >= options.warmup;
		std::map<std::string, double> ms;
		field_math::TLinkCounters counters;
		ms["block update"]      = time_ms([&] { counters = passes.block(iteration); });
		ms["species walk"]      = time_ms([&] { passes.species_walk(iteration); });
		ms["species density"]   = time_ms([&] { result.species_density = passes.species_density(); });
		ms["molecule walk"]     = time_ms([&] { passes.molecule_walk(iteration); });
		ms["molecule density"]  = time_ms([&] { result.molecule_density = passes.molecule_density(); });
		ms["posterior counter"] = time_ms([&] { passes.counter(iteration); });
		result.link_counters    = counts_of(counters);

		double total = 0.0;
		std::cout << "    iteration " << iteration << (timed ? "" : " (warm-up)") << ":";
		for (const auto &pass : PASSES) {
			total += ms[pass];
			if (timed) { result.pass_ms[pass].push_back(ms[pass]); }
		}
		std::printf(" %.0f ms\n", total);
	}
}

/// Every density is summed in clique order, so it is the same number under every backing.
double sum_in_order(const std::vector<double> &per_clique) {
	double sum = 0.0;
	for (const double d : per_clique) { sum += d; }
	return sum;
}

// -------------------------------------------------------------------------------------------------
// own: the project's dense storages through the production code paths
// -------------------------------------------------------------------------------------------------

/// `own` sizes both node states row-major. `own_clique` sizes each the way TTree sizes its own,
/// with the tree's own dimension contiguous (ADR-0013).
TResult run_own(const std::string &name, const TOptions &options, const TModelInputs &in,
                bool clique_contiguous) {
	TResult result{.name = name};
	const size_t n_s = in.species.n_leaves();
	const size_t n_m = in.molecule.n_leaves();

	const double before = resident_mib();
	TStorageYDense Y(options.warmup + options.iterations, IndexArray{n_s, n_m});
	TStorageZDense Z_s(IndexArray{in.species.n_nodes(), n_m}, clique_contiguous ? SPECIES : 1);
	TStorageZDense Z_m(IndexArray{n_s, in.molecule.n_nodes()}, MOLECULE);
	result.nominal_mib = static_cast<double>(Y.total_size_of_container_space() * sizeof(TStorageY) +
	                                         Z_s.total_size_of_container_space() +
	                                         Z_m.total_size_of_container_space()) /
	                     MIB;

	const TCliqueSpace<> species_cliques({n_s, n_m}, SPECIES);
	const TCliqueSpace<> molecule_cliques({n_s, n_m}, MOLECULE);
	TBenchBlockModel model(in);
	const field_math::TErrorProbability omega(OMEGA);

	auto species_view = [&](size_t c) {
		return TCliqueView<TStorageZDense>(Z_s, in.species, species_cliques.index_of(c), SPECIES);
	};
	auto molecule_view = [&](size_t c) {
		return TCliqueView<TStorageZDense>(Z_m, in.molecule, molecule_cliques.index_of(c),
		                                   MOLECULE);
	};

	TPasses passes;
	passes.block = [&](size_t iteration) {
		const TCellUniforms uniforms(SEED, TCellStream::field, iteration);
		return block_update::run<Policy>(Y, Z_s, Z_m, in.species, in.molecule, omega, model,
		                                 uniforms);
	};
	auto walk = [&](auto &view_factory, const TPhylogeny &tree, const TProcess &process,
	                size_t n_cliques, size_t dimension, size_t iteration) {
		const TCellUniforms uniforms(SEED, TCellStream::node_state, iteration, dimension);
		auto bins   = process.bin_of();
		auto kernel = [&](size_t /*c*/, TCliqueView<TStorageZDense> &view) {
			node_state_walk::update_clique(tree, process.grid, bins, uniforms, view);
		};
		(void)clique_traversal::run(n_cliques, view_factory, kernel);
	};
	auto density = [&](auto &view_factory, const TPhylogeny &tree, const TProcess &process,
	                   size_t n_cliques) {
		auto bins   = process.bin_of();
		auto kernel = [&](size_t /*c*/, TCliqueView<TStorageZDense> &view) {
			return node_state_density::log_density_of_clique(tree, process.grid, view, bins);
		};
		return sum_in_order(clique_traversal::run(n_cliques, view_factory, kernel).first);
	};
	passes.species_walk = [&](size_t iteration) {
		walk(species_view, in.species, in.species_process, n_m, SPECIES, iteration);
	};
	passes.species_density = [&] {
		return density(species_view, in.species, in.species_process, n_m);
	};
	passes.molecule_walk = [&](size_t iteration) {
		walk(molecule_view, in.molecule, in.molecule_process, n_s, MOLECULE, iteration);
	};
	passes.molecule_density = [&] {
		return density(molecule_view, in.molecule, in.molecule_process, n_s);
	};
	passes.counter = [&](size_t iteration) { Y.add_to_counter(iteration); };

	run_chain(options, passes, result);
	result.resident_mib = resident_mib() - before;

	// The fingerprints, by row-major linear index -- which names a cell under every layout.
	auto fingerprint_z = [](const TStorageZDense &Z) {
		TFingerprint f;
		uint64_t n = 0, h = 0;
		const size_t total = Z.total_size_of_container_space();
#pragma omp parallel for reduction(+ : n, h)
		for (size_t i = 0; i < total; ++i) {
			if (Z.is_one(i)) {
				++n;
				h += TFingerprint::mix(i);
			}
		}
		f.n_nonzero = n;
		f.hash      = h;
		return f;
	};
	{
		uint64_t n = 0, h = 0;
		const size_t total = Y.total_size_of_container_space();
#pragma omp parallel for reduction(+ : n, h)
		for (size_t i = 0; i < total; ++i) {
			const uint16_t value = Y[i].value();
			if (value != 0) {
				++n;
				h += TFingerprint::mix(i + (static_cast<uint64_t>(value) << 40U));
			}
		}
		result.y = {n, h};
	}
	result.z_species  = fingerprint_z(Z_s);
	result.z_molecule = fingerprint_z(Z_m);
	return result;
}

// -------------------------------------------------------------------------------------------------
// arma and arma_clique: arma::Mat, in either orientation
// -------------------------------------------------------------------------------------------------

/// A logical (row, column) matrix held in an arma::Mat. `Transposed` stores the logical (r, c) at
/// the matrix's (c, r), which in column-major order makes a logical row contiguous.
///
/// `name` is the row-major linear index of the logical cell: what `own` holds the cell at, and so
/// what its uniform is hashed from.
template<typename Elem, bool Transposed> class TArmaGrid {
	arma::Mat<Elem> _m;
	size_t _n_cols = 0;

public:
	TArmaGrid(size_t n_rows, size_t n_cols)
	    : _m(Transposed ? n_cols : n_rows, Transposed ? n_rows : n_cols, arma::fill::zeros),
	      _n_cols(n_cols) {}

	[[nodiscard]] Elem &at(size_t row, size_t col) {
		if constexpr (Transposed) {
			return _m.at(col, row);
		} else {
			return _m.at(row, col);
		}
	}
	[[nodiscard]] Elem at(size_t row, size_t col) const {
		if constexpr (Transposed) {
			return _m.at(col, row);
		} else {
			return _m.at(row, col);
		}
	}
	[[nodiscard]] size_t name(size_t row, size_t col) const { return row * _n_cols + col; }
	[[nodiscard]] size_t bytes() const { return _m.n_elem * sizeof(Elem); }
	[[nodiscard]] arma::Mat<Elem> &mat() { return _m; }
	[[nodiscard]] const arma::Mat<Elem> &mat() const { return _m; }

	/// The logical position of the element at the matrix's own (i, j).
	[[nodiscard]] size_t name_of_element(size_t i, size_t j) const {
		return Transposed ? name(j, i) : name(i, j);
	}
};

constexpr uint16_t Y_STATE   = 0x8000;
constexpr uint16_t Y_COUNTER = 0x7FFF;

/// One clique of a node state held in a TArmaGrid, as the walk and the density read it.
/// `NodeIsRow`: the tree's node index is the logical row (the species tree) or column (molecule).
template<typename Grid, bool NodeIsRow> class TArmaCliqueColumn {
	Grid *_z;
	size_t _fixed;

public:
	TArmaCliqueColumn(Grid &z, size_t fixed) : _z(&z), _fixed(fixed) {}

	[[nodiscard]] bool is_one(size_t node) const {
		return (NodeIsRow ? _z->at(node, _fixed) : _z->at(_fixed, node)) != 0;
	}
	[[nodiscard]] size_t linear_index(size_t node) const {
		return NodeIsRow ? _z->name(node, _fixed) : _z->name(_fixed, node);
	}
	/// Writes only a change, as TCliqueView does: an unconditional write dirties the line for
	/// nothing, and a dirty line shared with another thread's clique is taken from it.
	void set_state(size_t node, bool state) {
		uint8_t &cell = NodeIsRow ? _z->at(node, _fixed) : _z->at(_fixed, node);
		if ((cell != 0) != state) { cell = state ? 1 : 0; }
	}
	[[nodiscard]] std::vector<size_t> take_deferred_inserts() { return {}; }
};

template<bool YT, bool ZsT, bool ZmT>
TResult run_arma(const std::string &name, const TOptions &options, const TModelInputs &in) {
	TResult result{.name = name};
	const size_t n_s = in.species.n_leaves();
	const size_t n_m = in.molecule.n_leaves();

	using TY  = TArmaGrid<uint16_t, YT>;
	using TZs = TArmaGrid<uint8_t, ZsT>;
	using TZm = TArmaGrid<uint8_t, ZmT>;

	const double before = resident_mib();
	TY Y(n_s, n_m);
	TZs Z_s(in.species.n_nodes(), n_m);
	TZm Z_m(n_s, in.molecule.n_nodes());
	result.nominal_mib = static_cast<double>(Y.bytes() + Z_s.bytes() + Z_m.bytes()) / MIB;

	TBenchBlockModel model(in);
	const field_math::TErrorProbability omega(OMEGA);
	const size_t n_threads = ProgramOptions::NUMBER_OF_THREADS;

	TPasses passes;
	// block_update::run, with the storage reads and writes spelled against the matrices. The
	// schedule, the order inside a row, the draw and the tally are the same.
	passes.block = [&](size_t iteration) {
		const TCellUniforms uniforms(SEED, TCellStream::field, iteration);
		const TPhylogeny &species  = in.species;
		const TPhylogeny &molecule = in.molecule;
		std::vector<field_math::TLinkCounters> per_row(n_s);
#pragma omp parallel for num_threads(n_threads) schedule(static)
		for (size_t s = 0; s < n_s; ++s) {
			auto row                     = model.begin_row(s, s * n_m);
			const size_t species_parent  = species.parent_of(s);
			field_math::TLinkCounters counters;
			for (size_t m = 0; m < n_m; ++m) {
				const bool parent_s = Z_s.at(species_parent, m) != 0;
				const bool parent_m = Z_m.at(s, molecule.parent_of(m)) != 0;
				const auto factors =
				    model.factors(row, model.begin_column(m), m, parent_s, parent_m);

				uint16_t &y  = Y.at(s, m);
				uint8_t &z_s = Z_s.at(s, m);
				uint8_t &z_m = Z_m.at(s, m);
				const field_math::TBlockStates current{
				    .y = (y & Y_STATE) != 0, .z_s = z_s != 0, .z_m = z_m != 0};
				const auto draw = field_math::draw_block<Policy>(
				    factors.prob_z_s_is_one, factors.prob_z_m_is_one, omega, factors.lotus,
				    factors.simple_error, factors.mass_spec, current, uniforms.at(Y.name(s, m)));
				if (current.y != draw.drawn.y) {
					y = draw.drawn.y ? static_cast<uint16_t>(y | Y_STATE)
					                 : static_cast<uint16_t>(y & Y_COUNTER);
				}
				if (current.z_s != draw.drawn.z_s) { z_s = draw.drawn.z_s ? 1 : 0; }
				if (current.z_m != draw.drawn.z_m) { z_m = draw.drawn.z_m ? 1 : 0; }
				counters.add(draw.to.bucket, draw.to.y);
				model.record(row, m, factors, draw.drawn);
			}
			per_row[s].merge(counters);
			model.end_row(row);
		}
		field_math::TLinkCounters total;
		for (const auto &row : per_row) { total.merge(row); }
		return total;
	};

	auto species_view  = [&](size_t c) { return TArmaCliqueColumn<TZs, true>(Z_s, c); };
	auto molecule_view = [&](size_t c) { return TArmaCliqueColumn<TZm, false>(Z_m, c); };
	auto walk = [&](auto &view_factory, const TPhylogeny &tree, const TProcess &process,
	                size_t n_cliques, size_t dimension, size_t iteration) {
		const TCellUniforms uniforms(SEED, TCellStream::node_state, iteration, dimension);
		auto bins   = process.bin_of();
		auto kernel = [&](size_t /*c*/, auto &view) {
			node_state_walk::update_clique(tree, process.grid, bins, uniforms, view);
		};
		(void)clique_traversal::run(n_cliques, view_factory, kernel);
	};
	auto density = [&](auto &view_factory, const TPhylogeny &tree, const TProcess &process,
	                   size_t n_cliques) {
		auto bins   = process.bin_of();
		auto kernel = [&](size_t /*c*/, auto &view) {
			return node_state_density::log_density_of_clique(tree, process.grid, view, bins);
		};
		return sum_in_order(clique_traversal::run(n_cliques, view_factory, kernel).first);
	};
	passes.species_walk = [&](size_t iteration) {
		walk(species_view, in.species, in.species_process, n_m, SPECIES, iteration);
	};
	passes.species_density = [&] {
		return density(species_view, in.species, in.species_process, n_m);
	};
	passes.molecule_walk = [&](size_t iteration) {
		walk(molecule_view, in.molecule, in.molecule_process, n_s, MOLECULE, iteration);
	};
	passes.molecule_density = [&] {
		return density(molecule_view, in.molecule, in.molecule_process, n_s);
	};
	// TStorageYDense::add_to_counter, over the matrix's memory. The thinning factor is 1 for any
	// chain this benchmark runs, so every iteration counts.
	passes.counter = [&](size_t /*iteration*/) {
		uint16_t *cells    = Y.mat().memptr();
		const size_t total = Y.mat().n_elem;
#pragma omp parallel for num_threads(n_threads)
		for (size_t i = 0; i < total; ++i) {
			if ((cells[i] & Y_STATE) != 0) { cells[i] = static_cast<uint16_t>(cells[i] + 1); }
		}
	};

	run_chain(options, passes, result);
	result.resident_mib = resident_mib() - before;

	auto fingerprint = [](const auto &grid, bool with_value) {
		const auto &m = grid.mat();
		uint64_t n = 0, h = 0;
		const size_t n_rows = m.n_rows;
		const size_t n_cols = m.n_cols;
#pragma omp parallel for reduction(+ : n, h)
		for (size_t j = 0; j < n_cols; ++j) {
			for (size_t i = 0; i < n_rows; ++i) {
				const uint64_t value = m.at(i, j);
				if (value == 0) { continue; }
				++n;
				const uint64_t position = grid.name_of_element(i, j);
				h += TFingerprint::mix(with_value ? position + (value << 40U) : position);
			}
		}
		return TFingerprint{n, h};
	};
	result.y          = fingerprint(Y, true);
	result.z_species  = fingerprint(Z_s, false);
	result.z_molecule = fingerprint(Z_m, false);
	return result;
}

// -------------------------------------------------------------------------------------------------
// Report
// -------------------------------------------------------------------------------------------------

double mean(const std::vector<double> &xs) {
	double sum = 0.0;
	for (const double x : xs) { sum += x; }
	return sum / static_cast<double>(xs.size());
}

void report(const std::vector<TResult> &results) {
	std::printf("\n%-28s", "mean over timed iterations");
	for (const auto &r : results) { std::printf("%16s", r.name.c_str()); }
	std::printf("\n");
	std::vector<double> totals(results.size(), 0.0);
	for (const auto &pass : PASSES) {
		std::printf("%-28s", (pass + " [ms]").c_str());
		for (size_t i = 0; i < results.size(); ++i) {
			const double ms = mean(results[i].pass_ms.at(pass));
			totals[i] += ms;
			std::printf("%16.1f", ms);
		}
		std::printf("\n");
	}
	std::printf("%-28s", "iteration [ms]");
	for (const double t : totals) { std::printf("%16.1f", t); }
	std::printf("\n%-28s", "allocated [MiB]");
	for (const auto &r : results) { std::printf("%16.1f", r.nominal_mib); }
	std::printf("\n%-28s", "resident growth [MiB]");
	for (const auto &r : results) { std::printf("%16.1f", r.resident_mib); }
	std::printf("\n");
}

/// Every variant ran one chain, so every variant ends in one state. A difference means a backing
/// read or wrote a cell other than the one it named, and its timings measured something else.
bool same_chain(const std::vector<TResult> &results) {
	bool same = true;
	const TResult &first = results.front();
	std::printf("\nfinal state (non-zero cells / hash):\n");
	for (const auto &r : results) {
		std::printf("  %-12s Y %s  Z_s %s  Z_m %s  log p(Z_s) %.6f  log p(Z_m) %.6f\n",
		            r.name.c_str(), to_string(r.y).c_str(), to_string(r.z_species).c_str(),
		            to_string(r.z_molecule).c_str(), r.species_density, r.molecule_density);
		same = same && r.y == first.y && r.z_species == first.z_species &&
		       r.z_molecule == first.z_molecule && r.link_counters == first.link_counters &&
		       r.species_density == first.species_density &&
		       r.molecule_density == first.molecule_density;
	}
	std::printf(same ? "  every variant ran the same chain\n"
	                 : "  MISMATCH: the variants did not run the same chain\n");
	return same;
}

} // namespace

int main(int argc, char **argv) {
	const TOptions options = parse(argc, argv);
	ProgramOptions::NUMBER_OF_THREADS =
	    options.threads == 0 ? static_cast<size_t>(omp_get_max_threads()) : options.threads;
	omp_set_num_threads(static_cast<int>(ProgramOptions::NUMBER_OF_THREADS));

	auto scaled = [&](size_t n) {
		return std::max<size_t>(1, static_cast<size_t>(std::llround(options.scale *
		                                                            static_cast<double>(n))));
	};
	const bool real = !options.species_tree.empty();
	const TPhylogeny species =
	    real ? read_phylogeny(options.species_tree)
	         : synthetic_tree(scaled(LOTUS_SPECIES[0]), scaled(LOTUS_SPECIES[1]), LOTUS_SPECIES[2],
	                          SEED + 1);
	const TPhylogeny molecule =
	    real ? read_phylogeny(options.molecule_tree)
	         : synthetic_tree(scaled(LOTUS_MOLECULE[0]), scaled(LOTUS_MOLECULE[1]),
	                          LOTUS_MOLECULE[2], SEED + 2);

	const TBinGrid bin_grid(options.n_bins);
	const TProcess species_process(species, bin_grid, SEED + 3);
	const TProcess molecule_process(molecule, bin_grid, SEED + 4);
	const TModelInputs inputs{species, molecule, species_process, molecule_process};

	std::printf("trees: %s\n", real ? "from file" : "synthetic, LOTUS-shaped");
	std::printf("  species   %zu nodes, %zu leaves, %zu roots\n", species.n_nodes(),
	            species.n_leaves(), species.n_roots());
	std::printf("  molecules %zu nodes, %zu leaves, %zu roots\n", molecule.n_nodes(),
	            molecule.n_leaves(), molecule.n_roots());
	std::printf("field %zu x %zu = %.3g cells; %zu threads; %zu warm-up + %zu timed iterations\n\n",
	            species.n_leaves(), molecule.n_leaves(),
	            static_cast<double>(species.n_leaves()) * static_cast<double>(molecule.n_leaves()),
	            ProgramOptions::NUMBER_OF_THREADS, options.warmup, options.iterations);

	std::vector<TResult> results;
	for (const auto &variant : options.variants) {
		std::printf("%s\n", variant.c_str());
		std::fflush(stdout);
		if (variant == "own" || variant == "own_clique") {
			results.push_back(run_own(variant, options, inputs, variant == "own_clique"));
		} else if (variant == "arma") {
			results.push_back(run_arma<false, false, false>(variant, options, inputs));
		} else {
			results.push_back(run_arma<true, false, true>(variant, options, inputs));
		}
	}

	report(results);
	return same_chain(results) ? 0 : 1;
}
