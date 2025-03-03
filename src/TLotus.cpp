//
// Created by madleina on 03.03.25.
//
#include "TLotus.h"

TLotus::TLotus(
    std::vector<std::unique_ptr<TTree>> &trees, TypeParamGamma *gamma, size_t n_iterations,
    const std::vector<std::unique_ptr<stattools::TParameter<SpecMarkovField, TLotus>>> &markov_field_stattools_param,
    std::string prefix)
    : _trees(trees), _markov_field(n_iterations, trees), _markov_field_stattools_param(markov_field_stattools_param),
      _collapser(trees), _gamma(gamma), _tmp_state_along_last_dim(*trees.back().get(), 1), _prefix(std::move(prefix)) {
	this->addPriorParameter(_gamma);
	for (auto &it : _markov_field_stattools_param) { this->addPriorParameter(it.get()); }

	_oldLL = 0.0;
	_curLL = 0.0;
}

[[nodiscard]] std::string TLotus::name() const { return "lotus_likelihood"; }

void TLotus::initialize() {
	// initialize storage
	_gamma->initStorage(this, {_collapser.num_dim_to_keep()});
	if constexpr (UseSimpleErrorModel) {
		_epsilon = coretools::instances::parameters().get<double>("epsilon", 0.0001);
		coretools::instances::logfile().list("Using simple error model with epsilon = ", _epsilon);
	}
}

void TLotus::load_from_file(const std::string &filename) {
	coretools::instances::logfile().startIndent("Reading links from file '", filename, "' ...");
	coretools::TInputFile file(filename, coretools::FileType::Header);

	// initialize collapser: know which dimensions to keep and which to collapse
	std::vector<std::string> tree_names(_trees.size());
	for (const auto &tree : _trees) { tree_names.push_back(tree->get_tree_name()); }

	// check if all headers in file match a tree name
	for (const auto &header_name : file.header()) {
		if (std::find(tree_names.begin(), tree_names.end(), header_name) == tree_names.end()) {
			UERROR("Header '", header_name, "' in file '", filename, "' does not match any tree name!");
		}
	}

	// since lotus can also contain less dimensions than the trees, we need to check if the order of the headers
	// matches the order of the trees. To do so we can do a dequing of the tree names and the headers and check if
	// the headers are a subset of the tree names.
	std::deque<std::string> header_names_deque(file.header().size());
	for (const auto &header_name : file.header()) { header_names_deque.push_back(header_name); }
	for (const auto &tree_name : tree_names) {
		if (header_names_deque.empty()) { break; }
		if (header_names_deque.front() == tree_name) { header_names_deque.pop_front(); }
	}
	if (!header_names_deque.empty()) {
		UERROR("Headers in file '", filename, "' should be ordered in the same way as the trees.");
	}

	const auto len_per_dimension_lotus = _collapser.initialize(file.header(), "LOTUS");

	// initialize the size of L
	_L.initialize(0, len_per_dimension_lotus);

	_occurrence_counters.resize(_collapser.num_dim_to_keep()); // for example, size is 2 if keep molecules and species
	for (size_t i = 0; i < _collapser.num_dim_to_keep(); ++i) {
		_occurrence_counters[i].resize(_trees[_collapser.dim_to_keep(i)]->get_number_of_leaves(), 0);
	}

	std::vector<size_t> index_in_collapsed_space(_collapser.num_dim_to_keep());
	for (; !file.empty(); file.popFront()) {
		// loop over all columns
		for (size_t i = 0; i < _collapser.num_dim_to_keep(); ++i) {
			std::string node_name = std::string(file.get(i));

			const size_t tree_index = _collapser.dim_to_keep(i);
			if (!_trees[tree_index]->isLeaf(_trees[tree_index]->get_node_index(node_name))) {
				UERROR("Node '", node_name, "' in tree '", _trees[tree_index]->get_tree_name(), "' is not a leaf !");
			}

			const size_t ix             = _trees[tree_index]->get_index_within_leaves(node_name);
			index_in_collapsed_space[i] = ix;
			++_occurrence_counters[i][ix];
		}

		size_t linear_index_in_Y_space = _L.get_linear_index_in_container_space(index_in_collapsed_space);
		_L.insert_one(linear_index_in_Y_space);
	}
	coretools::instances::logfile().endIndent();
};

double TLotus::calculate_log_likelihood_of_L() const {
	if (_collapser.do_collapse()) { return _calculate_log_likelihood_of_L_do_collapse(); }
	return _calculate_log_likelihood_of_L_no_collapsing();
};

double TLotus::getSumLogPriorDensity(const Storage &) const { return _curLL; };

void TLotus::fill_tmp_state_along_last_dim(const std::vector<size_t> &start_index_clique_along_last_dim, size_t K) {
	// collapse start_index_in_leaves (this is the index in Y)
	if (_collapser.do_collapse()) {
		_tmp_state_along_last_dim.fill_Y_along_last_dim(_collapser.collapse(start_index_clique_along_last_dim), K, _L);
	} else { // no need to collapse
		_tmp_state_along_last_dim.fill_Y_along_last_dim(start_index_clique_along_last_dim, K, _L);
	}
};

/// This function will be used when we update Y.
void TLotus::calculate_LL_update_Y(const std::vector<size_t> &index_in_leaves_space, bool old_state,
                                   std::array<coretools::TSumLogProbability, 2> &sum_log) const {
	// function gets the old_state and needs to calculate LL for new_state = 0 and 1
	// for state 1, we know that the new x will always be 1 (at least one is a one)
	bool x_is_one_for_Y_0 = false; // Y = 0 -> x = 0 if we don't collapse
	if (_collapser.do_collapse()) { x_is_one_for_Y_0 = _collapser.x_is_one(index_in_leaves_space, old_state); }

	if (x_is_one_for_Y_0) { // Y=0 also results in x=1 (others in clique are a one)
		return;             // x is one for both states (due to collapsing) -> likelihood doesn't matter
	}

	const size_t leaf_index_last_dim    = index_in_leaves_space.back();
	const auto index_in_collapsed_space = _collapser.collapse(index_in_leaves_space);
	// new Y = 0 -> x_is_one_for_Y_0 will always be false here (because of the previous if-statement)
	// new Y = 1 -> x will always be true
	for (size_t i = 0; i < 2; ++i) {
		sum_log[i].add(_calculate_probability_of_L_given_x(i, _tmp_state_along_last_dim.get_Y(leaf_index_last_dim),
		                                                   index_in_collapsed_space));
	}
};

void TLotus::update_markov_field(size_t /*iteration*/) { _markov_field.update(*this); }

[[nodiscard]] double TLotus::calculateLLRatio(TypeParamGamma *, size_t /*Index*/, const Storage &) {
	_oldLL = _curLL;                          // store current likelihood
	_curLL = calculate_log_likelihood_of_L(); // calculate likelihood of new gamma
	return _curLL - _oldLL;
}

void TLotus::updateTempVals(TypeParamGamma *, size_t /*Index*/, bool Accepted) {
	if (!Accepted) {
		_curLL = _oldLL; // reset
	}
}

void TLotus::guessInitialValues() {
	// TODO: think of an initialization scheme
	_gamma->set(0.001);
}

const TStorageYVector &TLotus::get_Lotus() const { return _L; }

double TLotus::_calculate_research_effort(const std::vector<size_t> &index_in_collapsed_space) const {
	double prod = 1.0;
	for (size_t i = 0; i < _collapser.num_dim_to_keep(); ++i) {
		const size_t leaf_index = index_in_collapsed_space[i];
		const auto counts       = (double)_occurrence_counters[i][leaf_index];
		prod *= 1.0 - exp(-(double)_gamma->value(i) * counts);
	}
	return prod;
}

double TLotus::_calculate_probability_of_L_given_x(bool x, bool L,
                                                   const std::vector<size_t> &index_in_collapsed_space) const {
	if constexpr (UseSimpleErrorModel) {
		if (x && L) { return 1 - _epsilon; }
		if (x) { return _epsilon; }
		if (L) { return _epsilon; }
		return 1 - _epsilon;
	} else {
		if (x && L) { return _calculate_research_effort(index_in_collapsed_space); }
		if (x) { return 1.0 - _calculate_research_effort(index_in_collapsed_space); }
		if (L) { return 0.0; }
		return 1.0;
	}
}

double TLotus::_calculate_probability_of_L_given_x(bool x, bool L, size_t linear_index_in_collapsed_space) const {
	auto index_in_L_space = _L.get_multi_dimensional_index(linear_index_in_collapsed_space);
	return _calculate_probability_of_L_given_x(x, L, index_in_L_space);
};

/// Note that this function does not go through the whole vector L, but only through the 1's in L or Y. This is
/// justified because if L and Y are both 0, the probability is 1 and that does not matter.
double TLotus::_calculate_log_likelihood_of_L_no_collapsing() const {

	// Y = L
	size_t index_in_Y = 0;
	size_t index_in_L = 0;

	const auto length_Y = _markov_field.size_Y();
	const auto length_L = _L.size();

	coretools::TSumLogProbability sum_log;

	for (size_t i = 0; i < std::max(length_Y, length_L); ++i) {

		// if we reach the end of vector L but not the end of vector Y, that means that all the other Ls are 0.
		// But Y can still be 1 or 0.
		if (index_in_L == length_L) {
			const auto index_in_L_space = _markov_field.get_Y(index_in_Y).get_linear_index_in_container_space();
			sum_log.add(
			    _calculate_probability_of_L_given_x(_markov_field.get_Y(index_in_Y).is_one(), false, index_in_L_space));
			++index_in_Y;
		} else if (index_in_Y == length_Y) {
			// if we reach the end of vector Y but not the end of vector L, that means that all the other Y are 0.
			const auto index_in_L_space = _L[index_in_L].get_linear_index_in_container_space();
			sum_log.add(_calculate_probability_of_L_given_x(false, _L[index_in_L].is_one(), index_in_L_space));
			++index_in_L;
		} else if (_markov_field.get_Y(index_in_Y).get_linear_index_in_container_space() ==
		           _L[index_in_L].get_linear_index_in_container_space()) {
			const auto index_in_L_space = _L[index_in_L].get_linear_index_in_container_space();
			sum_log.add(_calculate_probability_of_L_given_x(_markov_field.get_Y(index_in_Y).is_one(),
			                                                _L[index_in_L].is_one(), index_in_L_space));
			++index_in_Y;
			++index_in_L;
		} else if (_markov_field.get_Y(index_in_Y).get_linear_index_in_container_space() <
		           _L[index_in_L].get_linear_index_in_container_space()) {
			const auto index_in_L_space = _markov_field.get_Y(index_in_Y).get_linear_index_in_container_space();
			sum_log.add(
			    _calculate_probability_of_L_given_x(_markov_field.get_Y(index_in_Y).is_one(), false, index_in_L_space));
			++index_in_Y;
		} else {
			const auto index_in_L_space = _L[index_in_L].get_linear_index_in_container_space();
			sum_log.add(_calculate_probability_of_L_given_x(false, _L[index_in_L].is_one(), index_in_L_space));
			++index_in_L;
		}
	}
	return sum_log.getSum();
}

double TLotus::_calculate_log_likelihood_of_L_do_collapse() const {
	const auto length_L = _L.size();

	coretools::TSumLogProbability sum_log;

	size_t previous_index_in_L_space = 0;
	for (size_t i = 0; i < length_L; ++i) { // loop over all 1's in L
		const auto linear_index_in_L_space = _L[i].get_linear_index_in_container_space();
		auto multi_dim_index_in_L_space    = _L.get_multi_dimensional_index(linear_index_in_L_space);
		for (size_t j = previous_index_in_L_space; j <= linear_index_in_L_space; ++j) {
			const bool L = (j == linear_index_in_L_space);
			const bool x = _collapser.x_is_one(multi_dim_index_in_L_space);
			sum_log.add(_calculate_probability_of_L_given_x(x, L, multi_dim_index_in_L_space));
		}
		previous_index_in_L_space = linear_index_in_L_space + 1;
	}
	return sum_log.getSum();
};

void TLotus::_simulateUnderPrior(Storage *) {
	// first simulate Markov random field
	_markov_field.simulate(*this);

	// by default, we keep all the trees
	std::vector<std::string> tree_names_to_keep_default;
	for (const auto &tree : _trees) { tree_names_to_keep_default.push_back(tree->get_tree_name()); }

	// else get the tree names to keep from CLI
	std::vector<std::string> tree_names_to_keep =
	    coretools::instances::parameters().get("tree_names_to_keep", tree_names_to_keep_default);

	const auto len_per_dimension_lotus = _collapser.initialize(tree_names_to_keep, "LOTUS");
	// initialize the size of L
	_L.initialize(1, len_per_dimension_lotus);

	std::vector<std::vector<std::string>> node_names;
	bool x;
	for (size_t i = 0; i < _L.total_size_of_container_space(); ++i) {
		auto multi_dim_index_in_L_space = _L.get_multi_dimensional_index(i);
		if (_collapser.do_collapse()) {
			x = _collapser.x_is_one(multi_dim_index_in_L_space);
		} else {
			// When we don't collapse, this means that the dimensions of Y are equal to the dimensions of Lotus.
			// Since we do simulations, speed in not that important. We could simply run a binary search for each
			// i in Y.
			auto result = _markov_field.get_Y_vector().binary_search(i);
			if (result.first) {
				x = _markov_field.get_Y_vector().is_one(result.second);
			} else {
				x = false;
			}
		}
		double proba = _calculate_probability_of_L_given_x(x, true, multi_dim_index_in_L_space);
		coretools::Probability p(proba);
		bool draw = coretools::instances::randomGenerator().pickOneOfTwo(p);
		if (draw) {
			// for each draw we need to get the node name of the leaf in the correct tree and write it a file
			std::vector<std::string> line(_collapser.num_dim_to_keep());
			for (size_t j = 0; j < _collapser.num_dim_to_keep(); ++j) {
				const size_t tree_index         = _collapser.dim_to_keep(j);
				const size_t leaf_index         = multi_dim_index_in_L_space[j];
				const size_t node_index_in_tree = _trees[tree_index]->get_node_index_from_leaf_index(leaf_index);
				line[j]                         = _trees[tree_index]->get_node_id(node_index_in_tree);
			}
			node_names.push_back(line);
		}
	}

	// write to file
	std::string file_name = _prefix + "_lotus.tsv";

	// we get the tree name for the header of the file.
	std::vector<std::string> header(_collapser.num_dim_to_keep());
	for (size_t j = 0; j < _collapser.num_dim_to_keep(); ++j) {
		const size_t tree_index = _collapser.dim_to_keep(j);
		header[j]               = _trees[tree_index]->get_tree_name();
	}

	coretools::TOutputFile file(file_name, header, "\t");
	for (const auto &line : node_names) { file.writeln(line); }
};
