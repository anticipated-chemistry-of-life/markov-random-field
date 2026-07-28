#pragma once

#include <string>
#include <vector>
#define CPPHTTPLIB_OPENSSL_SUPPORT

class TNtfyNotifier {
public:
	struct ParamStats {
		double mean;
		double var;
		double sd;
	};

	/// A scalar parameter reported under its own name. Which scalars exist depends on the data
	/// sources compiled in (epsilon comes with LOTUS, epsilon_simple_model with the simple error
	/// model), so they are reported as a named list rather than as fixed arguments.
	struct NamedStats {
		std::string name;
		ParamStats stats;
	};

private:
	std::string _topic;

	void _send(const std::string &title, const std::string &message, int priority = 3,
	           const std::string &tag = "") const;

	/// Body shared by the three parameter notifications: one line per kept dimension of gamma,
	/// then one line per scalar. Empty vectors simply produce no lines, which is what a build
	/// without LOTUS (no gamma, no epsilon) reports.
	[[nodiscard]] static std::string
	_format_param_stats(const std::vector<std::string> &dim_names,
	                    const std::vector<ParamStats> &gamma_stats,
	                    const std::vector<NamedStats> &scalar_stats);

public:
	TNtfyNotifier();

	[[nodiscard]] std::string get_topic_url() const;

	void notify_start(const std::vector<std::string> &tree_names,
	                  const std::vector<size_t> &leaf_counts, size_t n_iterations, size_t n_burnin,
	                  size_t n_burnin_iters) const;

	void notify_burnin_round(size_t round, size_t total_rounds,
	                         const std::vector<std::string> &dim_names,
	                         const std::vector<ParamStats> &gamma_stats,
	                         const std::vector<NamedStats> &scalar_stats) const;

	void notify_burnin_finished(const std::vector<std::string> &dim_names,
	                            const std::vector<ParamStats> &gamma_stats,
	                            const std::vector<NamedStats> &scalar_stats) const;

	void notify_mcmc_finished(const std::vector<std::string> &dim_names,
	                          const std::vector<ParamStats> &gamma_stats,
	                          const std::vector<NamedStats> &scalar_stats) const;
};
