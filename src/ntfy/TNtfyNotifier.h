#pragma once

#include <string>
#include <vector>
#define CPPHTTPLIB_OPENSSL_SUPPORT

class TNtfyNotifier {
private:
	std::string _topic;

	void _send(const std::string &title, const std::string &message, int priority = 3,
	           const std::string &tag = "") const;

public:
	TNtfyNotifier();

	[[nodiscard]] std::string get_topic_url() const;

	void notify_start(const std::vector<std::string> &tree_names,
	                  const std::vector<size_t> &leaf_counts, size_t n_iterations, size_t n_burnin,
	                  size_t n_burnin_iters) const;

	struct ParamStats {
		double mean;
		double var;
		double sd;
	};

	void notify_burnin_round(size_t round, size_t total_rounds,
	                         const std::vector<std::string> &dim_names,
	                         const std::vector<ParamStats> &gamma_stats,
	                         const ParamStats &epsilon_stats) const;

	void notify_burnin_finished(const std::vector<std::string> &dim_names,
	                            const std::vector<ParamStats> &gamma_stats,
	                            const ParamStats &epsilon_stats) const;

	void notify_mcmc_finished(const std::vector<std::string> &dim_names,
	                          const std::vector<ParamStats> &gamma_stats,
	                          const ParamStats &epsilon_stats) const;
};
