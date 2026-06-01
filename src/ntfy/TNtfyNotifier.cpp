#include "TNtfyNotifier.h"
#include "coretools/Main/TLog.h"

#include "./httplib.h"
#include <nlohmann/json.hpp>
#include <uuid.h>

#include <iomanip>
#include <random>
#include <sstream>

using json = nlohmann::json;

TNtfyNotifier::TNtfyNotifier()
    : _topic([] {
	      std::random_device rd;
	      auto seed = std::array<std::random_device::result_type, std::mt19937::state_size>{};
	      std::generate(seed.begin(), seed.end(), std::ref(rd));
	      std::seed_seq seq(seed.begin(), seed.end());
	      std::mt19937 engine(seq);
	      uuids::uuid_random_generator gen{engine};
	      return uuids::to_string(gen());
      }()) {
	coretools::instances::logfile().list("ntfy.sh topic: ", get_topic_url());
}

std::string TNtfyNotifier::get_topic_url() const { return "https://ntfy.sh/" + _topic; }

void TNtfyNotifier::_send(const std::string &title, const std::string &message, int priority,
                          const std::string &tag) const {
	json payload;
	payload["topic"]    = _topic;
	payload["title"]    = title;
	payload["message"]  = message;
	payload["priority"] = priority;
	if (!tag.empty()) { payload["tags"] = json::array({tag}); }

	try {
		httplib::Client cli("https://ntfy.sh");
		cli.set_connection_timeout(5);
		cli.set_read_timeout(10);
		cli.Post("/", payload.dump(), "application/json");
	} catch (...) {
		// silently ignore network errors — MCMC must not be interrupted
	}
}

void TNtfyNotifier::notify_start(const std::vector<std::string> &tree_names,
                                 const std::vector<size_t> &leaf_counts, size_t n_iterations,
                                 size_t n_burnin, size_t n_burnin_iters) const {
	std::string msg;
	for (size_t i = 0; i < tree_names.size(); ++i) {
		msg += tree_names[i] + ": " + std::to_string(leaf_counts[i]) + " leaves\n";
	}
	msg += std::to_string(n_burnin) + " burnins x " + std::to_string(n_burnin_iters) + " iters\n";
	msg += "MCMC: " + std::to_string(n_iterations) + " iterations";
	_send("ACOL run started", msg, 3, "rocket");
}

void TNtfyNotifier::notify_burnin_round(size_t round, size_t total_rounds,
                                        const std::vector<std::string> &dim_names,
                                        const std::vector<ParamStats> &gamma_stats,
                                        const ParamStats &epsilon_stats) const {
	std::ostringstream msg;
	msg << std::setprecision(6);
	for (size_t i = 0; i < dim_names.size(); ++i) {
		msg << "gamma_" << dim_names[i] << ": mean=" << gamma_stats[i].mean
		    << " var=" << gamma_stats[i].var << " sd=" << gamma_stats[i].sd << "\n";
	}
	msg << "epsilon: mean=" << epsilon_stats.mean << " var=" << epsilon_stats.var
	    << " sd=" << epsilon_stats.sd;
	std::string title =
	    "Burnin " + std::to_string(round) + "/" + std::to_string(total_rounds) + " finished";
	_send(title, msg.str(), 2, "white_check_mark");
}

void TNtfyNotifier::notify_burnin_finished(const std::vector<std::string> &dim_names,
                                           const std::vector<ParamStats> &gamma_stats,
                                           const ParamStats &epsilon_stats) const {
	std::ostringstream msg;
	msg << std::setprecision(6);
	for (size_t i = 0; i < dim_names.size(); ++i) {
		msg << "gamma_" << dim_names[i] << ": mean=" << gamma_stats[i].mean
		    << " var=" << gamma_stats[i].var << " sd=" << gamma_stats[i].sd << "\n";
	}
	msg << "epsilon: mean=" << epsilon_stats.mean << " var=" << epsilon_stats.var
	    << " sd=" << epsilon_stats.sd;
	_send("Burnin finished", msg.str(), 2, "white_check_mark");
}

void TNtfyNotifier::notify_mcmc_finished(const std::vector<std::string> &dim_names,
                                         const std::vector<ParamStats> &gamma_stats,
                                         const ParamStats &epsilon_stats) const {
	std::ostringstream msg;
	msg << std::setprecision(6);
	for (size_t i = 0; i < dim_names.size(); ++i) {
		msg << "gamma_" << dim_names[i] << ": mean=" << gamma_stats[i].mean
		    << " var=" << gamma_stats[i].var << " sd=" << gamma_stats[i].sd << "\n";
	}
	msg << "epsilon: mean=" << epsilon_stats.mean << " var=" << epsilon_stats.var
	    << " sd=" << epsilon_stats.sd;
	_send("MCMC finished", msg.str(), 4, "tada");
}
