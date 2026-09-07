#include "scc_paper_swap_ablation.hpp"

#include <algorithm>
#include <iomanip>
#include <limits>
#include <string>

#include "get_time.hpp"
#include "parseCommandLine.hpp"

using namespace std;

struct MethodConfig {
  string name;
  bool fake_links;
  bool bounded_pivot;
  bool gbbs_vgc_compatibility;
};

struct PartitionStats {
  size_t scc_count;
  NodeId largest_scc;
};

bool parse_method(const string& name, MethodConfig& config) {
  if (name == "GBBS-VGC") {
    config = {name, false, false, true};
  } else if (name == "Fast-SCC-FL") {
    config = {name, true, false, false};
  } else if (name == "Fast-SCC-PS") {
    config = {name, false, true, false};
  } else if (name == "Fast-SCC") {
    config = {name, true, true, false};
  } else {
    return false;
  }
  return true;
}

PartitionStats partition_stats(const sequence<size_t>& labels,
                               sequence<NodeId>& scratch) {
  parallel_for(0, scratch.size(), [&](size_t i) { scratch[i] = 0; });
  for (size_t i = 0; i < labels.size(); ++i) {
    if (labels[i] >= scratch.size()) return {0, 0};
    ++scratch[labels[i]];
  }
  size_t count = 0;
  NodeId largest = 0;
  for (size_t i = 0; i < scratch.size(); ++i) {
    if (scratch[i] != 0) {
      ++count;
      largest = std::max(largest, scratch[i]);
    }
  }
  return {count, largest};
}

bool equivalent_partitions(const sequence<size_t>& reference,
                           const sequence<size_t>& actual,
                           size_t reference_scc_count,
                           size_t actual_scc_count,
                           sequence<NodeId>& scratch) {
  if (reference.size() != actual.size() ||
      reference_scc_count != actual_scc_count) {
    return false;
  }
  parallel_for(0, scratch.size(),
               [&](size_t i) { scratch[i] = UINT_N_MAX; });
  for (size_t i = 0; i < reference.size(); ++i) {
    if (reference[i] >= scratch.size() || actual[i] >= scratch.size()) {
      return false;
    }
    NodeId& mapped = scratch[reference[i]];
    const NodeId actual_label = static_cast<NodeId>(actual[i]);
    if (mapped == UINT_N_MAX) {
      mapped = actual_label;
    } else if (mapped != actual_label) {
      return false;
    }
  }
  return true;
}

void print_result(const MethodConfig& config, int run, double seconds,
                  const AblationMetrics& m, const PartitionStats& stats,
                  bool correct) {
  const uint64_t forward_rounds =
      m.initial_forward_rounds + m.batch_forward_rounds;
  const uint64_t backward_rounds =
      m.initial_backward_rounds + m.batch_backward_rounds;
  const uint64_t classified_edge_checks =
      m.real_link_checks + m.fake_link_checks;
  const double fake_ratio = classified_edge_checks == 0
                                ? 0.0
                                : static_cast<double>(m.fake_link_checks) /
                                      classified_edge_checks;

  cout << fixed << setprecision(9)
       << "ABLATION_RESULT"
       << " revision=2"
       << " method=" << config.name << " run=" << run
       << " fake_link_identification=" << (config.fake_links ? "ON" : "OFF")
       << " bounded_pivot_selection=" << (config.bounded_pivot ? "ON" : "OFF")
       << " gbbs_vgc_compatibility="
       << (config.gbbs_vgc_compatibility ? "ON" : "OFF")
       << " runtime_seconds=" << seconds
       << " edge_checks_total=" << m.batch_edge_checks
       << " real_link_checks=" << m.real_link_checks
       << " fake_link_checks=" << m.fake_link_checks
       << " fake_link_ratio=" << fake_ratio
       << " partition_edge_checks=" << m.partition_edge_checks
       << " pivot_validation_edge_checks="
       << m.pivot_validation_edge_checks
       << " scope_reused_scans=" << m.scope_reused_scans
       << " scope_updates=" << m.scope_updates
       << " partitioned_lists=" << m.partitioned_lists
       << " swapped_entries=" << m.swapped_entries
       << " global_search_rounds=" << (forward_rounds + backward_rounds)
       << " global_forward_rounds=" << forward_rounds
       << " global_backward_rounds=" << backward_rounds
       << " initial_forward_rounds=" << m.initial_forward_rounds
       << " initial_backward_rounds=" << m.initial_backward_rounds
       << " batch_forward_rounds=" << m.batch_forward_rounds
       << " batch_backward_rounds=" << m.batch_backward_rounds
       << " executed_batches=" << m.executed_batches
       << " pivot_candidates_examined=" << m.pivot_candidates_examined
       << " pivots_selected=" << m.pivots_selected
       << " resolved_candidates_skipped=" << m.resolved_candidates_skipped
       << " singleton_candidates_finalized="
       << m.singleton_candidates_finalized
       << " refill_candidates_examined=" << m.refill_candidates_examined
       << " table_resizes=" << m.table_resizes
       << " baseline_invalid_candidates="
       << m.baseline_invalid_candidates
       << " scc_count=" << stats.scc_count
       << " largest_scc=" << stats.largest_scc
       << " correctness=" << (correct ? "PASS" : "FAIL") << endl;
}

int main(int argc, char** argv) {
  if (argc < 3) {
    cerr << "usage: ./scc_paper_swap_ablation METHOD GRAPH [-t N] "
            "[-local_reach] [-local_scc] [-large] [-beta B] [-tau N]"
         << endl;
    cerr << "METHOD: GBBS-VGC | Fast-SCC-FL | Fast-SCC-PS | Fast-SCC"
         << endl;
    return 1;
  }

  MethodConfig config;
  if (!parse_method(argv[1], config)) {
    cerr << "invalid METHOD: " << argv[1] << endl;
    return 2;
  }

  CommandLine P(argc, argv);
  char* file_name = argv[2];
  const double beta = P.getOptionDouble("-beta", 1.5);
  Graph graph = P.getOption("-large") ? read_large_graph(file_name)
                                      : read_graph(file_name);
  const bool local_reach = P.getOption("-local_reach");
  const bool local_scc = P.getOption("-local_scc");
  _min_bag_size_ = P.getOptionInt("-lambda", 1 << 10);
  const int default_tau = config.gbbs_vgc_compatibility
                              ? 512
                              : (graph.m >= 200000000 ? 2048 : 512);
  tau = P.getOptionInt("-tau", default_tau);
  const int repeat = P.getOptionInt("-t", 10);
  if (repeat <= 0) {
    cerr << "-t must be greater than zero" << endl;
    return 2;
  }

  cout << "ABLATION_CONFIG revision=2"
       << " method=" << config.name
       << " fake_link_identification=" << (config.fake_links ? "ON" : "OFF")
       << " bounded_pivot_selection=" << (config.bounded_pivot ? "ON" : "OFF")
       << " gbbs_vgc_compatibility="
       << (config.gbbs_vgc_compatibility ? "ON" : "OFF")
       << " graph=" << file_name << " vertices=" << graph.n
       << " edges=" << graph.m << " repeats=" << repeat << " beta=" << beta
       << " tau=" << tau << " local_reach=" << local_reach
       << " local_scc=" << local_scc << endl;
  cout << "ABLATION_NOTE edge_checks_total covers batched multi-search, "
          "in-place partition checks, and GBBS-VGC pivot-validation checks; "
          "the common initial single-pivot reach is reported in round metrics "
          "but not in edge-check metrics."
       << endl;

  sequence<size_t> labels(graph.n);
  sequence<size_t> reference_labels(graph.n);
  sequence<NodeId> scratch(graph.n);
  SCC solver(graph);

  // Generate one untimed OFF/OFF reference partition for every measured run.
  solver.set_ablation(false, false, config.gbbs_vgc_compatibility);
  solver.scc(reference_labels, beta, local_reach, local_scc);
  const PartitionStats reference_stats =
      partition_stats(reference_labels, scratch);
  cout << "ABLATION_REFERENCE method=GBBS-VGC scc_count="
       << reference_stats.scc_count
       << " largest_scc=" << reference_stats.largest_scc << endl;

  solver.set_ablation(config.fake_links, config.bounded_pivot,
                      config.gbbs_vgc_compatibility);
  solver.scc(labels, beta, local_reach, local_scc);
  solver.timer_reset();

  timer measured;
  bool all_correct = true;
  for (int i = 0; i < repeat; ++i) {
    measured.start();
    solver.scc(labels, beta, local_reach, local_scc);
    const double seconds = measured.stop();

    const AblationMetrics metrics = solver.metrics();
    const PartitionStats stats = partition_stats(labels, scratch);
    const bool correct = equivalent_partitions(
        reference_labels, labels, reference_stats.scc_count, stats.scc_count,
        scratch);
    print_result(config, i + 1, seconds, metrics, stats, correct);
    all_correct &= correct;
  }

  cout << "ABLATION_SUMMARY method=" << config.name
       << " repeats=" << repeat
       << " average_runtime_seconds=" << fixed << setprecision(9)
       << measured.get_total() / repeat
       << " correctness=" << (all_correct ? "PASS" : "FAIL") << endl;
  solver.breakdown_report(repeat);
  return all_correct ? 0 : 3;
}
