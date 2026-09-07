#ifndef SCC_PAPER_SWAP_H
#define SCC_PAPER_SWAP_H

#include <math.h>
#include "parlay/random.h"

#include <stdio.h>

#include <algorithm>
#include <iostream>
#include <queue>  // std::queue

#include "get_time.hpp"
#include "graph.hpp"
#include "hash_bag.h"
#include "reach.hpp"
#include "resizable_table.h"
#include "utilities.h"

using namespace std;

constexpr size_t TOP_BIT = size_t(UINT_N_MAX) + 1;
constexpr size_t VAL_MASK = UINT_N_MAX;
size_t _min_bag_size_;
// size_t multi_tau;
// size_t per_core;
using label_type = size_t;

using K = NodeId;
using V = NodeId;
using T = std::tuple<K, V>;
struct hash_kv {
  NodeId operator()(const K& k) { return _hash(k); }
};

// hash32 is sufficient

class SCC {
 private:
  Graph& graph;
  //------------------------------------------------------------
  gbbs::resizable_table<K, V, hash_kv> table_forw;
  gbbs::resizable_table<K, V, hash_kv> table_back;
  hashbag<NodeId> bag;
  sequence<NodeId> front;
  sequence<NodeId> substep_front;
  sequence<bool> bits;
  sequence<NodeId> centers;
  sequence<uint32_t> bag_stamp;
  sequence<EdgeId> out_scope;
  sequence<EdgeId> in_scope;
  sequence<EdgeId> out_scope_begin;
  sequence<EdgeId> in_scope_begin;
  sequence<EdgeId> out_scope_end;
  sequence<EdgeId> in_scope_end;
  sequence<uint32_t> out_scope_epoch;
  sequence<uint32_t> in_scope_epoch;
  sequence<uint8_t> pivot_valid;
  uint32_t current_stamp;
  uint32_t current_scope_epoch = 0;
  size_t n_front;
  size_t n;
  size_t label_offset;
  bool conservative_tables = false;
  timer bits_clear_timer;
  NodeId c =100;
  bool edge_map(NodeId cur_node, NodeId ngb_node,
                gbbs::resizable_table<K, V, hash_kv>& table);
  bool edge_map_first(NodeId cur_node, NodeId ngb_node,
                      gbbs::resizable_table<K, V, hash_kv>& table);
  void identify_fake_links(NodeId v, bool forward, EdgeId scan_begin,
                           EdgeId scan_end, EdgeId first_real,
                           EdgeId last_real, EdgeId real_count);
  bool is_potential_singleton(NodeId v, sequence<size_t>& labels);
  int multi_search(sequence<size_t>& labels,
                   gbbs::resizable_table<K, V, hash_kv>& table, bool forward,
                   bool local);
  int multi_search_safe(sequence<size_t>& labels,
                        gbbs::resizable_table<K, V, hash_kv>& table,
                        bool forward, bool local);
  void time_stamp_update() {
    current_stamp++;
    if (current_stamp == numeric_limits<uint32_t>::max() - 1) {
      parallel_for(0, n, [&](size_t i) {
        bag_stamp[i] =
            (bag_stamp[i] ==
             current_stamp);  // if current stamp will overflow, reset to 1
      });
      current_stamp = 1;
    }
  }

 public:
  size_t front_thresh;
  size_t num_round;
  void scc(sequence<size_t>& labels, double beta, bool local1, bool local2);
  SCC() = delete;
  SCC(Graph& G) : graph(G), bag(G.n, _min_bag_size_) {
    n = graph.n;
    num_round = 0;
    front = sequence<NodeId>(n);
    substep_front = sequence<NodeId>(n);
    bits = sequence<bool>(n);
    bag_stamp = sequence<uint32_t>(bag.get_bag_capacity());
    out_scope = sequence<EdgeId>(n);
    in_scope = sequence<EdgeId>(n);
    out_scope_begin = sequence<EdgeId>(n);
    in_scope_begin = sequence<EdgeId>(n);
    out_scope_end = sequence<EdgeId>(n);
    in_scope_end = sequence<EdgeId>(n);
    out_scope_epoch = sequence<uint32_t>(n);
    in_scope_epoch = sequence<uint32_t>(n);
    pivot_valid = sequence<uint8_t>(n);
    parallel_for(0, n, [&](size_t i) {
      out_scope_epoch[i] = 0;
      in_scope_epoch[i] = 0;
    });
    table_forw = gbbs::resizable_table<K, V, hash_kv>(
        n, {UINT_N_MAX, UINT_N_MAX}, hash_kv());
    table_back = gbbs::resizable_table<K, V, hash_kv>(
        n, {UINT_N_MAX, UINT_N_MAX}, hash_kv());
  }

  timer init_timer;
  timer first_round_timer;
  timer multi_search_timer;
  //timer trim_timer;
  timer multi_search_safe_timer;
  timer scc_timer;
  void timer_reset() {
    init_timer.reset();
    first_round_timer.reset();
    multi_search_timer.reset();
    //trim_timer.reset();
    multi_search_safe_timer.reset();
    scc_timer.reset();
  }
  void breakdown_report(int repeat) {
    cout << "average init_time " << init_timer.get_total() / repeat << endl;
    cout << "average first_round_time "
         << first_round_timer.get_total() / repeat << endl;
    cout << "average multi_search time "
         << multi_search_timer.get_total() / repeat << endl;
    //cout << "average trim time "
    //     << trim_timer.get_total() / repeat << endl;
    cout << "average hash table resize time "
        << (multi_search_safe_timer.get_total()-multi_search_timer.get_total())/repeat << endl;
  }
};

bool SCC::edge_map(NodeId cur_node, NodeId ngb_node,
                   gbbs::resizable_table<K, V, hash_kv>& table) {
  bool labels_changed = false;
  auto s_iter = table.get_iter(cur_node);
  if (s_iter.init()) {
    auto table_entry = s_iter.next();
    NodeId s_label = std::get<1>(table_entry);
    labels_changed |= table.insert(std::make_tuple(ngb_node, s_label));
    if (table.is_full()) return false;

    while (s_iter.has_next()) {
      auto table_entry = s_iter.next();
      s_label = std::get<1>(table_entry);
      labels_changed |= table.insert(std::make_tuple(ngb_node, s_label));
      if (table.is_full()) return false;
    }
  }

  return labels_changed;
}

bool SCC::edge_map_first(NodeId s_label, NodeId ngb_node,
                         gbbs::resizable_table<K, V, hash_kv>& table) {
  return table.insert(std::make_tuple(ngb_node, s_label));
}

// Records fake-link information discovered by the reachability scan that was
// already required by GBBS-VGC.  Edges outside [first_real,last_real] are
// permanently irrelevant because SCC subproblem labels only split and never
// merge.  Restricting later scans to this range recovers the paper's scope
// optimization without copying or concurrently mutating the input CSR.
void SCC::identify_fake_links(NodeId v, bool forward, EdgeId scan_begin,
                              EdgeId scan_end, EdgeId first_real,
                              EdgeId last_real, EdgeId real_count) {
  auto& scope = forward ? out_scope : in_scope;
  auto& scope_begin = forward ? out_scope_begin : in_scope_begin;
  auto& scope_end = forward ? out_scope_end : in_scope_end;
  auto& scope_epoch = forward ? out_scope_epoch : in_scope_epoch;

  const EdgeId new_begin =
      real_count == 0 ? scan_end : std::max(scan_begin, first_real);
  const EdgeId new_end = real_count == 0 ? scan_end : last_real + 1;
  __atomic_store_n(&scope[v], real_count, __ATOMIC_RELAXED);
  __atomic_store_n(&scope_begin[v], new_begin, __ATOMIC_RELAXED);
  __atomic_store_n(&scope_end[v], new_end, __ATOMIC_RELAXED);
  __atomic_store_n(&scope_epoch[v], current_scope_epoch, __ATOMIC_RELEASE);
}

bool SCC::is_potential_singleton(NodeId v, sequence<size_t>& labels) {
  if (__atomic_load_n(&labels[v], __ATOMIC_RELAXED) & TOP_BIT) return false;
  const bool no_out =
      __atomic_load_n(&out_scope_epoch[v], __ATOMIC_ACQUIRE) ==
          current_scope_epoch &&
                      __atomic_load_n(&out_scope[v], __ATOMIC_RELAXED) == 0;
  const bool no_in =
      __atomic_load_n(&in_scope_epoch[v], __ATOMIC_ACQUIRE) ==
          current_scope_epoch &&
                     __atomic_load_n(&in_scope[v], __ATOMIC_RELAXED) == 0;
  return no_out || no_in;
}

int SCC::multi_search(sequence<size_t>& labels,
                      gbbs::resizable_table<K, V, hash_kv>& table, bool forward,
                      bool local) {
  parallel_for(0, n_front, [&](size_t i) {
    table.insert(make_tuple(front[i], label_offset + i));
    substep_front[i] = front[i];
  });
  table.update_nelms();
  size_t sub_n_front = n_front;

  sequence<EdgeId>& offset = forward ? graph.offset : graph.in_offset;
  sequence<NodeId>& E = forward ? graph.E : graph.in_E;
  sequence<EdgeId>& scope = forward ? out_scope : in_scope;
  sequence<EdgeId>& scope_begin =
      forward ? out_scope_begin : in_scope_begin;
  sequence<EdgeId>& scope_end = forward ? out_scope_end : in_scope_end;
  sequence<uint32_t>& scope_epoch =
      forward ? out_scope_epoch : in_scope_epoch;
  size_t round = 0;
  size_t queue_size = tau;
  int sub_round = 0;
  while (sub_n_front > 0) {
    sub_round++;
    time_stamp_update();
    round++;
    #if defined(DEBUG)
    cout << "multi-search round " << round << " frontier " << sub_n_front
         << endl;
    cout << "resizable table m " << table.m << endl;
    #endif  
    parallel_for(
        0, sub_n_front, [&](size_t i) { bits[substep_front[i]] = 0; }, 2048);
    parallel_for(
        0, sub_n_front,
        [&](size_t i) {
          NodeId node = substep_front[i];
          const bool node_scope_known =
              __atomic_load_n(&scope_epoch[node], __ATOMIC_ACQUIRE) ==
              current_scope_epoch;
          const EdgeId node_begin = node_scope_known
                                        ? scope_begin[node]
                                        : offset[node];
          const EdgeId node_end = node_scope_known
                                      ? scope_end[node]
                                      : offset[node + 1];
          EdgeId degree = node_end - node_begin;
          EdgeId edge_map_cnt = 0;
          if (local && (degree < queue_size) &&
              (edge_map_cnt < queue_size)) {
            NodeId Q[queue_size];
            int head = 0, tail = 0;
            Q[tail++] = node;
            while (head < tail && edge_map_cnt < queue_size) {
              NodeId u = Q[head++];
              const bool scope_is_known =
                  __atomic_load_n(&scope_epoch[u], __ATOMIC_ACQUIRE) ==
                  current_scope_epoch;
              const EdgeId begin = scope_is_known
                                       ? scope_begin[u]
                                       : offset[u];
              const EdgeId end = scope_is_known
                                     ? scope_end[u]
                                     : offset[u + 1];
              const EdgeId deg = end - begin;
              if (deg >= queue_size) break;

              EdgeId first_real = end;
              EdgeId last_real = begin;
              EdgeId real_count = 0;
              const size_t source_label = labels[u];
              for (EdgeId j = begin; j < end; j++) {
                NodeId v = E[j];
                const size_t target_label = labels[v];
                if (!(target_label & TOP_BIT) &&
                    target_label == source_label) {
                  if (real_count == 0) first_real = j;
                  last_real = j;
                  ++real_count;
                  edge_map_cnt++;
                  bool edge_map_success =
                      (sub_round == 1 && local)
                          ? edge_map_first(label_offset + i, v, table)
                          : edge_map(node, v, table);
                  if (edge_map_success) {
                    if (edge_map_cnt < queue_size) {
                      Q[tail++] = v;
                    } else if (compare_and_swap(&bits[v], false, true)) {
                      bag.insert(v, make_slice(bag_stamp), current_stamp);
                    }
                  }
                }
              }
              if (real_count != deg || deg == 0) {
                identify_fake_links(u, forward, begin, end, first_real,
                                    last_real, real_count);
              }
            }

            if (head < tail) {
              for (int j = head; j < tail; j++) {
                NodeId u = Q[j];
                if (compare_and_swap(&bits[u], false, true)) {
                  bag.insert(u, make_slice(bag_stamp), current_stamp);
                }
              }
            }
          } else if (degree > 0) {
            parallel_for(
                0, degree,
                [&](size_t j) {
                  NodeId ngb_node = E[node_begin + j];
                    if (!(labels[ngb_node] & TOP_BIT) &&
                      (labels[ngb_node] == labels[node]) ) {
                    bool edge_map_success =
                        (sub_round == 1 && local)
                            ? edge_map_first(label_offset + i, ngb_node, table)
                            : edge_map(node, ngb_node, table);
                    if (edge_map_success) {
                      if (compare_and_swap(&bits[ngb_node], false, true)) {
                        bag.insert(ngb_node, make_slice(bag_stamp),
                                   current_stamp);
                      }
                    }
                  }
                },
                1024);
          }
        },
        1);
    if (local && sub_round == 1) {
      // Propagation above is read-only.  Once it completes, compact only
      // very-low-degree roots with at least two-thirds fake links.  This covers
      // the paper's motivating pattern without rescanning medium-degree lists.
      parallel_for(
          0, sub_n_front,
          [&](size_t i) {
            const NodeId node = substep_front[i];
            const EdgeId begin = offset[node];
            const EdgeId original_end = offset[node + 1];
            const EdgeId degree = original_end - begin;
            if (degree == 0 || degree > 4) return;
            if (__atomic_load_n(&scope_epoch[node], __ATOMIC_ACQUIRE) !=
                current_scope_epoch) {
              return;
            }
            const EdgeId real_count =
                __atomic_load_n(&scope[node], __ATOMIC_RELAXED);
            if (real_count * 3 > degree) return;

            const size_t source_label = labels[node];
            EdgeId end = original_end;
            EdgeId j = begin;
            while (j < end) {
              const NodeId v = E[j];
              const size_t target_label = labels[v];
              if ((target_label & TOP_BIT) ||
                  target_label != source_label) {
                --end;
                if (j != end) std::swap(E[j], E[end]);
              } else {
                ++j;
              }
            }
            identify_fake_links(
                node, forward, begin, end, begin,
                end == begin ? begin : end - 1, end - begin);
          },
          2048);
    }
    if (table.is_full()) {
      return -1;
    }
    sub_n_front = bag.pack(make_slice(substep_front), make_slice(bag_stamp),
                           current_stamp);
    table.update_nelms();
#if defined(DEBUG)
    cout << "table.ne = " << table.ne << endl;
#endif
  }
  return round;
}

/*
int SCC::multi_search(sequence<size_t>& labels,
                      gbbs::resizable_table<K, V, hash_kv>& table, bool forward,
                      bool local) {
  parallel_for(0, n_front, [&](size_t i) {
    table.insert(make_tuple(front[i], label_offset + i));
    substep_front[i] = front[i];
  });
  table.update_nelms();
  size_t sub_n_front = n_front;

  //sequence<EdgeId>& offset = forward ? graph.offset : graph.in_offset;
  sequence<EdgeId>& offset = forward ? valid_offset : in_valid_offset;
  sequence<NodeId>& E = forward ? graph.E : graph.in_E;
  sequence<NodeId>& valid_count = forward ? out_valid_count : in_valid_count;
 


  size_t round = 0;
  size_t queue_size = tau;
  int sub_round = 0;
  while (sub_n_front > 0) {
    sub_round++;
    time_stamp_update();
    round++;
    #if defined(DEBUG)
    cout << "multi-search round " << round << " frontier " << sub_n_front
         << endl;
    cout << "resizable table m " << table.m << endl;
    #endif  
    parallel_for(
        0, sub_n_front, [&](size_t i) { bits[substep_front[i]] = 0; }, 2048);
    parallel_for(
        0, sub_n_front,
        [&](size_t i) {
          NodeId node = substep_front[i];
          //EdgeId degree = offset[node + 1] - offset[node];
          NodeId valid_deg = valid_count[node];
          EdgeId edge_map_cnt = 0;
          if (local && valid_deg && (valid_deg < queue_size) && (edge_map_cnt < queue_size)) {
            NodeId Q[queue_size];
            int head = 0, tail = 0;
            Q[tail] = node;
            tail++;
            while (head < tail && (edge_map_cnt < queue_size)) {
              NodeId u = Q[head];
              //EdgeId deg = offset[u + 1] - offset[u];
              EdgeId deg = valid_count[u];
              //EdgeId deg_n = valid_count_n[u];
              if (deg >= queue_size ) {
                break;
              }
              head++;
              if(deg == 0 ){
                continue;
              }
              NodeId count = 0;
              NodeId u_end = 0;
              NodeId u_fir = 0;
              EdgeId u_off = offset[u];
              for (EdgeId j = 0; j <  deg; j++) {
                NodeId v = E[j + u_off];
                if (!(labels[v] & TOP_BIT) && (labels[v] == labels[u])) {
                  edge_map_cnt++;
                  count++;
                  if(count == 1){ u_fir = j;}
                  u_end = j;
                  bool edge_map_success =
                      (sub_round == 1 && local)
                          ? edge_map_first(label_offset + i, v, table)
                          : edge_map(node, v, table);
                  if (edge_map_success) {
                    if (edge_map_cnt < queue_size) {
                      Q[tail] = v;
                      tail++;
                    } else {
                      if (compare_and_swap(&bits[v], false, true)) {
                        bag.insert(v, make_slice(bag_stamp), current_stamp);
                      }
                    }
                  }
                  //if(count == deg){break;}//has visited all valid edges
                }
              }
              //if(count < deg){valid_count[u] = count;} //update valid 
              offset[u] = u_fir + u_off ;
              valid_count[u] = count!=0 ? u_end - u_fir + 1 : 0;
             

            }

            if (head < tail) {
              for (int j = head; j < tail; j++) {
                NodeId u = Q[j];
                if (compare_and_swap(&bits[u], false, true)) {
                  bag.insert(u, make_slice(bag_stamp), current_stamp);
                }
              }
            }

          } else if (valid_deg ) {
           
            NodeId count =0;
            //EdgeId first = 0 ;
            //EdgeId node_off = offset[node];
            EdgeId degree = offset[node + 1] - offset[node];
            bool count_f = false;
            parallel_for(
                0, degree,
                [&](size_t j) {
                  NodeId ngb_node = E[offset[node] + j];
                  if (!(labels[ngb_node] & TOP_BIT) &&
                      (labels[ngb_node] == labels[node]) ) {
                    //flag_n[j] = true;
                    
                    if(compare_and_swap(&count_f, false, true)){
                      count++;
                      //if(count == 1){first =  off;}
                    }
                    
                    bool edge_map_success =
                        (sub_round == 1 && local)
                            ? edge_map_first(label_offset + i, ngb_node, table)
                            : edge_map(node, ngb_node, table);
                    if (edge_map_success) {
                      if (compare_and_swap(&bits[ngb_node], false, true)) {
                        bag.insert(ngb_node, make_slice(bag_stamp),
                                   current_stamp);
                      }
                    }
                  }
                },
                1024);
          
          if(count < valid_count[node]){valid_count[node] = count;}
          //offset[node] = first ;
          //if(count == 0){valid_count[node] = 0;}
          
          
          }
        },
        1);
    if (table.is_full()) {
      return -1;
    }
    sub_n_front = bag.pack(make_slice(substep_front), make_slice(bag_stamp),
                           current_stamp);
    table.update_nelms();
#if defined(DEBUG)
    cout << "table.ne = " << table.ne << endl;
#endif
  }
  return round;
}
*/

/*
int SCC::multi_search(sequence<size_t>& labels,
                      gbbs::resizable_table<K, V, hash_kv>& table, bool forward,
                      bool local) {
  parallel_for(0, n_front, [&](size_t i) {
    table.insert(make_tuple(front[i], label_offset + i));
    substep_front[i] = front[i];
  });
  table.update_nelms();
  size_t sub_n_front = n_front;

  sequence<EdgeId>& offset = forward ? graph.offset : graph.in_offset;
  sequence<NodeId>& E = forward ? graph.E : graph.in_E;
  sequence<NodeId>& cut_off = forward ? cut_offset : in_cut_offset;
  sequence<NodeId>& valid_count = forward ? out_valid_count : in_valid_count;
  sequence<NodeId>& valid_count_n = forward ? in_valid_count : out_valid_count;


  size_t round = 0;
  size_t queue_size = tau;
  int sub_round = 0;
  while (sub_n_front > 0) {
    sub_round++;
    time_stamp_update();
    round++;
    #if defined(DEBUG)
    cout << "multi-search round " << round << " frontier " << sub_n_front
         << endl;
    cout << "resizable table m " << table.m << endl;
    #endif  
    parallel_for(
        0, sub_n_front, [&](size_t i) { bits[substep_front[i]] = 0; }, 2048);
    parallel_for(
        0, sub_n_front,
        [&](size_t i) {
          NodeId node = substep_front[i];
          NodeId valid_deg = valid_count[node];
          EdgeId edge_map_cnt = 0;
          if (local && valid_deg && valid_count_n[node]&& (valid_deg < queue_size) && (edge_map_cnt < queue_size)) {
            NodeId Q[queue_size];
            int head = 0, tail = 0;
            Q[tail] = node;
            tail++;
            while (head < tail && (edge_map_cnt < queue_size)) {
              NodeId u = Q[head];
              EdgeId deg = valid_count[u];
              EdgeId deg_n = valid_count_n[u];
              if (deg >= queue_size ) {
                break;
              }
              head++;
              if(deg == 0 ||deg_n == 0){
                continue;
              }
              NodeId count = 0;
              NodeId cut_count =0;
              EdgeId base_offset = offset[u];
              for (EdgeId j = 0; j < cut_off[u]; j++) {
                NodeId v = E[base_offset + j];
                if (!(labels[v] & TOP_BIT) && (labels[v] == labels[u])) {
                  edge_map_cnt++;
                  //if(valid_count[v] == 0){continue;}
                  count++;
                  cut_count =  j;
                  bool edge_map_success =
                      (sub_round == 1 && local)
                          ? edge_map_first(label_offset + i, v, table)
                          : edge_map(node, v, table);
                  if (edge_map_success) {
                    if (edge_map_cnt < queue_size) {
                      Q[tail] = v;
                      tail++;
                    } else {
                      if (compare_and_swap(&bits[v], false, true)) {
                        bag.insert(v, make_slice(bag_stamp), current_stamp);
                      }
                    }
                  }
                  if(count == deg){break;}//has visited all valid edges
                }
              }
              if(count < deg){valid_count[u] = count;} //update valid count
              if(cut_off[u] > cut_count && count){ cut_off[u] = cut_count + 1;}
              

            }

            if (head < tail) {
              for (int j = head; j < tail; j++) {
                NodeId u = Q[j];
                if (compare_and_swap(&bits[u], false, true)) {
                  bag.insert(u, make_slice(bag_stamp), current_stamp);
                }
              }
            }

          } else if (valid_deg && valid_count_n[node]) {
            
           
            NodeId count =0;
            bool count_f = false;
            NodeId cut_count = 0;
            parallel_for(
                0, cut_off[node],
                [&](size_t j) {
                  NodeId ngb_node = E[offset[node] + j];
                  if (!(labels[ngb_node] & TOP_BIT) &&
                      (labels[ngb_node] == labels[node]) ) {
                    
                    if(compare_and_swap(&count_f, false, true)){
                      count++;
                      cut_count = j;
                    }
                    
                    bool edge_map_success =
                        (sub_round == 1 && local)
                            ? edge_map_first(label_offset + i, ngb_node, table)
                            : edge_map(node, ngb_node, table);
                    if (edge_map_success) {
                      if (compare_and_swap(&bits[ngb_node], false, true)) {
                        bag.insert(ngb_node, make_slice(bag_stamp),
                                   current_stamp);
                      }
                    }
                  }
                },
                1024);
         
          if(count < valid_deg){valid_count[node] = count;}
          if(cut_count > cut_off[node] && count){cut_off[node] = cut_count+1;}
          
          }
        },
        1);
    if (table.is_full()) {
      return -1;
    }
    sub_n_front = bag.pack(make_slice(substep_front), make_slice(bag_stamp),
                           current_stamp);
    table.update_nelms();
#if defined(DEBUG)
    cout << "table.ne = " << table.ne << endl;
#endif
  }
  return round;
}
*/



int SCC::multi_search_safe(sequence<size_t>& labels,
                           gbbs::resizable_table<K, V, hash_kv>& table,
                           bool forward, bool local) {
  int round;
  multi_search_safe_timer.start();
  double multi_search_time=0;

  while (true){
    multi_search_timer.start();
    round = multi_search(labels, table, forward, local);
    multi_search_time=multi_search_timer.stop();
    if (round != -1){break;}
    cout << "trigger table resize" << endl;
    multi_search_timer.total_time -= multi_search_time;
    conservative_tables = true;
    parallel_for(0, graph.n, [&](size_t i) { bits[i] = 0; });
    table.double_size();
  }
  multi_search_safe_timer.stop();
  return round;
}

void SCC::scc(sequence<size_t>& labels, double beta, bool local_reach,
              bool local_scc) {
  init_timer.start();
  scc_timer.reset();
  scc_timer.start();
  current_stamp = 1;
  ++current_scope_epoch;
  if (current_scope_epoch == 0) {
    parallel_for(0, n, [&](size_t i) {
      out_scope_epoch[i] = 0;
      in_scope_epoch[i] = 0;
    });
    current_scope_epoch = 1;
  }

  sequence<bool> dist_1(n);
  sequence<bool> dist_2(n);

  sequence<NodeId> vertices = sequence<NodeId>(n);
  parallel_for(0, n, [&](size_t i) {
    labels[i] = 0;
    bits[i] = false;
    vertices[i] = i;
  });
  parallel_for(0, bag.get_bag_capacity(), [&](size_t i) { bag_stamp[i] = 0; });

  auto NON_ZEROS = parlay::filter(vertices, [&](NodeId v) {
    return ((graph.offset[v + 1] - graph.offset[v]) != 0) &&
           ((graph.in_offset[v + 1] - graph.in_offset[v]) != 0);
  });
  auto ZEROS = parlay::filter(vertices, [&](NodeId v) {
    return ((graph.offset[v + 1] - graph.offset[v]) == 0) ||
           ((graph.in_offset[v + 1] - graph.in_offset[v]) == 0);
  });
  auto P = parlay::random_shuffle(NON_ZEROS);
  parallel_for(0, ZEROS.size(),
               [&](NodeId i) { labels[ZEROS[i]] = 1 + (i | TOP_BIT); });
  init_timer.stop();

  // cout << "------------------------------------" << endl;
  // cout << "initial time " << scc_timer.stop() << endl;
  // cout << "Filtered: " << ZEROS.size()
  //      << " vertices. Num remaining = " << NON_ZEROS.size() << endl;
  scc_timer.start();
  NodeId step_size = 1, cur_offset = 0, cur_round = 0;
  label_offset = ZEROS.size() + 1;

  if (P.empty()) {
    parallel_for(0, graph.n,
                 [&](size_t i) { labels[i] = (labels[i] & VAL_MASK) - 1; });
    scc_timer.stop();
    return;
  }

  // ----------------first round, for the BIG SCC -------------------------
  NodeId source = P[0];
  
  // BFS BFS_P(graph);
  first_round_timer.start();
  REACH REACH_P(graph);
  // forward search
  #ifdef ROUND
  int fowd_depth = REACH_P.reach(source, dist_1, local_reach);
  #else
  REACH_P.reach(source, dist_1, local_reach);
  #endif
  // backward search
  REACH_P.swap_graph();
  #ifdef ROUND
  int bawd_depth = REACH_P.reach(source, dist_2, local_reach);
  #else
  REACH_P.reach(source, dist_2, local_reach);
  #endif
  #ifdef ROUND
  printf("Single Reach forward search depth %d\n", fowd_depth);
  printf("Single Reach backward search depth %d\n", bawd_depth);
  #endif

  REACH_P.swap_graph();
  NodeId label = label_offset;
  parallel_for(0, P.size(), [&](size_t i) {
    if ((dist_1[P[i]] != false) && (dist_2[P[i]] != false)) {
      labels[P[i]] = label | TOP_BIT;
    } else if ((dist_1[P[i]] != false) || (dist_2[P[i]] != false)) {
      labels[P[i]] = label;
      //node_bit[P[i]] = 0x04;//from set
    }
    /*
    else if((dist_1[P[i]] != false) || (dist_2[P[i]] == false)){
      labels[P[i]] = label;
      //node_bit[P[i]] = 0x02;//to set
    }
    */

  });

  centers = parlay::filter(P, [&](NodeId v) { return !(labels[v] & TOP_BIT); });
 
  first_round_timer.stop();
  #if defined(DEBUG)
  // cout << "First Round Time " << scc_timer.stop() << endl;
  // scc_timer.start();
  #endif
  label_offset++;
  // cout << "After first round, |FRONT| = " << centers.size()
  //      << " vertices remain. Total done = " << P.size() - centers.size()
  //      << endl;
  // cout << "centers: ";
  // for (int i = 0; i<10; i++){
  //     cout << centers[i] << " ";
  // }
  // cout << endl;

  timer t_search;
  timer t_round;
  NodeId n_remaining = centers.size();

  NodeId in_table_ne = 1;
  NodeId out_table_ne = 1;
  while (cur_offset < n_remaining) {
    t_round.start();
    const NodeId target_size =
        std::min(step_size, n_remaining - cur_offset);
    const NodeId refill_budget = std::min<NodeId>(
        4096, std::max<NodeId>(1, target_size / 64));
    const NodeId scan_limit =
        std::min(n_remaining, cur_offset + target_size + refill_budget);
    n_front = 0;

    // Paper pivot selection: resolved vertices and dynamically exposed
    // singleton SCCs do not consume a pivot slot.  Continue scanning the
    // permutation until the requested batch is full or no candidates remain.
    while (n_front < target_size && cur_offset < scan_limit) {
      const NodeId needed = target_size - n_front;
      const NodeId end = std::min(cur_offset + needed, scan_limit);
      auto candidates = centers.cut(cur_offset, end);

      parallel_for(0, candidates.size(), [&](size_t i) {
        const NodeId v = candidates[i];
        if (__atomic_load_n(&labels[v], __ATOMIC_RELAXED) & TOP_BIT) {
          pivot_valid[v] = 0;
        } else if (is_potential_singleton(v, labels)) {
          const size_t singleton_label =
              __atomic_fetch_add(&label_offset, size_t{1}, __ATOMIC_RELAXED);
          __atomic_store_n(&labels[v], singleton_label | TOP_BIT,
                           __ATOMIC_RELAXED);
          pivot_valid[v] = 0;
        } else {
          pivot_valid[v] = 1;
        }
      });

      auto selected = parlay::filter(
          candidates, [&](NodeId v) { return pivot_valid[v] != 0; });
      parallel_for(0, selected.size(),
                   [&](size_t i) { front[n_front + i] = selected[i]; });
      n_front += selected.size();
      cur_offset = end;
    }

    step_size = ceil(step_size * beta);
    cur_round++;
    // cout << "n_front = " << n_front << ", origin was " << vs_size
    //      << " Total vertices remaining = " << n_remaining - finished << endl;
    if (n_front == 0) continue;
    // cout << "QUEUE_SIZE " << tau << endl;
    t_search.start();
    T empty = std::make_tuple(UINT_N_MAX, UINT_N_MAX);
    // size_t out_table_m = max((size_t)min((NodeId)ceil(0.35*n_remaining) ,
    // (NodeId)6000000), (size_t)(beta)*out_table_ne);
    size_t out_table_base =
        max((size_t)min((NodeId)ceil(0.3 * n_remaining), (NodeId)6000000),
            (size_t)(beta)*out_table_ne);
    size_t out_table_m = conservative_tables
                             ? 2 * out_table_base
                             : (5 * out_table_base) / 4;
    table_forw =
        gbbs::resizable_table<K, V, hash_kv>(out_table_m, empty, hash_kv());
    #ifdef ROUND
    int out_depth = multi_search_safe(labels, table_forw, true, local_scc);
    #else
    multi_search_safe(labels, table_forw, true, local_scc);
    #endif
    t_search.stop();
    t_search.start();
    // size_t in_table_m = max((size_t)min((NodeId)ceil(0.35*n_remaining),
    // (NodeId)6000000), (size_t)(beta)*in_table_ne);
    size_t in_table_base =
        max((size_t)min((NodeId)ceil(0.3 * n_remaining), (NodeId)6000000),
            (size_t)(beta)*in_table_ne);
    size_t in_table_m = conservative_tables
                            ? 2 * in_table_base
                            : (5 * in_table_base) / 4;
    table_back =
        gbbs::resizable_table<K, V, hash_kv>(in_table_m, empty, hash_kv());
    #ifdef ROUND
    int in_depth = multi_search_safe(labels, table_back, false, local_scc);
    t_search.stop();
    #else
    multi_search_safe(labels, table_back, false, local_scc);
    t_search.stop();
    #endif

    label_offset += n_front;
    // std::cout << "in_table, m = " << table_back.m << " ne = " << table_back.ne
    //           << "\n";
    // std::cout << "out_table, m = " << table_forw.m << " ne = " << table_forw.ne
    //           << "\n";
    // cout << "In search time " << backward_time << endl;
    // cout << "Out search time " << forward_time << endl;
    #ifdef ROUND
    printf("round %d forward search depth %d\n", cur_round, out_depth);
    printf("round %d backward search depth %d\n", cur_round, in_depth);
    #endif

    auto& smaller_t = (table_forw.m <= table_back.m) ? table_forw : table_back;
    auto& larger_t = (table_forw.m > table_back.m) ? table_forw : table_back;

    // intersect the tables
    //step 1：get the node set to parallel
    //step 2: get the small_t size of node i
    //step 3: get the large_t size of node i
    //step 4: make intersecton.


    auto map_f = [&](const std::tuple<K, V>& kev) {
      NodeId v = std::get<0>(kev);
      size_t label = std::get<1>(kev);
      if (larger_t.contains(v, label)) {
        // in 'label' scc
        // Max visitor from this StronglyConnectedComponents acquires it.
        
        //if(in_valid_count[v]!=0 && out_valid_count[v]!=0){
          write_max(&labels[v], label | TOP_BIT,
                  [&](size_t a, size_t b) { return a < b; });
          
        //}
        //else{

          //        labels[v] = label | TOP_BIT;

        //}
      } else {
        write_max(&labels[v], label, [&](size_t a, size_t b) { return a < b; });
        //node_bit[v] = 0x02;
      }
    };
    smaller_t.map(map_f);

    // set the subproblems
    auto sp_map = [&](const std::tuple<K, V>& kev) {
      NodeId v = std::get<0>(kev);
      size_t label = std::get<1>(kev);
      // note that if v is already in an StronglyConnectedComponents (from (1)),
      // the pbbslib::write_max will read, compare and fail, as the top bit is
      // already set.
      write_max(&labels[v], label, [&](size_t a, size_t b) { return a < b; });
      //node_bit[v] = 0x04;
    };
    larger_t.map(sp_map);
  /*
    if(cur_round   == round_set ){
    trim_timer.start();
    auto trim = parlay::filter(centers.cut(first_vis,n_remaining), [&](NodeId v) { return (!(labels[v] & TOP_BIT))&&
                                                       ((out_valid_count[v] == 0)||(in_valid_count[v] == 0)); });
          
          parallel_for(0, trim.size(),
                  [&](NodeId i) { 
                    labels[trim[i]] = i+(label_offset| TOP_BIT); 
                  });
          label_offset += trim.size();
          cout<<trim.size()<<endl;

        //edge_set_after_trim(trim);
      
     trim_timer.stop();

    }
    */
  
    
    

    

    out_table_ne = table_forw.ne;
    in_table_ne = table_back.ne;

    // in_table.del();
    // out_table.del();
    t_round.stop();
    // cout << "Round time " << round_time << endl;
  }
  /*
  trim_timer.start();
   auto trim = parlay::filter(centers, [&](NodeId v) { return (!(labels[v] & TOP_BIT)); });
          
          parallel_for(0, trim.size(),
                  [&](NodeId i) { 
                    labels[trim[i]] = i+(label_offset| TOP_BIT); 
                  });
          label_offset += trim.size();
          cout<<trim.size()<<endl;

        //edge_set_after_trim(trim);
      
     trim_timer.stop();
  */

  parallel_for(0, graph.n,
               [&](NodeId i) { labels[i] = (labels[i] & VAL_MASK) - 1; });
  scc_timer.stop();
  // cout << "scc_time " << scc_timer.get_total() << endl;
}

#endif
