// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#include "tahpp/timing_criticality.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <unordered_set>

#include "db_sta/dbSta.hh"
#include "odb/db.h"
#include "sta/MinMax.hh"
#include "sta/PathEnd.hh"
#include "sta/PathExpanded.hh"
#include "sta/Search.hh"
#include "utl/Logger.h"

namespace mpl {
namespace tahpp {

namespace {

bool isMacroInst(odb::dbInst* inst)
{
  return inst && (inst->isBlock() || inst->isCore()) && !inst->isFixed();
}

double macroPairDistance(odb::dbInst* a, odb::dbInst* b)
{
  const odb::Rect ba = a->getBBox()->getBox();
  const odb::Rect bb = b->getBBox()->getBox();
  const int ax = ba.xCenter();
  const int ay = ba.yCenter();
  const int bx = bb.xCenter();
  const int by = bb.yCenter();
  return std::hypot(static_cast<double>(ax - bx), static_cast<double>(ay - by));
}

}  // namespace

TimingCriticalityMap TimingCriticalityMap::build(odb::dbBlock* block,
                                                 sta::dbSta* sta,
                                                 utl::Logger* logger,
                                                 const int nworst_paths)
{
  TimingCriticalityMap map;

  for (odb::dbInst* inst : block->getInsts()) {
    if (!isMacroInst(inst)) {
      continue;
    }
    map.macro_index[inst] = static_cast<int>(map.macros.size());
    map.macros.push_back(inst);
    map.path_criticality[inst] = 0.0;
  }

  if (map.macros.empty()) {
    return map;
  }

  if (!sta) {
    for (odb::dbInst* inst : map.macros) {
      map.path_criticality[inst] = 1.0;
    }
    return map;
  }

  sta->ensureGraph();
  sta->ensureLevelized();
  sta->updateTiming(false);

  const sta::MinMax* min_max = sta::MinMax::max();
  const int nworst = std::max(1, nworst_paths);

  sta::PathEndSeq path_ends
      = sta->findPathEnds(nullptr, nullptr, nullptr, min_max, nworst, true, -sta::INF, 0, true);

  std::set<std::pair<int, int>> seen_pairs;

  for (sta::PathEnd* path_end : path_ends) {
    if (!path_end) {
      continue;
    }
    sta::Path* path = path_end->path();
    if (!path) {
      continue;
    }

    const sta::Slack slack = path_end->slack();
    const double slack_value = static_cast<double>(slack);
    if (slack_value >= 0.0) {
      continue;
    }

    const double path_weight = -slack_value;
    sta::PathExpanded expand(path, sta);
    std::vector<odb::dbInst*> path_macros;

    sta::dbNetwork* network = sta->getDbNetwork();
    for (size_t i = 0; i < expand.size(); i++) {
      const sta::Path* path_pt = expand.path(i);
      sta::Pin* pin = path_pt->pin(sta);
      if (!pin) {
        continue;
      }
      sta::Instance* inst = network->instance(pin);
      if (!inst) {
        continue;
      }
      odb::dbInst* db_inst = nullptr;
      odb::dbITerm* db_iterm = nullptr;
      odb::dbBTerm* db_bterm = nullptr;
      network->staToDb(inst, db_inst, db_iterm, db_bterm);
      if (!db_inst || !isMacroInst(db_inst)) {
        continue;
      }
      map.path_criticality[db_inst] += path_weight;
      if (path_macros.empty() || path_macros.back() != db_inst) {
        path_macros.push_back(db_inst);
      }
    }

    for (size_t i = 0; i < path_macros.size(); ++i) {
      for (size_t j = i + 1; j < path_macros.size(); ++j) {
        odb::dbInst* a = path_macros[i];
        odb::dbInst* b = path_macros[j];
        const int ia = map.macro_index.at(a);
        const int ib = map.macro_index.at(b);
        const int lo = std::min(ia, ib);
        const int hi = std::max(ia, ib);
        if (seen_pairs.count({lo, hi})) {
          continue;
        }
        seen_pairs.insert({lo, hi});

        const double dist = macroPairDistance(a, b);
        MacroAffinityPair pair;
        pair.macro_idx_a = lo;
        pair.macro_idx_b = hi;
        pair.target_dist_dbu = static_cast<float>(dist * 0.85);
        pair.weight = static_cast<float>(path_weight);
        map.affinity_pairs.push_back(pair);
      }
    }
  }

  // Supplement with negative-slack nets for macros not on reported paths.
  for (odb::dbInst* inst : map.macros) {
    if (map.path_criticality[inst] > 0.0) {
      continue;
    }
    double score = 0.0;
    for (odb::dbITerm* iterm : inst->getITerms()) {
      odb::dbNet* net = iterm->getNet();
      if (!net || net->getITerms().size() < 2) {
        continue;
      }
      const sta::Slack net_slack = sta->slack(net, min_max);
      const double sv = static_cast<double>(net_slack);
      if (sv < 0.0) {
        const int degree = static_cast<int>(net->getITerms().size());
        score += (-sv) / std::max(1, degree - 1);
      }
    }
    map.path_criticality[inst] = score;
  }

  if (logger) {
    logger->info(utl::MPL,
                 110,
                 "TAHPP: {} macros, {} affinity pairs from {} worst paths",
                 map.macros.size(),
                 map.affinity_pairs.size(),
                 path_ends.size());
  }

  return map;
}

double TimingCriticalityMap::getMacroScore(odb::dbInst* inst) const
{
  const auto it = path_criticality.find(inst);
  if (it == path_criticality.end()) {
    return 0.0;
  }
  return it->second;
}

double TimingCriticalityMap::computeGlobalAffinityCost(odb::dbBlock* block) const
{
  (void) block;
  double cost = 0.0;
  for (const MacroAffinityPair& pair : affinity_pairs) {
    if (pair.macro_idx_a < 0 || pair.macro_idx_b < 0
        || pair.macro_idx_a >= static_cast<int>(macros.size())
        || pair.macro_idx_b >= static_cast<int>(macros.size())) {
      continue;
    }
    odb::dbInst* a = macros[pair.macro_idx_a];
    odb::dbInst* b = macros[pair.macro_idx_b];
    const double dist = macroPairDistance(a, b);
    const double excess = std::max(0.0, dist - pair.target_dist_dbu);
    cost += pair.weight * excess;
  }
  return cost;
}

double TimingCriticalityMap::computeNetTimingPressure(odb::dbBlock* block,
                                                    sta::dbSta* sta) const
{
  if (!sta) {
    return 0.0;
  }
  double total = 0.0;
  for (odb::dbInst* inst : macros) {
    const sta::MinMax* min_max = sta::MinMax::max();
    for (odb::dbITerm* iterm : inst->getITerms()) {
      odb::dbNet* net = iterm->getNet();
      if (!net || net->getITerms().size() < 2) {
        continue;
      }
      const sta::Slack slack = sta->slack(net, min_max);
      const double sv = static_cast<double>(slack);
      if (sv < 0.0) {
        const int degree = static_cast<int>(net->getITerms().size());
        total += (-sv) / std::max(1, degree - 1);
      }
    }
  }
  (void) block;
  return total;
}

}  // namespace tahpp
}  // namespace mpl
