// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#pragma once

#include <unordered_map>
#include <vector>

namespace odb {
class dbBlock;
class dbInst;
}  // namespace odb

namespace sta {
class dbSta;
}  // namespace sta

namespace utl {
class Logger;
}

namespace mpl {
namespace tahpp {

struct MacroAffinityPair
{
  int macro_idx_a = -1;
  int macro_idx_b = -1;
  float target_dist_dbu = 0.0f;
  float weight = 1.0f;
};

// Path-based timing criticality and macro affinity penalties for surrogate cost.
class TimingCriticalityMap
{
 public:
  std::vector<odb::dbInst*> macros;
  std::unordered_map<odb::dbInst*, int> macro_index;
  std::unordered_map<odb::dbInst*, double> path_criticality;
  std::vector<MacroAffinityPair> affinity_pairs;

  static TimingCriticalityMap build(odb::dbBlock* block,
                                    sta::dbSta* sta,
                                    utl::Logger* logger,
                                    int nworst_paths = 20);

  double getMacroScore(odb::dbInst* inst) const;
  double computeGlobalAffinityCost(odb::dbBlock* block) const;
  double computeNetTimingPressure(odb::dbBlock* block, sta::dbSta* sta) const;
};

}  // namespace tahpp
}  // namespace mpl
