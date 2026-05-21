// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#pragma once

#include <cstdint>
#include <vector>

namespace odb {
class dbBlock;
class dbInst;
class dbOrientType;
}  // namespace odb

namespace mpl {
namespace tahpp {

struct MacroPlacement
{
  int macro_idx = -1;
  int x_dbu = 0;
  int y_dbu = 0;
  int orient = 0;  // odb::dbOrientType::Value
};

struct PlacementSolution
{
  std::vector<MacroPlacement> macros;
  double fitness = 0.0;
  double hpwl = 0.0;
  double affinity_cost = 0.0;
  double overlap_penalty = 0.0;
  double outline_penalty = 0.0;
  double wns_proxy = 0.0;  // negative affinity (lower is better timing)
  double floorplan_area = 0.0;
  bool legal = true;

  void captureFromBlock(const std::vector<odb::dbInst*>& macro_list,
                        odb::dbBlock* block);
  void applyToBlock(const std::vector<odb::dbInst*>& macro_list,
                    odb::dbBlock* block) const;
};

struct TahppConfig
{
  bool use_gpu = false;
  int pop_size = 32;
  int num_islands = 4;
  int num_generations = 20;
  int init_population = 256;
  int elite_count = 4;
  int gpu_batch_size = 64;
  int migration_interval = 10;
  int sta_refine_iterations = 80;
  float timing_weight = 1.0f;
  float wns_weight = 0.5f;
  float area_weight = 0.25f;
  float wl_weight = 0.25f;
  float candidate_radius_um = 50.0f;
  int hybrid_iterations = 6;
  unsigned seed = 42;
};

}  // namespace tahpp
}  // namespace mpl
