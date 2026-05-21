// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#include "tahpp/gpu_batch_eval.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "odb/db.h"

#ifdef BUILD_TAHPP_GPU
#include <cuda_runtime.h>
#endif

namespace mpl {
namespace tahpp {

namespace {

double computeHpwlForSolution(const std::vector<odb::dbInst*>& macro_list,
                              const PlacementSolution& sol)
{
  double cost = 0.0;
  for (const MacroPlacement& mp : sol.macros) {
    if (mp.macro_idx < 0 || mp.macro_idx >= static_cast<int>(macro_list.size())) {
      continue;
    }
    odb::dbInst* inst = macro_list[mp.macro_idx];
    for (odb::dbITerm* iterm : inst->getITerms()) {
      odb::dbNet* net = iterm->getNet();
      if (!net) {
        continue;
      }
      int min_x = std::numeric_limits<int>::max();
      int min_y = std::numeric_limits<int>::max();
      int max_x = std::numeric_limits<int>::min();
      int max_y = std::numeric_limits<int>::min();
      bool valid = false;
      for (odb::dbITerm* niterm : net->getITerms()) {
        odb::dbInst* ninst = niterm->getInst();
        if (!ninst) {
          continue;
        }
        int cx;
        int cy;
        if (ninst->isBlock() && ninst != inst) {
          const auto mit = std::find_if(
              sol.macros.begin(),
              sol.macros.end(),
              [ninst, &macro_list](const MacroPlacement& m) {
                return m.macro_idx >= 0
                       && m.macro_idx < static_cast<int>(macro_list.size())
                       && macro_list[m.macro_idx] == ninst;
              });
          if (mit != sol.macros.end()) {
            odb::dbMaster* master = ninst->getMaster();
            cx = mit->x_dbu + master->getWidth() / 2;
            cy = mit->y_dbu + master->getHeight() / 2;
          } else {
            const odb::Rect bbox = ninst->getBBox()->getBox();
            cx = bbox.xCenter();
            cy = bbox.yCenter();
          }
        } else {
          const odb::Rect bbox = ninst->getBBox()->getBox();
          cx = bbox.xCenter();
          cy = bbox.yCenter();
        }
        if (ninst == inst) {
          odb::dbMaster* master = inst->getMaster();
          cx = mp.x_dbu + master->getWidth() / 2;
          cy = mp.y_dbu + master->getHeight() / 2;
        }
        min_x = std::min(min_x, cx);
        min_y = std::min(min_y, cy);
        max_x = std::max(max_x, cx);
        max_y = std::max(max_y, cy);
        valid = true;
      }
      if (valid) {
        cost += static_cast<double>(max_x - min_x + max_y - min_y);
      }
    }
  }
  return cost;
}

double computeAffinityForSolution(const TimingCriticalityMap& crit,
                                const std::vector<odb::dbInst*>& macro_list,
                                const PlacementSolution& sol)
{
  double cost = 0.0;
  for (const MacroAffinityPair& pair : crit.affinity_pairs) {
    if (pair.macro_idx_a < 0 || pair.macro_idx_b < 0
        || pair.macro_idx_a >= static_cast<int>(sol.macros.size())
        || pair.macro_idx_b >= static_cast<int>(sol.macros.size())) {
      continue;
    }
    const MacroPlacement& ma = sol.macros[pair.macro_idx_a];
    const MacroPlacement& mb = sol.macros[pair.macro_idx_b];
    odb::dbInst* inst_a = macro_list[ma.macro_idx];
    odb::dbInst* inst_b = macro_list[mb.macro_idx];
    const int ax = ma.x_dbu + inst_a->getMaster()->getWidth() / 2;
    const int ay = ma.y_dbu + inst_a->getMaster()->getHeight() / 2;
    const int bx = mb.x_dbu + inst_b->getMaster()->getWidth() / 2;
    const int by = mb.y_dbu + inst_b->getMaster()->getHeight() / 2;
    const double dist = std::hypot(static_cast<double>(ax - bx),
                                   static_cast<double>(ay - by));
    const double excess = std::max(0.0, dist - pair.target_dist_dbu);
    cost += pair.weight * excess;
  }
  return cost;
}

double computeOverlapPenalty(const std::vector<odb::dbInst*>& macro_list,
                             const PlacementSolution& sol)
{
  double penalty = 0.0;
  const size_t n = sol.macros.size();
  for (size_t i = 0; i < n; ++i) {
    const MacroPlacement& mi = sol.macros[i];
    odb::dbInst* inst_i = macro_list[mi.macro_idx];
    const int w_i = inst_i->getMaster()->getWidth();
    const int h_i = inst_i->getMaster()->getHeight();
    const odb::Rect bi(mi.x_dbu, mi.y_dbu, mi.x_dbu + w_i, mi.y_dbu + h_i);
    for (size_t j = i + 1; j < n; ++j) {
      const MacroPlacement& mj = sol.macros[j];
      odb::dbInst* inst_j = macro_list[mj.macro_idx];
      const int w_j = inst_j->getMaster()->getWidth();
      const int h_j = inst_j->getMaster()->getHeight();
      const odb::Rect bj(mj.x_dbu, mj.y_dbu, mj.x_dbu + w_j, mj.y_dbu + h_j);
      if (bi.overlaps(bj)) {
        penalty += 1e6;
      }
    }
  }
  return penalty;
}

double computeOutlinePenalty(const odb::Rect& core_bbox,
                             const std::vector<odb::dbInst*>& macro_list,
                             const PlacementSolution& sol)
{
  double penalty = 0.0;
  for (const MacroPlacement& mp : sol.macros) {
    odb::dbInst* inst = macro_list[mp.macro_idx];
    const int w = inst->getMaster()->getWidth();
    const int h = inst->getMaster()->getHeight();
    const odb::Rect bbox(mp.x_dbu, mp.y_dbu, mp.x_dbu + w, mp.y_dbu + h);
    if (!core_bbox.contains(bbox)) {
      penalty += 1e5;
    }
  }
  return penalty;
}

double computeFloorplanArea(const std::vector<odb::dbInst*>& macro_list,
                            const PlacementSolution& sol)
{
  int min_x = std::numeric_limits<int>::max();
  int min_y = std::numeric_limits<int>::max();
  int max_x = std::numeric_limits<int>::min();
  int max_y = std::numeric_limits<int>::min();
  for (const MacroPlacement& mp : sol.macros) {
    odb::dbInst* inst = macro_list[mp.macro_idx];
    const int w = inst->getMaster()->getWidth();
    const int h = inst->getMaster()->getHeight();
    min_x = std::min(min_x, mp.x_dbu);
    min_y = std::min(min_y, mp.y_dbu);
    max_x = std::max(max_x, mp.x_dbu + w);
    max_y = std::max(max_y, mp.y_dbu + h);
  }
  if (min_x > max_x) {
    return 0.0;
  }
  return static_cast<double>(max_x - min_x) * static_cast<double>(max_y - min_y);
}

}  // namespace

GpuBatchEvaluator::GpuBatchEvaluator(const bool use_gpu) : use_gpu_(use_gpu)
{
#ifdef BUILD_TAHPP_GPU
  int device_count = 0;
  if (cudaGetDeviceCount(&device_count) == cudaSuccess && device_count > 0) {
    gpu_available_ = true;
  }
#else
  (void) use_gpu_;
  gpu_available_ = false;
#endif
  if (use_gpu_ && !gpu_available_) {
    use_gpu_ = false;
  }
}

GpuBatchEvaluator::~GpuBatchEvaluator() = default;

void GpuBatchEvaluator::evaluateBatch(
    const BatchEvalInput& input,
    std::vector<PlacementSolution>& solutions)
{
  if (use_gpu_ && gpu_available_) {
    evaluateBatchCuda(input, solutions);
    return;
  }
  evaluateBatchCpu(input, solutions);
}

void GpuBatchEvaluator::evaluateBatchCuda(
    const BatchEvalInput& input,
    std::vector<PlacementSolution>& solutions)
{
#ifdef BUILD_TAHPP_GPU
  evaluateBatchCpu(input, solutions);
#else
  evaluateBatchCpu(input, solutions);
#endif
}

void GpuBatchEvaluator::evaluateBatchCpu(
    const BatchEvalInput& input,
    std::vector<PlacementSolution>& solutions)
{
  for (PlacementSolution& sol : solutions) {
    sol.hpwl = computeHpwlForSolution(*input.macro_list, sol);
    sol.affinity_cost = computeAffinityForSolution(
        *input.criticality, *input.macro_list, sol);
    sol.overlap_penalty = computeOverlapPenalty(*input.macro_list, sol);
    sol.outline_penalty
        = computeOutlinePenalty(input.core_bbox, *input.macro_list, sol);
    sol.floorplan_area
        = computeFloorplanArea(*input.macro_list, sol);
    sol.wns_proxy = sol.affinity_cost;
    sol.legal = (sol.overlap_penalty == 0.0 && sol.outline_penalty == 0.0);
    sol.fitness = input.wl_weight * sol.hpwl
                  + input.wns_weight * sol.affinity_cost
                  + input.area_weight * sol.floorplan_area
                  + sol.overlap_penalty + sol.outline_penalty;
  }
}

}  // namespace tahpp
}  // namespace mpl
