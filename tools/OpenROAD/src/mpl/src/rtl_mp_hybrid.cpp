// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#include "mpl/rtl_mp.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>
#include <vector>

#include "db_sta/dbSta.hh"
#include "odb/db.h"
#include "sta/MinMax.hh"
#include "tahpp/placement_solution.h"
#include "tahpp/tahpp_optimizer.h"
#include "tahpp/timing_criticality.h"
#include "utl/Logger.h"

namespace mpl {
namespace {

bool isMacroInst(odb::dbInst* inst)
{
  return inst && (inst->isBlock() || inst->isCore()) && !inst->isFixed();
}

bool isCandidateWithinBlock(const odb::Rect& candidate_bbox,
                            const odb::Rect& block_bbox)
{
  return candidate_bbox.xMin() >= block_bbox.xMin()
         && candidate_bbox.yMin() >= block_bbox.yMin()
         && candidate_bbox.xMax() <= block_bbox.xMax()
         && candidate_bbox.yMax() <= block_bbox.yMax();
}

double computeNetHPWLCost(odb::dbBlock* block)
{
  double cost = 0.0;
  for (odb::dbNet* net : block->getNets()) {
    const auto& iterms = net->getITerms();
    if (iterms.size() < 2) {
      continue;
    }

    int min_x = std::numeric_limits<int>::max();
    int min_y = std::numeric_limits<int>::max();
    int max_x = std::numeric_limits<int>::min();
    int max_y = std::numeric_limits<int>::min();
    bool valid = false;

    for (odb::dbITerm* iterm : iterms) {
      auto* inst = iterm->getInst();
      if (!inst) {
        continue;
      }
      const odb::Rect bbox = inst->getBBox()->getBox();
      const int x = bbox.xCenter();
      const int y = bbox.yCenter();
      min_x = std::min(min_x, x);
      min_y = std::min(min_y, y);
      max_x = std::max(max_x, x);
      max_y = std::max(max_y, y);
      valid = true;
    }

    if (valid) {
      cost += static_cast<double>(max_x - min_x) + static_cast<double>(max_y - min_y);
    }
  }
  return cost;
}

std::vector<odb::dbOrientType::Value> legalOrientations()
{
  return {odb::dbOrientType::R0,
          odb::dbOrientType::R180,
          odb::dbOrientType::MY,
          odb::dbOrientType::MX};
}

double computeGlobalHybridCost(odb::dbBlock* block,
                               const tahpp::TimingCriticalityMap& criticality,
                               const double timing_scale)
{
  return computeNetHPWLCost(block)
         + timing_scale * criticality.computeGlobalAffinityCost(block);
}

void runImprovedHybridLocalSearch(
    MacroPlacer* placer,
    odb::dbBlock* block,
    sta::dbSta* sta,
    utl::Logger* logger,
    const tahpp::TimingCriticalityMap& criticality,
    const float timing_weight,
    const int hybrid_iterations,
    int radius_dbu)
{
  const odb::Rect block_bbox = block->getBBox();
  const double scale = timing_weight > 0.0 ? timing_weight : 1.0;
  const auto orients = legalOrientations();

  std::vector<std::tuple<double, odb::dbInst*>> macro_scores;
  for (odb::dbInst* inst : criticality.macros) {
    macro_scores.emplace_back(criticality.getMacroScore(inst), inst);
  }
  if (macro_scores.empty()) {
    return;
  }

  std::sort(macro_scores.begin(),
            macro_scores.end(),
            [](const auto& a, const auto& b) {
              return std::get<0>(a) > std::get<0>(b);
            });

  double current_cost = computeGlobalHybridCost(block, criticality, scale);
  const int active_macros = std::min(
      static_cast<int>(macro_scores.size()),
      std::max(8, static_cast<int>(macro_scores.size() / 2)));

  for (int iteration = 0; iteration < hybrid_iterations; ++iteration) {
    bool improved = false;

    for (int macro_index = 0; macro_index < active_macros; ++macro_index) {
      odb::dbInst* inst = std::get<1>(macro_scores[macro_index]);
      const odb::Rect original_bbox = inst->getBBox()->getBox();
      const int original_x = original_bbox.xMin();
      const int original_y = original_bbox.yMin();
      const odb::dbOrientType original_orient = inst->getOrient();

      std::vector<std::pair<int, int>> moves = {
          {0, 0},
          {radius_dbu, 0},
          {-radius_dbu, 0},
          {0, radius_dbu},
          {0, -radius_dbu},
          {radius_dbu, radius_dbu},
          {-radius_dbu, radius_dbu},
          {radius_dbu, -radius_dbu},
          {-radius_dbu, -radius_dbu}};

      double best_cost = current_cost;
      int best_x = original_x;
      int best_y = original_y;
      odb::dbOrientType best_orient = original_orient;

      for (const auto& [dx, dy] : moves) {
        for (const auto orient_val : orients) {
          const odb::dbOrientType orient(orient_val);
          const odb::Rect candidate_bbox(original_bbox.xMin() + dx,
                                         original_bbox.yMin() + dy,
                                         original_bbox.xMax() + dx,
                                         original_bbox.yMax() + dy);

          if (!isCandidateWithinBlock(candidate_bbox, block_bbox)) {
            continue;
          }

          inst->setLocation(original_x + dx, original_y + dy);
          inst->setOrient(orient);

          if (!placer->findOverlappedMacros(inst).empty()) {
            continue;
          }

          const double candidate_cost
              = computeGlobalHybridCost(block, criticality, scale);

          if (candidate_cost < best_cost) {
            best_cost = candidate_cost;
            best_x = original_x + dx;
            best_y = original_y + dy;
            best_orient = orient;
          }
        }
      }

      inst->setLocation(original_x, original_y);
      inst->setOrient(original_orient);

      if (best_x != original_x || best_y != original_y
          || best_orient != original_orient) {
        inst->setLocation(best_x, best_y);
        inst->setOrient(best_orient);
        if (placer->findOverlappedMacros(inst).empty()) {
          current_cost = best_cost;
          improved = true;
          radius_dbu = std::min(radius_dbu * 2, radius_dbu + radius_dbu / 2 + 1);
          logger->info(utl::MPL,
                       100,
                       "Hybrid macro move improved cost for {} to {:.3f}",
                       inst->getName(),
                       current_cost);
          break;
        }
        inst->setLocation(original_x, original_y);
        inst->setOrient(original_orient);
      } else {
        radius_dbu = std::max(1, radius_dbu / 2);
      }
    }

    if (!improved) {
      break;
    }
  }

  if (sta) {
    sta->updateTiming(false);
    const sta::Slack worst_slack = sta->worstSlack(sta::MinMax::max());
    logger->info(utl::MPL,
                 101,
                 "Hybrid macro placement complete. Worst slack = {:.3f}",
                 static_cast<double>(worst_slack));
  }
}

std::vector<odb::dbInst*> collectMacros(odb::dbBlock* block)
{
  std::vector<odb::dbInst*> macros;
  for (odb::dbInst* inst : block->getInsts()) {
    if (isMacroInst(inst)) {
      macros.push_back(inst);
    }
  }
  return macros;
}

}  // namespace

bool MacroPlacer::placeHybrid(const int num_threads,
                              const int max_num_macro,
                              const int min_num_macro,
                              const int max_num_inst,
                              const int min_num_inst,
                              const float tolerance,
                              const int max_num_level,
                              const float coarsening_ratio,
                              const int large_net_threshold,
                              const int halo_width,
                              const int halo_height,
                              const odb::Rect global_fence,
                              const float area_weight,
                              const float outline_weight,
                              const float wirelength_weight,
                              const float guidance_weight,
                              const float fence_weight,
                              const float boundary_weight,
                              const float notch_weight,
                              const float macro_blockage_weight,
                              const float target_util,
                              const float min_ar,
                              const float timing_weight,
                              const int hybrid_iterations,
                              const float candidate_radius,
                              const char* report_directory,
                              const bool keep_clustering_data)
{
  if (!place(num_threads,
             max_num_macro,
             min_num_macro,
             max_num_inst,
             min_num_inst,
             tolerance,
             max_num_level,
             coarsening_ratio,
             large_net_threshold,
             halo_width,
             halo_height,
             global_fence,
             area_weight,
             outline_weight,
             wirelength_weight,
             guidance_weight,
             fence_weight,
             boundary_weight,
             notch_weight,
             macro_blockage_weight,
             target_util,
             min_ar,
             report_directory,
             keep_clustering_data)) {
    return false;
  }

  auto* block = db_->getChip()->getBlock();
  tahpp::TimingCriticalityMap criticality
      = tahpp::TimingCriticalityMap::build(block, sta_, logger_, 20);

  if (criticality.macros.empty()) {
    return true;
  }

  const int radius_dbu = std::max(1, block->micronsToDbu(candidate_radius));
  runImprovedHybridLocalSearch(this,
                               block,
                               sta_,
                               logger_,
                               criticality,
                               timing_weight,
                               hybrid_iterations,
                               radius_dbu);
  return true;
}

bool MacroPlacer::placeTahpp(const int num_threads,
                             const int max_num_macro,
                             const int min_num_macro,
                             const int max_num_inst,
                             const int min_num_inst,
                             const float tolerance,
                             const int max_num_level,
                             const float coarsening_ratio,
                             const int large_net_threshold,
                             const int halo_width,
                             const int halo_height,
                             const odb::Rect global_fence,
                             const float area_weight,
                             const float outline_weight,
                             const float wirelength_weight,
                             const float guidance_weight,
                             const float fence_weight,
                             const float boundary_weight,
                             const float notch_weight,
                             const float macro_blockage_weight,
                             const float target_util,
                             const float min_ar,
                             const float timing_weight,
                             const int hybrid_iterations,
                             const float candidate_radius,
                             const bool use_gpu,
                             const int pop_size,
                             const int num_islands,
                             const int num_generations,
                             const int init_population,
                             const int elite_count,
                             const int gpu_batch_size,
                             const float wns_weight,
                             const float area_weight_fitness,
                             const float wl_weight,
                             const char* report_directory,
                             const bool keep_clustering_data)
{
  if (!place(num_threads,
             max_num_macro,
             min_num_macro,
             max_num_inst,
             min_num_inst,
             tolerance,
             max_num_level,
             coarsening_ratio,
             large_net_threshold,
             halo_width,
             halo_height,
             global_fence,
             area_weight,
             outline_weight,
             wirelength_weight,
             guidance_weight,
             fence_weight,
             boundary_weight,
             notch_weight,
             macro_blockage_weight,
             target_util,
             min_ar,
             report_directory,
             keep_clustering_data)) {
    return false;
  }

  auto* block = db_->getChip()->getBlock();
  std::vector<odb::dbInst*> macro_list = collectMacros(block);
  if (macro_list.empty()) {
    return true;
  }

  tahpp::TimingCriticalityMap criticality
      = tahpp::TimingCriticalityMap::build(block, sta_, logger_, 20);

  tahpp::TahppConfig config;
  config.use_gpu = use_gpu;
  config.pop_size = std::max(4, pop_size);
  config.num_islands = std::max(1, num_islands);
  config.num_generations = std::max(1, num_generations);
  config.init_population = std::max(config.pop_size, init_population);
  config.elite_count = std::max(1, elite_count);
  config.gpu_batch_size = std::max(1, gpu_batch_size);
  config.timing_weight = timing_weight;
  config.hybrid_iterations = hybrid_iterations;
  config.candidate_radius_um = candidate_radius;
  config.wns_weight = wns_weight;
  config.area_weight = area_weight_fitness;
  config.wl_weight = wl_weight;

  tahpp::TahppOptimizer optimizer(this, db_, sta_, logger_);
  return optimizer.run(config, macro_list, block, criticality);
}

}  // namespace mpl
