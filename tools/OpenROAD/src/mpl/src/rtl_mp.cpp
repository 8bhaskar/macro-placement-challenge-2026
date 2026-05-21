// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#include "mpl/rtl_mp.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "MplObserver.h"
#include "db_sta/dbSta.hh"
#include "hier_rtlmp.h"
#include "object.h"
#include "odb/db.h"
#include "odb/geom.h"
#include "sta/MinMax.hh"
#include "utl/Logger.h"

namespace mpl {
using std::string;
using utl::MPL;

class Snapper;

MacroPlacer::MacroPlacer(sta::dbNetwork* network,
                         odb::dbDatabase* db,
                         sta::dbSta* sta,
                         utl::Logger* logger,
                         par::PartitionMgr* tritonpart)
{
  hier_rtlmp_ = std::make_unique<HierRTLMP>(network, db, logger, tritonpart);
  logger_ = logger;
  db_ = db;
  sta_ = sta;
}

MacroPlacer::~MacroPlacer() = default;

bool MacroPlacer::place(const int num_threads,
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
                        const char* report_directory,
                        const bool keep_clustering_data)
{
  hier_rtlmp_->init();
  hier_rtlmp_->setClusterSize(
      max_num_macro, min_num_macro, max_num_inst, min_num_inst);
  hier_rtlmp_->setClusterSizeTolerance(tolerance);
  hier_rtlmp_->setMaxNumLevel(max_num_level);
  hier_rtlmp_->setClusterSizeRatioPerLevel(coarsening_ratio);
  hier_rtlmp_->setLargeNetThreshold(large_net_threshold);
  hier_rtlmp_->setDefaultHalo(halo_width, halo_height);
  hier_rtlmp_->setGlobalFence(global_fence);
  hier_rtlmp_->setAreaWeight(area_weight);
  hier_rtlmp_->setOutlineWeight(outline_weight);
  hier_rtlmp_->setWirelengthWeight(wirelength_weight);
  hier_rtlmp_->setGuidanceWeight(guidance_weight);
  hier_rtlmp_->setFenceWeight(fence_weight);
  hier_rtlmp_->setBoundaryWeight(boundary_weight);
  hier_rtlmp_->setNotchWeight(notch_weight);
  hier_rtlmp_->setMacroBlockageWeight(macro_blockage_weight);
  hier_rtlmp_->setTargetUtil(target_util);
  hier_rtlmp_->setMinAR(min_ar);
  hier_rtlmp_->setReportDirectory(report_directory);
  hier_rtlmp_->setNumThreads(num_threads);
  hier_rtlmp_->setKeepClusteringData(keep_clustering_data);
  hier_rtlmp_->setGuidanceRegions(guidance_regions_);

  hier_rtlmp_->run();

  return true;
}

void MacroPlacer::placeMacro(odb::dbInst* inst,
                             const float& x_origin,
                             const float& y_origin,
                             const odb::dbOrientType& orientation,
                             const bool exact,
                             const bool allow_overlap)
{
  odb::dbBlock* block = inst->getBlock();

  const int x1 = block->micronsToDbu(x_origin);
  const int y1 = block->micronsToDbu(y_origin);
  const int x2 = x1 + inst->getBBox()->getDX();
  const int y2 = y1 + inst->getBBox()->getDY();

  odb::Rect macro_new_bbox(x1, y1, x2, y2);
  odb::Rect core_area = inst->getBlock()->getCoreArea();

  if (!core_area.contains(macro_new_bbox)) {
    logger_->error(MPL,
                   34,
                   "Cannot place {} at ({}, {}) ({}, {}), outside of the core "
                   "({}, {}) ({}, {}).",
                   inst->getName(),
                   block->dbuToMicrons(macro_new_bbox.xMin()),
                   block->dbuToMicrons(macro_new_bbox.yMin()),
                   block->dbuToMicrons(macro_new_bbox.xMax()),
                   block->dbuToMicrons(macro_new_bbox.yMax()),
                   block->dbuToMicrons(core_area.xMin()),
                   block->dbuToMicrons(core_area.yMin()),
                   block->dbuToMicrons(core_area.xMax()),
                   block->dbuToMicrons(core_area.yMax()));
  }

  // Orientation must be set before location so we don't end up flipping
  // and misplacing the macro.
  inst->setOrient(orientation);
  inst->setLocation(x1, y1);

  if (orientation.isRightAngleRotation()) {
    logger_->warn(MPL,
                  36,
                  "Orientation {} specified for macro {} is a right angle "
                  "rotation. Snapping is not possible.",
                  orientation.getString(),
                  inst->getName());
  } else if (!exact) {
    Snapper snapper(logger_, inst);
    snapper.snapMacro();
  }

  if (!allow_overlap) {
    std::vector<odb::dbInst*> overlapped_macros = findOverlappedMacros(inst);
    if (!overlapped_macros.empty()) {
      std::string overlapped_macros_names;
      for (odb::dbInst* overlapped_macro : overlapped_macros) {
        overlapped_macros_names
            += fmt::format(" {}", overlapped_macro->getName());
      }

      logger_->error(MPL,
                     41,
                     "Couldn't place {}. Found overlap with other macros:{}.",
                     inst->getName(),
                     overlapped_macros_names);
    }
  }

  inst->setPlacementStatus(odb::dbPlacementStatus::LOCKED);

  logger_->info(MPL,
                35,
                "Macro {} placed. Bounding box ({:.3f}um, {:.3f}um), "
                "({:.3f}um, {:.3f}um). Orientation {}",
                inst->getName(),
                block->dbuToMicrons(inst->getBBox()->xMin()),
                block->dbuToMicrons(inst->getBBox()->yMin()),
                block->dbuToMicrons(inst->getBBox()->xMax()),
                block->dbuToMicrons(inst->getBBox()->yMax()),
                orientation.getString());
}

std::vector<odb::dbInst*> MacroPlacer::findOverlappedMacros(odb::dbInst* macro)
{
  std::vector<odb::dbInst*> overlapped_macros;
  odb::dbBlock* block = macro->getBlock();
  const odb::Rect& source_macro_bbox = macro->getBBox()->getBox();

  for (odb::dbInst* inst : block->getInsts()) {
    if (!inst->isBlock() || !inst->isPlaced()) {
      continue;
    }

    const odb::Rect& target_macro_bbox = inst->getBBox()->getBox();
    if (source_macro_bbox.overlaps(target_macro_bbox)) {
      overlapped_macros.push_back(inst);
    }
  }

  return overlapped_macros;
}

void MacroPlacer::addGuidanceRegion(odb::dbInst* macro, odb::Rect region)
{
  odb::dbBlock* block = db_->getChip()->getBlock();
  const odb::Rect& core = block->getCoreArea();

  if (!core.contains(region)) {
    logger_->error(MPL,
                   42,
                   "Specified guidance region ({}, {}) ({}, {}) for the macro "
                   "{} is outside of the core ({}, {}) ({}, {}).",
                   region.xMin(),
                   region.yMin(),
                   region.xMax(),
                   region.yMax(),
                   macro->getName(),
                   block->dbuToMicrons(core.xMin()),
                   block->dbuToMicrons(core.yMin()),
                   block->dbuToMicrons(core.xMax()),
                   block->dbuToMicrons(core.yMax()));
  }

  if (guidance_regions_.find(macro) != guidance_regions_.end()) {
    logger_->warn(
        MPL, 44, "Overwriting guidance region for macro {}", macro->getName());
  }

  guidance_regions_[macro] = region;
}

void MacroPlacer::setMacroHalo(odb::dbInst* macro,
                               int left,
                               int bottom,
                               int right,
                               int top)
{
  hier_rtlmp_->setMacroHalo(macro, left, bottom, right, top);
}

void MacroPlacer::setMacroPlacementFile(const std::string& file_name)
{
  hier_rtlmp_->setMacroPlacementFile(file_name);
}

void MacroPlacer::setDebug(std::unique_ptr<MplObserver>& graphics)
{
  hier_rtlmp_->setDebug(graphics);
}

void MacroPlacer::setDebugShowBundledNets(bool show_bundled_nets)
{
  hier_rtlmp_->setDebugShowBundledNets(show_bundled_nets);
}
void MacroPlacer::setDebugShowClustersIds(bool show_clusters_ids)
{
  hier_rtlmp_->setDebugShowClustersIds(show_clusters_ids);
}

void MacroPlacer::setDebugSkipSteps(bool skip_steps)
{
  hier_rtlmp_->setDebugSkipSteps(skip_steps);
}

void MacroPlacer::setDebugOnlyFinalResult(bool only_final_result)
{
  hier_rtlmp_->setDebugOnlyFinalResult(only_final_result);
}

void MacroPlacer::setDebugTargetClusterId(const int target_cluster_id)
{
  hier_rtlmp_->setDebugTargetClusterId(target_cluster_id);
}

namespace {

double computeMacroNetCriticality(odb::dbInst* inst)
{
  double score = 0.0;
  for (odb::dbITerm* iterm : inst->getITerms()) {
    auto* net = iterm->getNet();
    if (!net) {
      continue;
    }
    const int degree = static_cast<int>(net->getITerms().size());
    score += 1.0 / std::max(1, degree - 1);
  }
  return score;
}

bool isCandidateWithinBlock(const odb::Rect& candidate_bbox,
                            const odb::Rect& block_bbox)
{
  return candidate_bbox.xMin() >= block_bbox.xMin()
      && candidate_bbox.yMin() >= block_bbox.yMin()
      && candidate_bbox.xMax() <= block_bbox.xMax()
      && candidate_bbox.yMax() <= block_bbox.yMax();
}

double computeMacroTimingPressure(odb::dbInst* inst, sta::dbSta* sta)
{
  if (!sta) {
    return 0.0;
  }

  const sta::MinMax* min_max = sta::MinMax::max();
  double score = 0.0;

  for (odb::dbITerm* iterm : inst->getITerms()) {
    auto* net = iterm->getNet();
    if (!net || net->getITerms().size() < 2) {
      continue;
    }

    const sta::Slack slack = sta->slack(net, min_max);
    const double slack_value = static_cast<double>(slack);
    if (slack_value >= 0.0) {
      continue;
    }

    const int degree = static_cast<int>(net->getITerms().size());
    const double criticality = 1.0 / std::max(1, degree - 1);
    score += (-slack_value) * criticality;
  }

  return score;
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
  const odb::Rect block_bbox = block->getBBox();
  const int radius_dbu = std::max(1, block->micronsToDbu(candidate_radius));
  const double scale = timing_weight > 0.0 ? timing_weight : 1.0;

  if (sta_) {
    sta_->updateTiming(false);
  }

  std::vector<std::tuple<double, double, odb::dbInst*>> macro_scores;
  for (odb::dbInst* inst : block->getInsts()) {
    if (!inst->isCore() || inst->isFixed()) {
      continue;
    }
    const double net_crit = computeMacroNetCriticality(inst);
    const double timing_pressure = computeMacroTimingPressure(inst, sta_);
    macro_scores.emplace_back(timing_pressure + (net_crit * 0.1), timing_pressure, inst);
  }

  if (macro_scores.empty()) {
    return true;
  }

  std::sort(macro_scores.begin(), macro_scores.end(),
            [](const auto& a, const auto& b) {
              return std::get<0>(a) > std::get<0>(b);
            });

  double current_hpwl = computeNetHPWLCost(block);
  double current_timing = 0.0;
  if (sta_) {
    for (odb::dbInst* inst : block->getInsts()) {
      if (inst->isCore() && !inst->isFixed()) {
        current_timing += computeMacroTimingPressure(inst, sta_);
      }
    }
  }
  double current_cost = current_hpwl + (scale * current_timing);
  const int active_macros = std::min(static_cast<int>(macro_scores.size()),
                                     std::max(8, static_cast<int>(macro_scores.size() / 2)));

  for (int iteration = 0; iteration < hybrid_iterations; ++iteration) {
    bool improved = false;

    for (int macro_index = 0; macro_index < active_macros; ++macro_index) {
      odb::dbInst* inst = std::get<2>(macro_scores[macro_index]);
      const odb::Rect original_bbox = inst->getBBox()->getBox();
      const int original_x = original_bbox.xMin();
      const int original_y = original_bbox.yMin();
      const odb::dbOrientType original_orient = inst->getOrient();

      const std::vector<std::pair<int, int>> moves = {
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

      for (const auto& [dx, dy] : moves) {
        const odb::Rect candidate_bbox(
            original_bbox.xMin() + dx,
            original_bbox.yMin() + dy,
            original_bbox.xMax() + dx,
            original_bbox.yMax() + dy);

        if (!isCandidateWithinBlock(candidate_bbox, block_bbox)) {
          continue;
        }

        inst->setLocation(original_x + dx, original_y + dy);
        inst->setOrient(original_orient);

        if (!findOverlappedMacros(inst).empty()) {
          continue;
        }

        const double candidate_hpwl = computeNetHPWLCost(block);
        const double candidate_timing = computeMacroTimingPressure(inst, sta_);
        const double candidate_cost = candidate_hpwl + (scale * candidate_timing);

        if (candidate_cost < best_cost) {
          best_cost = candidate_cost;
          best_x = original_x + dx;
          best_y = original_y + dy;
        }
      }

      inst->setLocation(original_x, original_y);
      inst->setOrient(original_orient);

      if (best_x != original_x || best_y != original_y) {
        inst->setLocation(best_x, best_y);
        if (findOverlappedMacros(inst).empty()) {
          current_cost = best_cost;
          improved = true;
          logger_->info(MPL,
                        100,
                        "Hybrid macro move improved cost for {} to {:.3f}",
                        inst->getName(),
                        current_cost);
          break;
        }
        inst->setLocation(original_x, original_y);
      }
    }

    if (!improved) {
      break;
    }
  }

  if (sta_) {
    sta::Slack worst_slack = sta_->worstSlack(sta::MinMax::max());
    logger_->info(MPL,
                  101,
                  "Hybrid macro placement complete. Worst slack = {:.3f}",
                  static_cast<double>(worst_slack));
  }

  return true;
}

}  // namespace mpl
