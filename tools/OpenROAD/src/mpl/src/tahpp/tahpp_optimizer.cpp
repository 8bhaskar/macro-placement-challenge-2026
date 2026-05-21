// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#include "tahpp/tahpp_optimizer.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <tuple>

#include "db_sta/dbSta.hh"
#include "odb/db.h"
#include "sta/MinMax.hh"
#include "tahpp/gpu_batch_eval.h"
#include "utl/Logger.h"

namespace mpl {
namespace tahpp {

namespace {

bool isMacroInst(odb::dbInst* inst)
{
  return inst && (inst->isBlock() || inst->isCore()) && !inst->isFixed();
}

std::vector<odb::dbOrientType::Value> legalOrientations()
{
  return {odb::dbOrientType::R0,
          odb::dbOrientType::R180,
          odb::dbOrientType::MY,
          odb::dbOrientType::MX};
}

bool solutionLegalOnBlock(MacroPlacer* placer,
                          const std::vector<odb::dbInst*>& macro_list,
                          const PlacementSolution& sol,
                          odb::dbBlock* block)
{
  PlacementSolution tmp = sol;
  tmp.applyToBlock(macro_list, block);
  for (const MacroPlacement& mp : sol.macros) {
    if (mp.macro_idx < 0 || mp.macro_idx >= static_cast<int>(macro_list.size())) {
      return false;
    }
    if (!placer->findOverlappedMacros(macro_list[mp.macro_idx]).empty()) {
      return false;
    }
  }
  const odb::Rect core = block->getCoreArea();
  for (const MacroPlacement& mp : sol.macros) {
    odb::dbInst* inst = macro_list[mp.macro_idx];
    const int w = inst->getMaster()->getWidth();
    const int h = inst->getMaster()->getHeight();
    const odb::Rect bbox(mp.x_dbu, mp.y_dbu, mp.x_dbu + w, mp.y_dbu + h);
    if (!core.contains(bbox)) {
      return false;
    }
  }
  return true;
}

}  // namespace

TahppOptimizer::TahppOptimizer(MacroPlacer* placer,
                               odb::dbDatabase* db,
                               sta::dbSta* sta,
                               utl::Logger* logger)
    : placer_(placer), db_(db), sta_(sta), logger_(logger)
{
}

bool TahppOptimizer::run(const TahppConfig& config,
                         const std::vector<odb::dbInst*>& macro_list,
                         odb::dbBlock* block,
                         TimingCriticalityMap& criticality)
{
  if (macro_list.empty()) {
    return true;
  }

  PlacementSolution baseline;
  baseline.captureFromBlock(macro_list, block);

  std::vector<PlacementSolution> population;
  seedPopulation(config, baseline, macro_list, block, population);

  runIslandGA(config, criticality, macro_list, block, population);

  std::vector<PlacementSolution> elites;
  selectParetoElites(config, population, elites);

  if (!refineElitesWithSta(config, criticality, macro_list, block, elites)) {
    logger_->warn(utl::MPL, 115, "TAHPP: STA refinement did not improve elites");
  }

  if (elites.empty()) {
    baseline.applyToBlock(macro_list, block);
    return true;
  }

  std::sort(elites.begin(),
            elites.end(),
            [](const PlacementSolution& a, const PlacementSolution& b) {
              return a.fitness < b.fitness;
            });
  elites.front().applyToBlock(macro_list, block);
  runHybridLocalSearch(config, criticality, macro_list, block, elites.front());

  if (sta_) {
    sta_->updateTiming(false);
    const sta::Slack wns = sta_->worstSlack(sta::MinMax::max());
    logger_->info(utl::MPL,
                  116,
                  "TAHPP complete. Best fitness {:.3f}, WNS {:.3f}",
                  elites.front().fitness,
                  static_cast<double>(wns));
  }

  return true;
}

void TahppOptimizer::seedPopulation(
    const TahppConfig& config,
    const PlacementSolution& baseline,
    const std::vector<odb::dbInst*>& macro_list,
    odb::dbBlock* block,
    std::vector<PlacementSolution>& population)
{
  std::mt19937 rng(config.seed);
  const odb::Rect core = block->getCoreArea();
  const int radius_dbu
      = std::max(1, block->micronsToDbu(config.candidate_radius_um));
  const auto orients = legalOrientations();

  population.clear();
  population.reserve(config.init_population);
  population.push_back(baseline);

  std::uniform_int_distribution<int> macro_dist(
      0, static_cast<int>(macro_list.size()) - 1);
  std::uniform_int_distribution<int> shift_dist(-radius_dbu, radius_dbu);
  std::uniform_int_distribution<int> orient_dist(
      0, static_cast<int>(orients.size()) - 1);
  std::uniform_real_distribution<double> prob(0.0, 1.0);

  while (static_cast<int>(population.size()) < config.init_population) {
    PlacementSolution sol = baseline;
    const int n_perturb = std::max(1, static_cast<int>(sol.macros.size() / 4));
    for (int p = 0; p < n_perturb; ++p) {
      const int midx = macro_dist(rng);
      MacroPlacement& mp = sol.macros[midx];
      mp.x_dbu += shift_dist(rng);
      mp.y_dbu += shift_dist(rng);
      if (prob(rng) < 0.3) {
        mp.orient = static_cast<int>(orients[orient_dist(rng)]);
      }
      odb::dbInst* inst = macro_list[mp.macro_idx];
      const int w = inst->getMaster()->getWidth();
      const int h = inst->getMaster()->getHeight();
      mp.x_dbu = std::clamp(mp.x_dbu, core.xMin(), core.xMax() - w);
      mp.y_dbu = std::clamp(mp.y_dbu, core.yMin(), core.yMax() - h);
    }
    if (solutionLegalOnBlock(placer_, macro_list, sol, block)) {
      population.push_back(sol);
    }
  }

  logger_->info(utl::MPL,
                111,
                "TAHPP: seeded {} legal solutions",
                population.size());
}

void TahppOptimizer::runIslandGA(
    const TahppConfig& config,
    const TimingCriticalityMap& criticality,
    const std::vector<odb::dbInst*>& macro_list,
    odb::dbBlock* block,
    std::vector<PlacementSolution>& population)
{
  GpuBatchEvaluator evaluator(config.use_gpu);
  BatchEvalInput batch_input;
  batch_input.macro_list = &macro_list;
  batch_input.criticality = &criticality;
  batch_input.block_bbox = block->getBBox();
  batch_input.core_bbox = block->getCoreArea();
  batch_input.wns_weight = config.wns_weight;
  batch_input.area_weight = config.area_weight;
  batch_input.wl_weight = config.wl_weight;

  std::mt19937 rng(config.seed + 1);
  const int radius_dbu
      = std::max(1, block->micronsToDbu(config.candidate_radius_um));
  std::uniform_int_distribution<int> shift_dist(-radius_dbu, radius_dbu);
  std::uniform_int_distribution<int> macro_dist(
      0, std::max(0, static_cast<int>(macro_list.size()) - 1));
  std::uniform_real_distribution<double> prob(0.0, 1.0);
  const auto orients = legalOrientations();

  const int island_count = std::max(1, config.num_islands);
  const int per_island = std::max(4, config.pop_size);

  for (int gen = 0; gen < config.num_generations; ++gen) {
    std::vector<PlacementSolution> offspring;
    offspring.reserve(population.size());

    for (int island = 0; island < island_count; ++island) {
      const int start = (island * per_island) % population.size();
      for (int k = 0; k < per_island; ++k) {
        PlacementSolution child = population[(start + k) % population.size()];
        MacroPlacement& mp = child.macros[macro_dist(rng)];
        mp.x_dbu += shift_dist(rng);
        mp.y_dbu += shift_dist(rng);
        if (prob(rng) < 0.25) {
          std::uniform_int_distribution<int> od(
              0, static_cast<int>(orients.size()) - 1);
          mp.orient = static_cast<int>(orients[od(rng)]);
        }
        odb::dbInst* inst = macro_list[mp.macro_idx];
        const odb::Rect core = block->getCoreArea();
        const int w = inst->getMaster()->getWidth();
        const int h = inst->getMaster()->getHeight();
        mp.x_dbu = std::clamp(mp.x_dbu, core.xMin(), core.xMax() - w);
        mp.y_dbu = std::clamp(mp.y_dbu, core.yMin(), core.yMax() - h);
        if (solutionLegalOnBlock(placer_, macro_list, child, block)) {
          offspring.push_back(child);
        }
      }
    }

    population.insert(population.end(), offspring.begin(), offspring.end());

    for (size_t offset = 0; offset < population.size();
         offset += config.gpu_batch_size) {
      const size_t end = std::min(offset + config.gpu_batch_size,
                                  population.size());
      std::vector<PlacementSolution> batch(population.begin() + offset,
                                           population.begin() + end);
      evaluator.evaluateBatch(batch_input, batch);
      for (size_t i = 0; i < batch.size(); ++i) {
        population[offset + i] = batch[i];
      }
    }

    std::sort(population.begin(),
              population.end(),
              [](const PlacementSolution& a, const PlacementSolution& b) {
                return a.fitness < b.fitness;
              });
    if (static_cast<int>(population.size()) > config.init_population) {
      population.resize(config.init_population);
    }

    if (config.migration_interval > 0
        && (gen + 1) % config.migration_interval == 0
        && population.size() >= 2) {
      PlacementSolution migrant = population.front();
      population.back() = migrant;
    }
  }
}

void TahppOptimizer::selectParetoElites(
    const TahppConfig& config,
    std::vector<PlacementSolution>& population,
    std::vector<PlacementSolution>& elites)
{
  elites.clear();
  const int k = std::max(1, config.elite_count);

  std::sort(population.begin(),
            population.end(),
            [](const PlacementSolution& a, const PlacementSolution& b) {
              if (a.wns_proxy != b.wns_proxy) {
                return a.wns_proxy < b.wns_proxy;
              }
              return a.floorplan_area < b.floorplan_area;
            });

  for (const PlacementSolution& sol : population) {
    if (!sol.legal) {
      continue;
    }
    bool dominated = false;
    for (const PlacementSolution& e : elites) {
      if (e.wns_proxy <= sol.wns_proxy && e.floorplan_area <= sol.floorplan_area
          && (e.wns_proxy < sol.wns_proxy
              || e.floorplan_area < sol.floorplan_area)) {
        dominated = true;
        break;
      }
    }
    if (!dominated) {
      elites.push_back(sol);
      elites.erase(
          std::remove_if(elites.begin(),
                         elites.end(),
                         [&sol](const PlacementSolution& e) {
                           return sol.wns_proxy <= e.wns_proxy
                                  && sol.floorplan_area <= e.floorplan_area
                                  && (sol.wns_proxy < e.wns_proxy
                                      || sol.floorplan_area < e.floorplan_area);
                         }),
          elites.end());
    }
    if (static_cast<int>(elites.size()) >= k * 2) {
      break;
    }
  }

  if (elites.size() > static_cast<size_t>(k)) {
    elites.resize(k);
  }
  if (elites.empty() && !population.empty()) {
    elites.push_back(population.front());
  }
}

bool TahppOptimizer::refineElitesWithSta(
    const TahppConfig& config,
    const TimingCriticalityMap& criticality,
    const std::vector<odb::dbInst*>& macro_list,
    odb::dbBlock* block,
    std::vector<PlacementSolution>& elites)
{
  if (!sta_ || elites.empty()) {
    return false;
  }

  const odb::Rect block_bbox = block->getBBox();
  const int radius_dbu
      = std::max(1, block->micronsToDbu(config.candidate_radius_um));
  const double scale = config.timing_weight > 0.0f ? config.timing_weight : 1.0;

  PlacementSolution* best = &elites[0];
  best->applyToBlock(macro_list, block);
  sta_->updateTiming(false);
  sta::Slack best_wns = sta_->worstSlack(sta::MinMax::max());

  std::vector<std::pair<double, int>> macro_order;
  for (size_t i = 0; i < macro_list.size(); ++i) {
    macro_order.emplace_back(criticality.getMacroScore(macro_list[i]),
                             static_cast<int>(i));
  }
  std::sort(macro_order.begin(),
            macro_order.end(),
            [](const auto& a, const auto& b) { return a.first > b.first; });

  const int top_macros
      = std::min(5, static_cast<int>(macro_order.size()));

  for (int iter = 0; iter < config.sta_refine_iterations; ++iter) {
    bool improved = false;
    for (int mi = 0; mi < top_macros; ++mi) {
      const int midx = macro_order[mi].second;
      MacroPlacement& mp = best->macros[midx];
      odb::dbInst* inst = macro_list[midx];
      const int ox = mp.x_dbu;
      const int oy = mp.y_dbu;
      const int oo = mp.orient;

      const std::vector<std::pair<int, int>> moves = {
          {0, 0},
          {radius_dbu, 0},
          {-radius_dbu, 0},
          {0, radius_dbu},
          {0, -radius_dbu}};

      sta::Slack local_best_wns = best_wns;
      int bx = ox;
      int by = oy;

      for (const auto& [dx, dy] : moves) {
        mp.x_dbu = ox + dx;
        mp.y_dbu = oy + dy;
        best->applyToBlock(macro_list, block);
        if (!placer_->findOverlappedMacros(inst).empty()) {
          continue;
        }
        odb::dbMaster* master = inst->getMaster();
        const odb::Rect cand(mp.x_dbu,
                             mp.y_dbu,
                             mp.x_dbu + master->getWidth(),
                             mp.y_dbu + master->getHeight());
        if (!block_bbox.contains(cand)) {
          continue;
        }
        sta_->updateTiming(false);
        const sta::Slack wns = sta_->worstSlack(sta::MinMax::max());
        if (wns > local_best_wns) {
          local_best_wns = wns;
          bx = mp.x_dbu;
          by = mp.y_dbu;
        }
      }

      mp.x_dbu = bx;
      mp.y_dbu = by;
      mp.orient = oo;
      best->applyToBlock(macro_list, block);

      if (local_best_wns > best_wns) {
        best_wns = local_best_wns;
        improved = true;
      }
    }
    if (!improved) {
      break;
    }
  }

  (void) scale;
  return true;
}

void TahppOptimizer::runHybridLocalSearch(
    const TahppConfig& config,
    const TimingCriticalityMap& criticality,
    const std::vector<odb::dbInst*>& macro_list,
    odb::dbBlock* block,
    PlacementSolution& solution)
{
  const odb::Rect block_bbox = block->getBBox();
  int radius_dbu = std::max(1, block->micronsToDbu(config.candidate_radius_um));
  const double scale = config.timing_weight > 0.0 ? config.timing_weight : 1.0;
  const auto orients = legalOrientations();

  solution.applyToBlock(macro_list, block);

  auto globalCost = [&]() {
    double hpwl = 0.0;
    for (odb::dbNet* net : block->getNets()) {
      const auto& iterms = net->getITerms();
      if (iterms.size() < 2) {
        continue;
      }
      int min_x = std::numeric_limits<int>::max();
      int min_y = std::numeric_limits<int>::max();
      int max_x = std::numeric_limits<int>::min();
      int max_y = std::numeric_limits<int>::min();
      for (odb::dbITerm* iterm : iterms) {
        odb::dbInst* inst = iterm->getInst();
        if (!inst) {
          continue;
        }
        const odb::Rect bbox = inst->getBBox()->getBox();
        min_x = std::min(min_x, bbox.xCenter());
        min_y = std::min(min_y, bbox.yCenter());
        max_x = std::max(max_x, bbox.xCenter());
        max_y = std::max(max_y, bbox.yCenter());
      }
      hpwl += static_cast<double>(max_x - min_x + max_y - min_y);
    }
    return hpwl + scale * criticality.computeGlobalAffinityCost(block);
  };

  double current_cost = globalCost();

  std::vector<std::pair<double, odb::dbInst*>> order;
  for (odb::dbInst* inst : macro_list) {
    order.emplace_back(criticality.getMacroScore(inst), inst);
  }
  std::sort(order.begin(), order.end(), [](const auto& a, const auto& b) {
    return a.first > b.first;
  });

  const int active = std::min(static_cast<int>(order.size()),
                              std::max(8, static_cast<int>(order.size() / 2)));

  for (int iter = 0; iter < config.hybrid_iterations; ++iter) {
    bool improved = false;
    for (int i = 0; i < active; ++i) {
      odb::dbInst* inst = order[i].second;
      const auto it = criticality.macro_index.find(inst);
      if (it == criticality.macro_index.end()) {
        continue;
      }
      const int midx = it->second;
      MacroPlacement& mp = solution.macros[midx];
      const int ox = mp.x_dbu;
      const int oy = mp.y_dbu;
      const int oo = mp.orient;

      double best_cost = current_cost;
      int bx = ox;
      int by = oy;
      int bo = oo;

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

      for (const auto& [dx, dy] : moves) {
        for (const auto orient : orients) {
          mp.x_dbu = ox + dx;
          mp.y_dbu = oy + dy;
          mp.orient = static_cast<int>(orient);
          solution.applyToBlock(macro_list, block);
          if (!placer_->findOverlappedMacros(inst).empty()) {
            continue;
          }
          odb::dbMaster* master = inst->getMaster();
          const odb::Rect cand(mp.x_dbu,
                               mp.y_dbu,
                               mp.x_dbu + master->getWidth(),
                               mp.y_dbu + master->getHeight());
          if (!block_bbox.contains(cand)) {
            continue;
          }
          const double cost = globalCost();
          if (cost < best_cost) {
            best_cost = cost;
            bx = mp.x_dbu;
            by = mp.y_dbu;
            bo = mp.orient;
          }
        }
      }

      mp.x_dbu = bx;
      mp.y_dbu = by;
      mp.orient = bo;
      solution.applyToBlock(macro_list, block);

      if (best_cost < current_cost) {
        current_cost = best_cost;
        improved = true;
        radius_dbu = std::min(radius_dbu * 2,
                              block->micronsToDbu(config.candidate_radius_um * 2));
      } else {
        radius_dbu = std::max(1, radius_dbu / 2);
      }
    }
    if (!improved) {
      break;
    }
  }
}

}  // namespace tahpp
}  // namespace mpl
