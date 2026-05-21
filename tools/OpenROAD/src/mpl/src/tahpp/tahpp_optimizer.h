// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#pragma once

#include "mpl/rtl_mp.h"
#include "tahpp/placement_solution.h"
#include "tahpp/timing_criticality.h"

namespace mpl {
namespace tahpp {

class TahppOptimizer
{
 public:
  TahppOptimizer(MacroPlacer* placer,
                 odb::dbDatabase* db,
                 sta::dbSta* sta,
                 utl::Logger* logger);

  bool run(const TahppConfig& config,
           const std::vector<odb::dbInst*>& macro_list,
           odb::dbBlock* block,
           TimingCriticalityMap& criticality);

 private:
  void seedPopulation(const TahppConfig& config,
                      const PlacementSolution& baseline,
                      const std::vector<odb::dbInst*>& macro_list,
                      odb::dbBlock* block,
                      std::vector<PlacementSolution>& population);

  void runIslandGA(const TahppConfig& config,
                   const TimingCriticalityMap& criticality,
                   const std::vector<odb::dbInst*>& macro_list,
                   odb::dbBlock* block,
                   std::vector<PlacementSolution>& population);

  void selectParetoElites(const TahppConfig& config,
                          std::vector<PlacementSolution>& population,
                          std::vector<PlacementSolution>& elites);

  bool refineElitesWithSta(const TahppConfig& config,
                           const TimingCriticalityMap& criticality,
                           const std::vector<odb::dbInst*>& macro_list,
                           odb::dbBlock* block,
                           std::vector<PlacementSolution>& elites);

  void runHybridLocalSearch(const TahppConfig& config,
                            const TimingCriticalityMap& criticality,
                            const std::vector<odb::dbInst*>& macro_list,
                            odb::dbBlock* block,
                            PlacementSolution& solution);

  MacroPlacer* placer_ = nullptr;
  odb::dbDatabase* db_ = nullptr;
  sta::dbSta* sta_ = nullptr;
  utl::Logger* logger_ = nullptr;
};

}  // namespace tahpp
}  // namespace mpl
