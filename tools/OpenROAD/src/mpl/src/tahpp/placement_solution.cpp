// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#include "tahpp/placement_solution.h"

#include "odb/db.h"

namespace mpl {
namespace tahpp {

void PlacementSolution::captureFromBlock(
    const std::vector<odb::dbInst*>& macro_list,
    odb::dbBlock* block)
{
  (void) block;
  macros.clear();
  macros.reserve(macro_list.size());
  for (size_t i = 0; i < macro_list.size(); ++i) {
    odb::dbInst* inst = macro_list[i];
    const odb::Rect bbox = inst->getBBox()->getBox();
    MacroPlacement mp;
    mp.macro_idx = static_cast<int>(i);
    mp.x_dbu = bbox.xMin();
    mp.y_dbu = bbox.yMin();
    mp.orient = static_cast<int>(inst->getOrient().getValue());
    macros.push_back(mp);
  }
}

void PlacementSolution::applyToBlock(
    const std::vector<odb::dbInst*>& macro_list,
    odb::dbBlock* block) const
{
  (void) block;
  for (const MacroPlacement& mp : macros) {
    if (mp.macro_idx < 0
        || mp.macro_idx >= static_cast<int>(macro_list.size())) {
      continue;
    }
    odb::dbInst* inst = macro_list[mp.macro_idx];
    inst->setOrient(odb::dbOrientType(static_cast<odb::dbOrientType::Value>(mp.orient)));
    inst->setLocation(mp.x_dbu, mp.y_dbu);
    inst->setPlacementStatus(odb::dbPlacementStatus::PLACED);
  }
}

}  // namespace tahpp
}  // namespace mpl
