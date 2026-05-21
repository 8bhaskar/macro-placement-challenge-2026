// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#pragma once

#include <vector>

#include "tahpp/placement_solution.h"
#include "tahpp/timing_criticality.h"

namespace odb {
class dbBlock;
class dbInst;
class Rect;
}  // namespace odb

namespace mpl {
namespace tahpp {

struct BatchEvalInput
{
  const std::vector<odb::dbInst*>* macro_list = nullptr;
  const TimingCriticalityMap* criticality = nullptr;
  odb::Rect block_bbox;
  odb::Rect core_bbox;
  float wns_weight = 0.5f;
  float area_weight = 0.25f;
  float wl_weight = 0.25f;
};

// Evaluates fitness for a batch of placement solutions (CPU or CUDA).
class GpuBatchEvaluator
{
 public:
  explicit GpuBatchEvaluator(bool use_gpu);
  ~GpuBatchEvaluator();

  void evaluateBatch(const BatchEvalInput& input,
                     std::vector<PlacementSolution>& solutions);

  bool gpuAvailable() const { return gpu_available_; }

 private:
  void evaluateBatchCpu(const BatchEvalInput& input,
                        std::vector<PlacementSolution>& solutions);
  void evaluateBatchCuda(const BatchEvalInput& input,
                         std::vector<PlacementSolution>& solutions);

  bool use_gpu_ = false;
  bool gpu_available_ = false;
};

}  // namespace tahpp
}  // namespace mpl
