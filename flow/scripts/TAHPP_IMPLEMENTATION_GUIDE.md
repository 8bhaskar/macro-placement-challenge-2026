# TAHPP Implementation Guide

## Table of Contents
1. Data Structures
2. Phase 1: Timing Analysis Implementation
3. Phase 2: GPU Initialization Implementation
4. Phase 3: Island GA Implementation
5. Phase 4: Slack-Flow Refinement Implementation
6. Tcl Integration & Command Registration

---

## 1. Core Data Structures

### 1.1 Solution Representation

```cpp
// Unified solution representation
struct PlacementSolution {
  std::vector<int> macro_ids;           // ordering of macros
  std::map<int, Point> positions;       // macro_id -> (x, y)
  std::map<int, Orientation> orientations;  // macro_id -> rotation
  
  // Fitness components
  float hpwl = 0.0f;                    // half-perimeter wirelength
  float area_penalty = 0.0f;
  float timing_penalty = 0.0f;          // WNS-related
  float outline_penalty = 0.0f;
  float combined_fitness = 0.0f;
  
  // Metadata
  uint64_t generation_created = 0;
  int island_id = -1;                   // which island created this
  bool is_legal = false;
  
  // For GPU transfer
  void* gpu_handle = nullptr;           // CUDA device memory reference
};

struct PopulationState {
  std::vector<PlacementSolution> solutions;
  std::vector<float> fitness_scores;
  int generation = 0;
  int island_id = -1;
  
  // For tracking diversity and convergence
  float avg_fitness = 0.0f;
  float max_fitness = 0.0f;
  float min_fitness = 1e9f;
  float fitness_std_dev = 0.0f;
};

struct ParetoFront {
  std::vector<PlacementSolution> solutions;
  std::vector<std::pair<float,float>> objectives;  // (wns, area) pairs
  
  // Check if solution dominates existing
  bool isDominated(const PlacementSolution& sol, 
                   float wns, float area);
  
  // Add if non-dominated
  void addIfNonDominated(const PlacementSolution& sol);
};
```

### 1.2 Timing Criticality

```cpp
struct MacroCriticality {
  int macro_id;
  float criticality_score;              // 0.0 to 1.0
  CriticalityLevel level;               // CRITICAL, MODERATE, NON_CRITICAL
  
  // Timing connections
  std::vector<int> critical_net_ids;    // nets on critical paths
  std::vector<int> sink_macro_ids;      // macros sinks of critical nets
  
  float slack_contribution;             // how much this macro impacts WNS
};

struct TimingAffinity {
  int macro_a, macro_b;
  float affinity_score;                 // -1.0 (push apart) to +1.0 (pull close)
  float preferred_distance;             // optimal distance in microns
};

class TimingAnalyzer {
  std::vector<MacroCriticality> macro_criticalities;
  std::vector<TimingAffinity> affinity_matrix;
  
  // Main entry
  void analyze(Database* db, Sta* sta) {
    computeCriticalities();
    extractTimingPaths();
    computeAffinities();
  }
  
  void computeCriticalities() {
    // For each macro:
    //   1. Get connected nets
    //   2. Query slack on each net from STA
    //   3. Compute: criticality = sum(abs(min_slack)) / max_path_delay
    //   4. Normalize and classify
  }
  
  void computeAffinities() {
    // For each critical path:
    //   1. Extract macro chain (source -> intermediate -> sink)
    //   2. Set affinity = POSITIVE for consecutive macros
    //   3. Affinity distance = typical macro dimension
  }
  
  float getCriticalityScore(int macro_id);
  bool isCritical(int macro_id);
  float getAffinityScore(int macro_a, int macro_b);
};
```

### 1.3 GPU Data Buffer

```cpp
struct GPUBatchBuffer {
  // For batched evaluation of 64 solutions
  
  // Host-side
  std::vector<PlacementSolution> host_solutions;
  std::vector<float> host_fitness_scores;
  
  // Device-side (CUDA)
  int* d_macro_positions_x;             // [batch_size * max_macros]
  int* d_macro_positions_y;
  float* d_fitness_output;              // [batch_size]
  
  int batch_size = 64;
  int max_macros = 512;  // Ariane133 worst case
  
  void uploadBatch(const std::vector<PlacementSolution>& solutions);
  void computeHPWLBatch();              // GPU kernel
  void computeLegalizationBatch();      // GPU kernel
  void downloadResults();
};
```

---

## 2. Phase 1: Timing Analysis Implementation

### Algorithm Pseudocode

```
function analyzeTimingAndExtractCriticality():
  // Input: placed_db (before macro placement)
  // Output: criticality_scores, timing_affinity_matrix
  
  sta = openSTA(db)
  sta.runPowerFlow()
  
  // Extract macro timing impact
  for each unfixed_macro in design:
    macro.criticality_score = 0.0
    
    for each net connected_to macro:
      slack_on_net = sta.getSlack(net)
      path_delay = sta.getPathDelay(net)
      
      // Criticality = how much this connection affects WNS
      contribution = max(0, -slack_on_net) / (path_delay + 1e-6)
      macro.criticality_score += contribution
    
    // Normalize
    macro.criticality_score /= (1 + macro.getNetCount())
    
    // Classify
    if criticality_score > CRITICAL_THRESHOLD:
      macro.level = CRITICAL
    else if criticality_score > MODERATE_THRESHOLD:
      macro.level = MODERATE
    else:
      macro.level = NON_CRITICAL
  
  // Extract critical paths for affinity
  critical_paths = sta.getCriticalPaths(num_paths=50)
  
  for each path in critical_paths:
    macro_sequence = extractMacroSequence(path)
    
    // Create positive affinity between consecutive macros
    for i = 0 to macro_sequence.length - 2:
      macro_a = macro_sequence[i]
      macro_b = macro_sequence[i+1]
      affinity = computePositiveAffinity(macro_a, macro_b)
      affinity_matrix[macro_a][macro_b] = affinity
      affinity_matrix[macro_b][macro_a] = affinity
  
  return criticality_scores, affinity_matrix
```

### CPP Implementation Sketch

```cpp
class TimingPhaseHandler {
 private:
  TimingAnalyzer* timing_;
  std::map<int, MacroCriticality> criticalities_;
  std::vector<std::vector<float>> affinity_matrix_;
  utl::Logger* logger_;
  
 public:
  void run(odb::dbDatabase* db, sta::dbSta* sta) {
    logger_->info(MPL, 100, "Phase 1: Timing Analysis Started");
    
    // STA computation (2 min runtime on Ariane133)
    timing_->computeCriticalities(db, sta);
    timing_->extractTimingPaths(50);  // Top 50 critical paths
    timing_->computeAffinities();
    
    // Extract into local maps
    for (const auto& [macro_id, crit] : timing_->getCriticalities()) {
      criticalities_[macro_id] = crit;
    }
    affinity_matrix_ = timing_->getAffinityMatrix();
    
    logger_->info(MPL, 101, "Phase 1: Found {} critical macros",
                  countCriticalMacros());
    
    // Save to file for later phases
    saveTimingData("rtlmp_timing_analysis.bin");
  }
  
  float getCriticalityScore(int macro_id) const {
    auto it = criticalities_.find(macro_id);
    return (it != criticalities_.end()) ? it->second.criticality_score : 0.0f;
  }
  
  bool isCritical(int macro_id) const {
    return criticalities_.at(macro_id).level == CRITICAL;
  }
};
```

---

## 3. Phase 2: GPU Initialization Implementation

### GPU CUDA Kernel (Legalization)

```cuda
// CUDA kernel: parallel sequence packing (place 16 macros per warp)
__global__ void parallelSequencePack(
    int* d_macro_x, int* d_macro_y,  // [batch_size * max_macros]
    int* d_widths, int* d_heights,
    float* d_fitness,
    int batch_size, int num_macros, int die_width, int die_height
) {
  int thread_idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (thread_idx >= batch_size) return;
  
  int batch_id = thread_idx;
  
  // Each thread handles one solution
  // Place macros left-to-right, bottom-to-top
  int current_x = 0, current_y = 0, max_y_in_row = 0;
  
  for (int i = 0; i < num_macros; i++) {
    int w = d_widths[batch_id * num_macros + i];
    int h = d_heights[batch_id * num_macros + i];
    
    // Check if macro fits in current row
    if (current_x + w > die_width) {
      // Move to next row
      current_x = 0;
      current_y += max_y_in_row;
      max_y_in_row = 0;
    }
    
    // Check if fits in die
    if (current_y + h > die_height) {
      // Illegal placement; mark as invalid
      d_fitness[batch_id] = 1e9f;
      return;
    }
    
    // Place macro
    d_macro_x[batch_id * num_macros + i] = current_x;
    d_macro_y[batch_id * num_macros + i] = current_y;
    
    current_x += w;
    max_y_in_row = max(max_y_in_row, h);
  }
  
  // Compute HPWL for this solution
  __shared__ float local_hpwl[256];
  local_hpwl[threadIdx.x] = computeHPWLForBatch(
      d_macro_x, d_macro_y, batch_id, num_macros);
  
  __syncthreads();
  
  // Store fitness
  d_fitness[batch_id] = local_hpwl[threadIdx.x];
}
```

### Host-Side Initialization

```cpp
class InitializationPhase {
 private:
  GPUBatchBuffer* gpu_buffer_;
  std::vector<PlacementSolution> initial_population_;
  
 public:
  void generateInitialPopulation(
      int population_size,
      const TimingPhaseHandler& timing) {
    
    logger_->info(MPL, 102, "Phase 2: Generating {} initial solutions",
                  population_size);
    
    // Generate random orderings
    std::vector<std::vector<int>> macro_orderings;
    for (int i = 0; i < population_size; i++) {
      auto ordering = generateRandomOrdering(db_->getUnfixedMacros());
      macro_orderings.push_back(ordering);
    }
    
    // Batch process on GPU (512 at a time)
    int num_batches = (population_size + 63) / 64;
    
    for (int batch = 0; batch < num_batches; batch++) {
      int batch_start = batch * 64;
      int batch_size = std::min(64, population_size - batch_start);
      
      // Upload to GPU
      gpu_buffer_->uploadBatch(
          std::vector<PlacementSolution>(
              macro_orderings.begin() + batch_start,
              macro_orderings.begin() + batch_start + batch_size));
      
      // Legalize on GPU
      gpu_buffer_->computeLegalizationBatch();
      
      // Download results
      auto batch_results = gpu_buffer_->downloadResults();
      
      // For each result, compute timing fitness on CPU
      for (const auto& sol : batch_results) {
        float timing_score = computeTimingFitness(sol, timing);
        sol.timing_penalty = timing_score;
        sol.combined_fitness = 
            0.5 * sol.hpwl + 0.3 * sol.area_penalty + 0.2 * timing_score;
        initial_population_.push_back(sol);
      }
      
      // Progress indicator
      if ((batch + 1) % 4 == 0) {
        logger_->info(MPL, 103, "  Generated {} / {} solutions",
                      (batch + 1) * 64, population_size);
      }
    }
    
    // Sort and keep top 128
    std::sort(initial_population_.begin(), initial_population_.end(),
              [](const auto& a, const auto& b) {
                return a.combined_fitness < b.combined_fitness;
              });
    initial_population_.resize(128);
    
    logger_->info(MPL, 104, "Phase 2: Top 128 solutions selected");
  }
  
 private:
  std::vector<int> generateRandomOrdering(
      const std::vector<odb::dbInst*>& macros) {
    std::vector<int> ordering;
    for (const auto* macro : macros) {
      ordering.push_back(macro->getId());
    }
    std::random_shuffle(ordering.begin(), ordering.end());
    return ordering;
  }
};
```

---

## 4. Phase 3: Island GA Implementation

### Multi-Island Orchestrator

```cpp
class IslandGAOptimizer {
 private:
  static const int NUM_ISLANDS = 16;
  std::vector<PopulationState> islands;
  std::vector<std::thread> worker_threads;
  ParetoFront pareto_front;
  GPUBatchBuffer* gpu_buffer_;
  
 public:
  void optimize(
      const std::vector<PlacementSolution>& initial_pop,
      const TimingPhaseHandler& timing) {
    
    logger_->info(MPL, 105, "Phase 3: Starting Island GA");
    
    // Initialize islands with seeded diversity
    for (int i = 0; i < NUM_ISLANDS; i++) {
      islands[i].island_id = i;
      islands[i].solutions = initial_pop;  // All start from same top-128
      
      // Diversify by random perturbations
      for (auto& sol : islands[i].solutions) {
        perturb(sol, 5 + i);  // Different random seeds per island
      }
    }
    
    // Main GA loop
    for (int gen = 0; gen < 40; gen++) {
      logger_->info(MPL, 106, "GA Generation {}/40", gen + 1);
      
      // Evaluate all islands in parallel
      #pragma omp parallel for num_threads(NUM_ISLANDS)
      for (int island_id = 0; island_id < NUM_ISLANDS; island_id++) {
        evaluateIsland(island_id, timing, gen);
      }
      
      // Migration every 10 generations
      if (gen > 0 && gen % 10 == 0) {
        performMigration();
      }
      
      // Update global Pareto front
      for (int i = 0; i < NUM_ISLANDS; i++) {
        for (const auto& sol : islands[i].solutions) {
          if (sol.is_legal) {
            pareto_front.addIfNonDominated(sol);
          }
        }
      }
    }
    
    logger_->info(MPL, 107, "Phase 3: Complete. Pareto front size: {}",
                  pareto_front.solutions.size());
  }
  
 private:
  void evaluateIsland(int island_id,
                      const TimingPhaseHandler& timing,
                      int generation) {
    auto& island = islands[island_id];
    auto& solutions = island.solutions;
    
    // Selection: tournament
    std::vector<PlacementSolution> parents;
    for (int i = 0; i < 32; i++) {
      auto winner = tournamentSelect(solutions, 3);
      parents.push_back(winner);
    }
    
    // Mutation & crossover
    std::vector<PlacementSolution> offspring;
    for (int i = 0; i < 32; i++) {
      PlacementSolution child = parents[i];
      
      // ML-guided mutation
      float rand = distribution_(generator_);
      if (rand < 0.60) {
        // Critical swap
        mutateCriticalSwap(child, timing, island_id);
      } else if (rand < 0.85) {
        // Position targeting
        mutatePositionTarget(child, timing, island_id);
      } else if (rand < 0.95) {
        // Rotation
        mutateRotate(child);
      } else {
        // Local search
        mutateLocalSearch(child);
      }
      
      offspring.push_back(child);
    }
    
    // Crossover
    for (int i = 0; i < 16; i++) {
      auto [p1, p2] = tournamentSelectPair(solutions, 2);
      auto child = crossover(p1, p2);
      offspring.push_back(child);
    }
    
    // Batch evaluate offspring on GPU
    auto batch_fitness = gpuBatchEvaluate(offspring);
    for (size_t i = 0; i < offspring.size(); i++) {
      offspring[i].hpwl = batch_fitness[i].hpwl;
      offspring[i].area_penalty = batch_fitness[i].area;
      offspring[i].is_legal = batch_fitness[i].is_legal;
    }
    
    // CPU timing evaluation on critical subset
    for (auto& sol : offspring) {
      if (sol.is_legal) {
        float timing_score = computeTimingFitness(sol, timing);
        sol.timing_penalty = timing_score;
      }
    }
    
    // Replace worst in population
    std::sort(offspring.begin(), offspring.end(),
              [](const auto& a, const auto& b) {
                return a.combined_fitness < b.combined_fitness;
              });
    
    solutions.erase(solutions.begin() + solutions.size() - offspring.size(),
                    solutions.end());
    solutions.insert(solutions.end(), offspring.begin(), offspring.end());
    
    // Recompute stats
    updateIslandStats(island);
  }
  
  void performMigration() {
    // Ring topology migration
    for (int i = 0; i < NUM_ISLANDS; i++) {
      int next = (i + 1) % NUM_ISLANDS;
      
      // Send top 4 solutions to next island
      auto top4 = getTopSolutions(islands[i].solutions, 4);
      
      // Replace bottom 4 in next island
      std::sort(islands[next].solutions.begin(),
                islands[next].solutions.end(),
                [](const auto& a, const auto& b) {
                  return a.combined_fitness > b.combined_fitness;
                });
      islands[next].solutions.erase(
          islands[next].solutions.end() - 4,
          islands[next].solutions.end());
      islands[next].solutions.insert(
          islands[next].solutions.end(),
          top4.begin(), top4.end());
    }
  }
};
```

---

## 5. Phase 4: Slack-Flow Refinement

### Slack-Flow Local Search

```cpp
class SlackFlowRefinement {
 private:
  sta::dbSta* sta_;
  utl::Logger* logger_;
  
 public:
  PlacementSolution refineTopSolutions(
      const ParetoFront& pareto_front,
      odb::dbDatabase* db) {
    
    logger_->info(MPL, 108, "Phase 4: Slack-Flow Refinement");
    
    PlacementSolution best_solution;
    float best_wns = -1e9;
    
    // Refine each of top 4 solutions
    for (size_t idx = 0; idx < std::min(size_t(4), pareto_front.solutions.size()); idx++) {
      auto candidate = pareto_front.solutions[idx];
      
      logger_->info(MPL, 109, "Refining solution {}/4", idx + 1);
      
      for (int iter = 0; iter < 100; iter++) {
        // Run STA
        sta_->run();
        float current_wns = sta_->worstNegativeSlack();
        
        // Extract timing critical macros (top 5)
        auto critical_macros = extractCriticalMacros(5);
        
        // Compute slack gradient for each
        for (const auto& macro_id : critical_macros) {
          Point gradient = computeSlackGradient(macro_id);
          
          // Move macro in direction of gradient
          float step_size = 0.1;  // microns
          candidate.positions[macro_id].x += step_size * gradient.x;
          candidate.positions[macro_id].y += step_size * gradient.y;
        }
        
        // Legalize
        legalizeOnGPU(candidate);
        
        // Check if improved
        sta_->run();
        float new_wns = sta_->worstNegativeSlack();
        
        if (new_wns > current_wns) {
          logger_->info(MPL, 110, "  Iter {}: WNS improved to {}", iter, new_wns);
        } else {
          // Revert
          for (const auto& macro_id : critical_macros) {
            candidate.positions[macro_id] = getPreviousPosition(macro_id);
          }
        }
      }
      
      // Track best
      sta_->run();
      float final_wns = sta_->worstNegativeSlack();
      if (final_wns > best_wns) {
        best_wns = final_wns;
        best_solution = candidate;
      }
    }
    
    logger_->info(MPL, 111, "Phase 4: Best WNS achieved: {} ns", best_wns);
    return best_solution;
  }
  
 private:
  std::vector<int> extractCriticalMacros(int top_k) {
    // Extract macros that are on critical paths
    std::vector<int> critical_macros;
    auto critical_paths = sta_->getCriticalPaths(top_k);
    
    for (const auto& path : critical_paths) {
      for (const auto& arc : path.getArcs()) {
        int macro_id = arcToMacroId(arc);
        if (macro_id != -1) {
          critical_macros.push_back(macro_id);
        }
      }
    }
    
    // Deduplicate and return
    std::sort(critical_macros.begin(), critical_macros.end());
    critical_macros.erase(
        std::unique(critical_macros.begin(), critical_macros.end()),
        critical_macros.end());
    
    return critical_macros;
  }
  
  Point computeSlackGradient(int macro_id) {
    // Estimate delay change if macro moves by (dx, dy)
    // gradient = -∇(slack) = ∇(delay)
    
    float dx_sensitivity = 0.0f, dy_sensitivity = 0.0f;
    
    // Sample: move macro by small amount, recompute slack
    Point original = getCurrentPosition(macro_id);
    
    for (float delta : {-1.0f, 1.0f}) {  // +/- 1 micron
      moveToPosition(macro_id, original.x + delta, original.y);
      sta_->run();
      float slack_at_delta = sta_->worstNegativeSlack();
      dx_sensitivity += (slack_at_delta - original_slack) * delta;
      
      moveToPosition(macro_id, original.x, original.y + delta);
      sta_->run();
      float slack_at_delta_y = sta_->worstNegativeSlack();
      dy_sensitivity += (slack_at_delta_y - original_slack) * delta;
    }
    
    // Restore
    moveToPosition(macro_id, original.x, original.y);
    
    // Gradient points toward improving slack
    return Point(dx_sensitivity, dy_sensitivity).normalize();
  }
};
```

---

## 6. Tcl Integration & Command Registration

### Tcl Command Registration

```cpp
// In mpl/src/object.cpp or similar

extern "C" {
  
  // Callback for rtl_macro_placer_hybrid command
  int rtlMacroPlacerHybrid_cmd(ClientData clientData,
                               Tcl_Interp* interp,
                               int argc,
                               const char* argv[]) {
    // Parse arguments
    int num_threads = 16;
    bool use_gpu = true;
    float weight_timing = 0.5;
    float weight_area = 0.25;
    float weight_wirelength = 0.25;
    
    // (Parse TCL args...)
    
    // Run algorithm
    TAHPPOptimizer optimizer(
        ord::get_db(),
        ord::get_sta(),
        num_threads,
        use_gpu);
    
    optimizer.optimize(weight_timing, weight_area, weight_wirelength);
    
    return TCL_OK;
  }
  
  // Register command
  int Mpl_Hybrid_Init(Tcl_Interp* interp) {
    Tcl_CreateCommand(interp, "rtl_macro_placer_hybrid",
                      rtlMacroPlacerHybrid_cmd,
                      nullptr, nullptr);
    return TCL_OK;
  }
}
```

### Tcl Flow Integration

```tcl
# In macro_place_util.tcl, add variant:

if { [env_var_exists_and_non_empty USE_HYBRID_PLACEMENT] && $::env(USE_HYBRID_PLACEMENT) == 1 } {
  log_cmd rtl_macro_placer_hybrid \
    -num_threads $num_threads \
    -use_gpu 1 \
    -weight_timing $::env(HYBRID_WEIGHT_TIMING) \
    -weight_area $::env(HYBRID_WEIGHT_AREA) \
    -weight_wirelength $::env(HYBRID_WEIGHT_WIRELENGTH) \
    -max_iterations 40
} else {
  log_cmd rtl_macro_placer {*}$all_args
}
```

---

## Performance Targets & Profiling

### Expected Runtime Breakdown (20 minutes total)

| Phase | Duration | Parallelism |
|-------|----------|-------------|
| Phase 1: Timing Analysis | 2 min | CPU 2 cores |
| Phase 2: GPU Initialization | 3 min | GPU + CPU 2 cores |
| Phase 3: Island GA | 8 min | CPU 16 cores + GPU |
| Phase 4: Slack-Flow Refinement | 5 min | CPU 8 cores + GPU |
| Phase 5: Post-processing | 2 min | CPU 2 cores |

### GPU Profiling Points

```cuda
// In GPU kernels, measure:
cudaEvent_t start, stop;
cudaEventCreate(&start);
cudaEventCreate(&stop);

cudaEventRecord(start);
// ... kernel launch ...
cudaEventRecord(stop);
cudaEventSynchronize(stop);
float ms;
cudaEventElapsedTime(&ms, start, stop);

logger_->info(MPL, 200, "GPU kernel runtime: {} ms", ms);
```

### CPU Thread Analysis

```cpp
// Use OpenMP profiling
#pragma omp parallel for schedule(dynamic) num_threads(16)
for (int i = 0; i < NUM_ISLANDS; i++) {
  // Each thread tracks its own work
  auto t_start = std::chrono::high_resolution_clock::now();
  
  evaluateIsland(i, ...);
  
  auto t_end = std::chrono::high_resolution_clock::now();
  auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);
  
  logger_->info(MPL, 201, "Island {} evaluation: {} ms", i, duration_ms.count());
}
```

---

## Validation Checklist

- [ ] Phase 1: Criticality scores match STA slack ordering
- [ ] Phase 2: Initial population diversity > 50% unique solutions
- [ ] Phase 3: Pareto front grows monotonically
- [ ] Phase 3: Island migration reduces local optima traps
- [ ] Phase 4: Slack-flow iterations improve WNS by > 1% per iteration initially
- [ ] Phase 5: Final solution is legal (no overlaps, within core area)
- [ ] Runtime: < 25 minutes on Ariane133/nangate45

---

## Testing on Ariane133

```bash
# Enable hybrid placer in flow
export USE_HYBRID_PLACEMENT=1
export HYBRID_WEIGHT_TIMING=0.50
export HYBRID_WEIGHT_AREA=0.25
export HYBRID_WEIGHT_WIRELENGTH=0.25

# Run flow
cd flow
make run_rtl DESIGN=ariane133 PLATFORM=nangate45 SYNTH_HIERARCHICAL=0

# Compare results
diff <(grep "TNS\|WNS" results/nangate45/ariane133/base/6_report.txt) \
     <(grep "TNS\|WNS" reference_results.txt)
```

