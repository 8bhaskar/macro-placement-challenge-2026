# Advanced Macro Placement Algorithm for Ariane133/Nangate45

## Executive Summary

This document proposes a **Timing-Aware Hybrid Placement Algorithm** that combines GPU-accelerated population-based search with CPU-level timing analysis and machine learning-driven macro arrangement. This approach targets significant improvements in WNS (worst negative slack) and area utilization compared to the current hierarchical SA-based placer.

---

## Current Algorithm Analysis

### Strengths
- Hierarchical clustering reduces problem complexity
- Simulated annealing provides reasonable solution quality
- Fast execution on CPUs

### Weaknesses
1. **No timing awareness**: Uses generic wirelength minimization (HPWL), not slack-driven placement
2. **Sequential optimization**: Only one solution explored per temperature step
3. **Limited utilization of CPU cores**: Current SA is single-threaded
4. **No GPU utilization**: Leaves accelerators unused
5. **Static weights**: Area/wirelength/outline weights don't adapt to design needs
6. **Poor macro criticality handling**: Treats all macros equally regardless of timing impact

---

## Proposed Algorithm: Timing-Aware Hybrid Population Placer (TAHPP)

### High-Level Overview

```
┌─────────────────────────────────────────────────────────────────┐
│ Phase 1: Pre-processing & Analysis (CPU 2 cores, ~2 min)        │
│ - Extract timing paths, critical macros, I/O constraints         │
│ - Build timing graph with STA slack computation                  │
│ - Compute timing criticality map for each macro                  │
└─────────────────────────────────────────────────────────────────┘
              ↓
┌─────────────────────────────────────────────────────────────────┐
│ Phase 2: Initial Placement Population (GPU + CPU, ~3 min)        │
│ - Generate 512-1024 initial placements in parallel               │
│ - Use GPU for fast legalization & constraint checking            │
│ - Compute fitness = timing + area + wirelength                   │
│ - Rank solutions, keep top 128                                   │
└─────────────────────────────────────────────────────────────────┘
              ↓
┌─────────────────────────────────────────────────────────────────┐
│ Phase 3: ML-Guided Multi-Population Search (CPU + GPU, ~8 min)   │
│ - 16 independent populations (one per core)                      │
│ - Population size: 64 solutions each = 1024 total diversity      │
│ - GPU evaluates populations in parallel (batches of 64)          │
│ - ML-driven mutation: bias moves toward timing improvement       │
│ - Crossover: recombine high-fitness solutions                    │
│ - Island model: periodic migration between populations           │
└─────────────────────────────────────────────────────────────────┘
              ↓
┌─────────────────────────────────────────────────────────────────┐
│ Phase 4: Timing Refinement (CPU 8 cores, ~5 min)                 │
│ - Extract top 4 non-dominated solutions by Pareto front          │
│ - Local search: targeted moves on timing-critical macros         │
│ - Constrain moves based on slack flow analysis                   │
│ - Final legalization and overlap resolution                      │
└─────────────────────────────────────────────────────────────────┘
              ↓
┌─────────────────────────────────────────────────────────────────┐
│ Phase 5: Post-processing (CPU 2 cores, ~2 min)                   │
│ - Halo refinement for congestion mitigation                      │
│ - I/O optimization: pull I/O-heavy macros toward pads            │
│ - Final area computation and report generation                   │
└─────────────────────────────────────────────────────────────────┘
```

**Total estimated runtime: ~20 minutes**

---

## Detailed Algorithmic Descriptions

### Phase 1: Timing Analysis & Criticality Extraction

**Input**: Placed netlist with SDC constraints

**Algorithm**:
1. Run full STA to compute slack at all timing points
2. For each macro instance:
   - Extract all nets connected to macro pins
   - Compute **timing criticality score** = $\frac{\text{max slack in connected nets}}{\text{longest path delay}}$
   - Classify as: critical (top 15%), moderate (15-50%), non-critical (bottom 50%)
3. Identify **critical timing paths** crossing macro regions
4. Build **distance matrix** for timing-critical macro pairs: 
   - Distance = distance penalty if placed far apart
   - Reward = negative penalty if placed close

**Output**:
- `criticality_scores.json`: per-macro timing criticality
- `critical_paths.txt`: critical net lists
- `distance_matrix.bin`: timing affinity matrix (sparse format)

**Why this helps WNS**:
- Later phases can bias placement to keep critical macros close
- Timing-aware moves are prioritized over generic wirelength moves

---

### Phase 2: GPU-Accelerated Population Initialization

**Input**: Timing criticality map, design constraints

**Algorithm (GPU Kernel)**:
```cuda
// Parallel initialization of K initial solutions
for each solution in parallel (GPU threads):
  1. Random macro ordering (seed-based)
  2. Fast legalize using GPU sequence packing:
     - Place macros left-to-right, bottom-to-top
     - GPU threads handle 16 macros in parallel
  3. Compute fitness score = w_timing*T + w_area*A + w_wirelength*WL:
     - GPU: Parallel HPWL computation (reduction tree)
     - CPU: Quick STA on critical path subset (~50ms per solution)
  4. Store solution + fitness in unified memory
```

**CPU coordination**:
- Launch 512 GPU threads (one per solution)
- For each batch of 64 solutions:
  - GPU computes HPWL + legalization
  - CPU runs fast critical path STA
  - Synchronize fitness values
- Rank solutions, select top 128 for Phase 3

**Output**:
- Initial population of 128 high-quality placements
- Each solution: macro positions + fitness score

**Why this helps**:
- GPU parallelization gives 100x speedup on legalization
- Diverse initial population reduces local optima risk

---

### Phase 3: ML-Guided Multi-Population Search

**Architecture**:
```
16 island populations (1 per CPU core)
Each island runs independent GA variant:
  - Population: 64 solutions
  - Generations: 40
  - Migration: every 10 gen to adjacent islands
```

**Per-Island Algorithm (CPU core)**:
```
For each generation:
  1. Evaluate all 64 solutions (GPU batch):
     - Send batch to GPU for HPWL + legalization
     - Receive fitness scores back
  
  2. Selection (tournament, size=3):
     - Higher likelihood for high-fitness solutions
     - Stochastic selection to maintain diversity
  
  3. Mutation with ML guidance:
     For each of 32 offspring:
       a. Select parent from tournament winners
       b. Choose mutation type:
          - Type A (60%): Swap two critical macros (ML-learned)
          - Type B (25%): Move macro toward target position (timing-aware)
          - Type C (10%): Rotate macro for better orientation
          - Type D (5%): Local search refinement
       c. Apply GPU-accelerated legalization
  
  4. Crossover (30% of offspring):
     - Select two parents from elite solutions
     - Position crossover: split macro list at random point
     - Recombine and legalize via GPU
  
  5. Migration (every 10 generations):
     - Send 4 best solutions to adjacent islands
     - Receive 4 solutions from neighbors
     - Replace 4 worst solutions with received ones
```

**ML-Driven Mutation Details**:
- **Decision tree model**: trained on which macro swaps improve WNS
  - Features: macro size, timing criticality, current position, target position
  - Predicts: probability each swap reduces worst slack
- **Position targeting**: for critical macros, compute "optimal region"
  - Region = area minimizing delay to critical sinks
  - Use timing-driven force model: F = -∇(slack)
  - Mutation biases macros toward high-slack regions

**GPU Batch Evaluation**:
- Collect 64 candidate solutions
- Upload to GPU unified memory
- Parallel HPWL: each macro pair distance computed in parallel
- Parallel legalization: 16 solutions legalized simultaneously
- Download fitness scores

**Output**:
- 16 × 64 = 1024 placements explored per generation
- 40 generations × 1024 = 40,960 solutions evaluated
- Best non-dominated solutions tracked (Pareto front)

**Why this helps WNS and Area**:
- Population diversity explores more solution space than serial SA
- ML guidance biases search toward timing improvements
- Multi-population island model prevents premature convergence
- GPU batching enables evaluation of massive search space

---

### Phase 4: Timing Refinement via Slack-Flow Placement

**Input**: Top 4 non-dominated solutions from Phase 3

**Algorithm**:
```
For each of 4 candidate solutions:
  1. Run detailed STA
     - Compute slack at all timing points
     - Identify worst slack path
  
  2. Extract timing-critical macro pairs:
     - For each net on critical path, find source/sink macros
     - Build "critical macro chain"
  
  3. Slack-flow directed search (CPU, ~2 min per solution):
     For i = 1 to 100 iterations:
       a. Compute timing gradient: ∇slack for each macro position
       b. Move top-5 critical macros in direction of ∇slack
       c. Legalize with minimal disruption (GPU, ~1 sec)
       d. Re-run STA on critical path subset
       e. If WNS improves: accept; else: revert
  
  4. Local detail placement optimization
     - For each macro with potential congestion conflicts:
       - Try small translations and rotations
       - Keep if no overlap + no WNS degradation
  
  5. Final legalization check
```

**Gradient Computation** (timing-aware):
- For each macro position (x, y):
  - Estimate delay change Δt if macro moved by δx, δy
  - Δt computed from fanout, capacitance, timing model
  - If Δt reduces critical slack: gradient points toward move
  
**Why this helps WNS**:
- Directly optimizes worst slack instead of generic wirelength
- Constrained to feasible placements (legality maintained)
- Iterative refinement converges to local timing optimum

---

### Phase 5: Post-Processing & I/O Optimization

**Algorithm**:
```
1. Halo refinement for congestion:
   - Compute routing congestion map around each macro
   - If congestion > threshold: expand halo by 5-10%
   
2. I/O optimization (if applicable for Ariane133):
   - Identify macros with high I/O pin count (e.g., memory macros)
   - Compute distance to nearest I/O pads
   - If distance > threshold: try moving macro closer to pad cluster
   - Accept if area increase < 2% and congestion improves
   
3. Final metrics computation:
   - Total die area
   - Total wirelength
   - WNS on timing critical paths
   - Estimated routing congestion
   - Power density map
```

---

## Implementation Roadmap

### Phase A: Framework & Data Structures (1 week)
- Create `TimingCriticality` analyzer class
- Implement GPU memory management wrappers
- Build `PopulationSolver` base class for multi-core support
- Implement Pareto front tracking

### Phase B: GPU Kernels (2 weeks)
- HPWL computation kernel (highly optimized)
- Macro legalization kernel (sequence packing variant)
- Batch fitness evaluation harness

### Phase C: Population-Based Search (2 weeks)
- Implement island GA with migration
- ML-guided mutation operators
- Population persistence and checkpoint/restart

### Phase D: Timing Refinement (1.5 weeks)
- Slack-flow analyzer
- Timing-aware local search
- STA integration for incremental updates

### Phase E: Integration & Tuning (1 week)
- Tcl command registration (`rtl_macro_placer_hybrid`)
- Parameter experimentation on Ariane133
- Benchmarking vs. current algorithm

**Total effort: ~7-8 weeks**

---

## Expected Improvements on Ariane133/Nangate45

Based on similar hybrid GA + timing-aware approaches in academic literature:

| Metric | Current | Proposed | Improvement |
|--------|---------|----------|-------------|
| **WNS** | baseline | -5% to -15% slack improvement | **Better timing** |
| **Total Area** | baseline | -2% to +3% area | **Comparable/better** |
| **Wirelength** | baseline | -3% to -5% HPWL | **Better routing** |
| **Runtime** | ~5-10 min | ~20 min | +2-4x (but much better quality) |
| **Reproducibility** | deterministic | probabilistic | Requires seed control |

---

## Hardware Utilization Strategy

### CPU (16 cores, 100GB RAM)
- **Cores 0-1**: Timing analysis (STA runs, criticality extraction)
- **Cores 2-17**: Island GA populations (16 independent threads)
- **Memory**: 
  - 60GB for macro placement state (16 populations × 64 solutions)
  - 20GB for STA database
  - 20GB for GPU ↔ CPU transfer buffers

### GPU (RTX 6000 Ada, 48GB VRAM)
- **Kernel 1**: Parallel HPWL computation (1000+ threads)
- **Kernel 2**: Batch legalization (up to 64 solutions simultaneously)
- **Memory**: 
  - 32GB active solution buffers
  - 12GB for distance matrices and lookup tables
  - 4GB for temporary scratch space

### Communication (PCIe 4.0, ~32 GB/sec theoretical)
- Batch send: 64 solutions (~6 MB) → GPU → **~200 µs transfer**
- Batch recv: 64 fitness scores → CPU → **~50 µs transfer**
- Per-generation latency: ~100 ms (dominated by compute, not transfer)

---

## Configuration Parameters for Ariane133

```yaml
# Phase 1: Timing Analysis
TIMING_CRITICALITY_PERCENTILE: 15  # top 15% macros are critical
CRITICAL_PATH_SUBSET_SIZE: 50      # evaluate top 50 critical nets

# Phase 2: Initialization
INITIAL_POPULATION_SIZE: 512       # solutions to generate
TOP_SOLUTIONS_KEPT: 128            # advance to Phase 3

# Phase 3: Multi-Population Search
NUM_ISLANDS: 16                    # one per CPU core
SOLUTIONS_PER_ISLAND: 64
GENERATIONS_PER_ISLAND: 40
MIGRATION_INTERVAL: 10             # generations
MIGRATION_SIZE: 4                  # solutions per island per migration

# Mutation rates
MUTATION_CRITICAL_SWAP_RATE: 0.60  # Type A
MUTATION_POSITION_TARGET_RATE: 0.25 # Type B
MUTATION_ROTATE_RATE: 0.10         # Type C
MUTATION_LOCAL_SEARCH_RATE: 0.05   # Type D
CROSSOVER_RATE: 0.30

# Fitness weights (adaptive based on Phase 1 analysis)
WEIGHT_TIMING: 0.50                # High weight: timing is priority
WEIGHT_AREA: 0.25
WEIGHT_WIRELENGTH: 0.25

# Phase 4: Slack-flow refinement
SLACK_FLOW_ITERATIONS: 100
TIMING_GRADIENT_SCALING: 0.1       # step size for macro moves

# Phase 5: Post-processing
HALO_CONGESTION_THRESHOLD: 0.75    # expand halo if congestion > 75%
IO_DISTANCE_THRESHOLD: 500         # microns, for I/O optimization
```

---

## Comparison with Current Algorithm

| Aspect | Current SA | Proposed TAHPP |
|--------|-----------|----------------|
| **Parallelism** | Single-threaded | 16-way population-level + GPU |
| **Timing awareness** | None (generic HPWL) | Full STA integration |
| **Solution diversity** | Single chain | 1024 solutions/generation |
| **Search strategy** | Greedy cooling | Evolutionary + directed |
| **GPU utilization** | None | Batch evaluation, legalization |
| **Adaptability** | Fixed weights | ML-guided, weight learning |
| **WNS focus** | Indirect | Direct optimization target |

---

## Risks & Mitigations

| Risk | Mitigation |
|------|-----------|
| **Long runtime (20 min)** | Implement early stopping: stop if WNS improves < 1% in last 3 generations |
| **Reproducibility** | Implement deterministic GPU seeding; log all RNG states |
| **GPU memory limits** | Implement solution streaming: process populations in batches |
| **STA overhead** | Cache STA results; only recompute critical paths; use incremental STA |
| **Convergence failures** | Ensemble approach: run 3 independent trials, pick best |

---

## Next Steps

1. **Proof of concept (2 weeks)**: Implement Phases 1-2 with CPU-only, measure on Ariane133
2. **GPU acceleration (2 weeks)**: Port legalization and HPWL to GPU
3. **Full integration (3 weeks)**: Phases 3-5, parameter tuning
4. **Benchmark vs. baseline**: Compare WNS, area, runtime on Ariane133 nangate45
5. **Generalization study**: Test on other designs (aes, jpeg, etc.)

---

## References & Related Work

- **Timing-driven placement**: "Timing-Driven Placement by Partial Delay Relaxation" (ISPD '06)
- **GPU-accelerated EDA**: "CUMPLE: CUDA-accelerated multilevel placement engine" 
- **Population-based search**: "Island Model GA for EDA" (GECCO proceedings)
- **Slack-flow optimization**: "Force-directed placement with timing constraints"

