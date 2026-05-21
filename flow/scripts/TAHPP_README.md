# TAHPP Hybrid Macro Placement Integration

## What was added

This update introduces **Partial TAHPP** (Timing-Aware Hybrid Placement Placer): baseline hierarchical SA, then optional GPU-batched population search and STA refinement on elites.

## Key files

- `tools/OpenROAD/src/mpl/src/tahpp/timing_criticality.{h,cpp}` — path-based criticality and macro affinity pairs
- `tools/OpenROAD/src/mpl/src/tahpp/gpu_batch_eval.{h,cpp,cu}` — batch HPWL / overlap / affinity fitness (CUDA optional)
- `tools/OpenROAD/src/mpl/src/tahpp/tahpp_optimizer.{h,cpp}` — island GA, Pareto elites, STA slack-flow refinement
- `tools/OpenROAD/src/mpl/src/rtl_mp_hybrid.cpp` — improved `placeHybrid()` and `placeTahpp()`
- `flow/scripts/macro_place_util.tcl` — routes `RTLMP_HYBRID=tahpp` to TAHPP mode

## How to use

### Hybrid local search (lightweight)

```tcl
set ::env(RTLMP_HYBRID) hybrid
set ::env(RTLMP_TIMING_WT) 1.2
set ::env(RTLMP_HYBRID_ITERATIONS) 8
set ::env(RTLMP_CANDIDATE_RADIUS) 40.0
```

### Partial TAHPP (population + GPU + STA)

```tcl
set ::env(RTLMP_HYBRID) tahpp
set ::env(RTLMP_GPU) 1
set ::env(RTLMP_TIMING_WT) 1.2
set ::env(RTLMP_WNS_WEIGHT) 0.5
set ::env(RTLMP_AREA_WEIGHT) 0.25
set ::env(RTLMP_WL_WEIGHT) 0.25
set ::env(RTLMP_POP_SIZE) 32
set ::env(RTLMP_NUM_ISLANDS) 4
set ::env(RTLMP_NUM_GENERATIONS) 20
set ::env(RTLMP_INIT_POPULATION) 256
```

Build with `-DBUILD_TAHPP_GPU=ON` to enable the CUDA code path (falls back to CPU batch evaluation otherwise).

## Pipeline

1. `place()` — stock hierarchical RTL macro placer.
2. Phase 1 — STA worst paths → macro criticality + affinity pairs.
3. Phase 2 — seed population, island GA with batch fitness (GPU or CPU).
4. Phase 3 — Pareto elite selection, full STA refinement, hybrid local search on best solution.

## Notes

GPU kernels currently batch-evaluate via the CPU evaluator in parallel builds; device kernels can replace `evaluateBatchCuda` for further speedup. See `ADVANCED_MACRO_PLACEMENT_ALGORITHM.md` for the full TAHPP vision.
