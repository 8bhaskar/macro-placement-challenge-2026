# Macro Placement Challenge 2026

This repository contains our **timing-aware hybrid macro placement** algorithm integrated into the [OpenROAD-flow-scripts](https://github.com/The-OpenROAD-Project/OpenROAD-flow-scripts) (ORFS) flow. The layout here mirrors the paths used in a standard ORFS checkout so you can copy files into your existing setup without restructuring your tree.

The algorithm extends OpenROAD’s `mpl` macro placer with a hybrid refinement pass (`rtl_macro_placer_hybrid`) that runs after the baseline hierarchical placer and applies timing-weighted local search on critical macros.

---

## What is in this repo

| Path in this repo | Typical location in your ORFS tree | Purpose |
|-------------------|-------------------------------------|---------|
| `flow/scripts/macro_place.tcl` | `flow/scripts/macro_place.tcl` | Macro placement stage (loads design, invokes placer, writes results) |
| `flow/scripts/macro_place_util.tcl` | `flow/scripts/macro_place_util.tcl` | Flow glue: env vars → `rtl_macro_placer` / `rtl_macro_placer_hybrid` |
| `flow/scripts/variables.yaml` | `flow/scripts/variables.yaml` | Documents hybrid tuning variables (`RTLMP_TIMING_WT`, etc.) |
| `tools/OpenROAD/src/mpl/include/mpl/rtl_mp.h` | `tools/OpenROAD/src/mpl/include/mpl/rtl_mp.h` | C++ API for `placeHybrid()` |
| `tools/OpenROAD/src/mpl/src/rtl_mp.cpp` | `tools/OpenROAD/src/mpl/src/rtl_mp.cpp` | Hybrid algorithm implementation |
| `tools/OpenROAD/src/mpl/src/mpl.tcl` | `tools/OpenROAD/src/mpl/src/mpl.tcl` | Tcl command `rtl_macro_placer_hybrid` |
| `tools/OpenROAD/src/mpl/src/mpl.i` | `tools/OpenROAD/src/mpl/src/mpl.i` | SWIG bindings for the hybrid command |
| `tools/OpenROAD/src/mpl/src/rtl_mp_hybrid.cpp` | `tools/OpenROAD/src/mpl/src/rtl_mp_hybrid.cpp` | Improved hybrid + TAHPP orchestration |
| `tools/OpenROAD/src/mpl/src/tahpp/*` | `tools/OpenROAD/src/mpl/src/tahpp/` | Partial TAHPP (timing map, GPU batch eval, island GA, STA refine) |
| `tools/OpenROAD/src/mpl/CMakeLists.txt` | merge into `tools/OpenROAD/src/mpl/CMakeLists.txt` | Optional `BUILD_TAHPP_GPU` CUDA support |

Additional design notes live under `flow/scripts/` (for example `TAHPP_README.md`, `MACRO_PLACEMENT_README.md`).

---

## Prerequisites

- A working **OpenROAD-flow-scripts** installation (local clone or container image).
- Ability to **rebuild OpenROAD** after updating the `mpl` sources (the flow Tcl changes alone are not enough).
- A design that reaches the **floorplan macro placement** stage (`2_2_floorplan_macro` in the default ORFS numbering).

---

## Step 1: Back up your current files

From the root of your ORFS checkout, save copies of the files you will replace:

```bash
ORFS_ROOT=/path/to/OpenROAD-flow-scripts   # adjust to your install
PATCH_ROOT=/path/to/macro-placement-challenge-2026

cd "$ORFS_ROOT"
mkdir -p .backup-macro-placement-2026/flow/scripts
mkdir -p .backup-macro-placement-2026/tools/OpenROAD/src/mpl

cp flow/scripts/macro_place.tcl flow/scripts/macro_place_util.tcl flow/scripts/variables.yaml \
   .backup-macro-placement-2026/flow/scripts/ 2>/dev/null || true

cp tools/OpenROAD/src/mpl/include/mpl/rtl_mp.h \
   tools/OpenROAD/src/mpl/src/rtl_mp.cpp \
   tools/OpenROAD/src/mpl/src/mpl.tcl \
   tools/OpenROAD/src/mpl/src/mpl.i \
   .backup-macro-placement-2026/tools/OpenROAD/src/mpl/ 2>/dev/null || true
```

---

## Step 2: Copy flow scripts into ORFS

Copy the Tcl/YAML files from this repo into your ORFS `flow/scripts` directory:

```bash
ORFS_ROOT=/path/to/OpenROAD-flow-scripts
PATCH_ROOT=/path/to/macro-placement-challenge-2026

cp "$PATCH_ROOT/flow/scripts/macro_place.tcl" \
   "$PATCH_ROOT/flow/scripts/macro_place_util.tcl" \
   "$ORFS_ROOT/flow/scripts/"

# variables.yaml is large; replace only if you have not customized it locally.
# If you have local edits, merge the RTLMP_* hybrid entries from our variables.yaml
# (RTLMP_VIRTUAL_PLACE, RTLMP_TIMING_WT, RTLMP_HYBRID_ITERATIONS, RTLMP_CANDIDATE_RADIUS).
cp "$PATCH_ROOT/flow/scripts/variables.yaml" "$ORFS_ROOT/flow/scripts/"
```

These scripts are picked up automatically when you run the normal ORFS floorplan flow; `SCRIPTS_DIR` already points at `flow/scripts`.

---

## Step 3: Copy OpenROAD `mpl` sources and rebuild

The hybrid placer is implemented inside OpenROAD’s `mpl` tool. Copy the patched files into the OpenROAD submodule (or standalone OpenROAD tree) under ORFS:

```bash
ORFS_ROOT=/path/to/OpenROAD-flow-scripts
PATCH_ROOT=/path/to/macro-placement-challenge-2026
MPL="$ORFS_ROOT/tools/OpenROAD/src/mpl"

cp "$PATCH_ROOT/tools/OpenROAD/src/mpl/include/mpl/rtl_mp.h" \
   "$MPL/include/mpl/"

cp "$PATCH_ROOT/tools/OpenROAD/src/mpl/src/rtl_mp.cpp" \
   "$PATCH_ROOT/tools/OpenROAD/src/mpl/src/mpl.tcl" \
   "$PATCH_ROOT/tools/OpenROAD/src/mpl/src/mpl.i" \
   "$MPL/src/"
```

Then merge the TAHPP sources into the mpl CMake target (see `tools/OpenROAD/src/mpl/CMakeLists.txt` in this repo) and rebuild OpenROAD:

```bash
cd "$ORFS_ROOT"
# In tools/OpenROAD/src/mpl/CMakeLists.txt add:
#   src/tahpp/timing_criticality.cpp
#   src/tahpp/placement_solution.cpp
#   src/tahpp/gpu_batch_eval.cpp
#   src/tahpp/tahpp_optimizer.cpp
#   src/rtl_mp_hybrid.cpp
# Optional GPU: -DBUILD_TAHPP_GPU=ON and src/tahpp/gpu_batch_eval.cu

./etc/Build.sh
# or: make openroad
```

Use the same build method you normally use for ORFS so the `openroad` binary on your `PATH` includes the new `rtl_macro_placer_hybrid` command.

**Verify the install** (optional):

```bash
openroad -exit -python <<'EOF'
# In an OpenROAD Tcl shell you can also run: help rtl_macro_placer_hybrid
EOF
```

Or start `openroad` interactively and run `help rtl_macro_placer_hybrid`. If the command is missing, the `mpl` sources were not rebuilt or the wrong binary is on your `PATH`.

---

## Step 4: Enable hybrid or TAHPP

By default, ORFS still calls the stock `rtl_macro_placer`. Set `RTLMP_HYBRID` to select the algorithm:

| `RTLMP_HYBRID` | Algorithm |
|----------------|-----------|
| `hybrid` / `1` / `true` | Baseline SA + path-based affinity local search (orientations, adaptive radius) |
| `tahpp` | Baseline SA + GPU/CPU population GA + Pareto elites + STA slack refinement |

### Option A — design `config.mk` or `config.tcl`

Add to your design configuration (for example `designs/<platform>/<design>/config.mk` or the Tcl/env file you already use):

```tcl
export RTLMP_HYBRID=1
export RTLMP_TIMING_WT=1.2
export RTLMP_HYBRID_ITERATIONS=8
export RTLMP_CANDIDATE_RADIUS=40.0
# Optional warm start:
# export RTLMP_VIRTUAL_PLACE=1
```

### Option B — `env.sh` or shell before `make`

```bash
export RTLMP_HYBRID=1
export RTLMP_TIMING_WT=1.2
export RTLMP_HYBRID_ITERATIONS=8
export RTLMP_CANDIDATE_RADIUS=40.0
cd "$ORFS_ROOT/designs/nangate45/ariane133"   # example
make
```

### Option C — one-off `make` argument

```bash
make RTLMP_HYBRID=1 RTLMP_TIMING_WT=1.2 RTLMP_HYBRID_ITERATIONS=8
```

### Partial TAHPP (recommended for WNS + GPU)

```bash
export RTLMP_HYBRID=tahpp
export RTLMP_GPU=1              # requires BUILD_TAHPP_GPU
export RTLMP_TIMING_WT=1.2
export RTLMP_WNS_WEIGHT=0.5
export RTLMP_AREA_WEIGHT=0.25
export RTLMP_WL_WEIGHT=0.25
export RTLMP_POP_SIZE=32
export RTLMP_NUM_ISLANDS=4
export RTLMP_NUM_GENERATIONS=20
export RTLMP_INIT_POPULATION=256
export RTLMP_GPU_BATCH=64
make
```

---

## Step 5: Run macro placement

Run the flow as you normally would. Macro placement is the floorplan substage that produces `2_2_floorplan_macro.odb`.

Examples:

```bash
# Full flow from floorplan through macro place
make

# Run only through macro placement (stage name may vary slightly by ORFS version)
make floorplan
# or, depending on your Makefile targets:
make do-2_2_floorplan_macro
```

With `RTLMP_HYBRID=1`, `macro_place_util.tcl` invokes **`rtl_macro_placer_hybrid`** instead of **`rtl_macro_placer`**. The hybrid flow:

1. Runs the standard hierarchical RTL macro placer for a baseline legal placement.
2. Updates timing (STA) and scores macros by timing pressure on connected nets.
3. Performs a limited local search (iterations and radius controlled by env vars) to improve timing while respecting overlap and die boundaries.

Outputs are unchanged from stock ORFS: `results/.../2_2_floorplan_macro.odb` and `2_2_floorplan_macro.tcl`, plus RTLMP reports under `objects/rtlmp/` when configured.

---

## Environment variables (quick reference)

| Variable | Default (in our `variables.yaml`) | Effect |
|----------|-----------------------------------|--------|
| `RTLMP_HYBRID` | off | `hybrid`/`1` for local search; `tahpp` for population + STA refinement |
| `RTLMP_GPU` | `0` | Enable GPU batch path when `RTLMP_HYBRID=tahpp` |
| `RTLMP_POP_SIZE` | `32` | Island population size (TAHPP) |
| `RTLMP_NUM_ISLANDS` | `4` | GA islands (TAHPP) |
| `RTLMP_NUM_GENERATIONS` | `20` | GA generations (TAHPP) |
| `RTLMP_INIT_POPULATION` | `256` | Initial seeded population (TAHPP) |
| `RTLMP_WNS_WEIGHT` | `0.5` | Timing affinity fitness weight |
| `RTLMP_AREA_WEIGHT` | `0.25` | Floorplan area fitness weight |
| `RTLMP_WL_WEIGHT` | `0.25` | HPWL fitness weight |
| `RTLMP_TIMING_WT` | `1.0` | Timing penalty weight (`-timing_weight`) |
| `RTLMP_HYBRID_ITERATIONS` | `4` | Refinement iterations after baseline placement |
| `RTLMP_CANDIDATE_RADIUS` | `40.0` | Local search radius in microns |
| `RTLMP_VIRTUAL_PLACE` | `0` | If `1`, run lightweight `global_placement` before macro place |
| `RTLMP_ARGS` | — | If set, overrides all constructed placer arguments |
| `RTLMP_*` (max level, weights, fence, etc.) | see `variables.yaml` | Passed through to both baseline and hybrid placers |

Standard ORFS macro variables (`MACRO_PLACE_HALO`, `MACRO_PLACEMENT_TCL`, `MACRO_WRAPPERS`, …) still apply.

---

## Reverting to stock OpenROAD macro placement

1. Restore files from `.backup-macro-placement-2026/` (Step 1).
2. Rebuild OpenROAD (`./etc/Build.sh` or `make openroad`).
3. Unset `RTLMP_HYBRID` or set it to `0`.

---

## Troubleshooting

| Symptom | Likely cause | What to do |
|---------|----------------|------------|
| `invalid command name "rtl_macro_placer_hybrid"` | OpenROAD not rebuilt or wrong binary | Rebuild after Step 3; confirm `which openroad` |
| Flow runs but placement looks unchanged | `RTLMP_HYBRID` not set | Export `RTLMP_HYBRID=1` in config or shell |
| `variables.yaml` merge conflicts | Local ORFS customizations | Merge only the new `RTLMP_*` keys from our file |
| Build errors in `mpl` | ORFS/OpenROAD version skew | Align ORFS submodule with the version you developed against, then re-apply patches |

---

## Further reading

- `flow/scripts/TAHPP_README.md` — hybrid integration overview  
- `flow/scripts/MACRO_PLACEMENT_README.md` — how Tcl stages connect to `mpl`  
- `flow/scripts/TAHPP_DECISION_GUIDE.md` — when to use hybrid vs stock placer  

---

## License

See [LICENSE](LICENSE). OpenROAD-derived files retain their original SPDX headers where present.
