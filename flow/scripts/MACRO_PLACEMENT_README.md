# Macro Placement README

This document explains the files involved in the macro placement stage of the OpenROAD flow and how they interact.

## Overview

Macro placement in this flow is driven by a Tcl stage in `flow/scripts` and implemented by the OpenROAD macro placer extension located under `dashboard/gui/resources/tools/OpenROAD/src/mpl`.

The high-level flow is:
1. `flow/scripts/macro_place.tcl` loads the stage and sources helper logic.
2. `flow/scripts/macro_place_util.tcl` configures macro placer environment variables and invokes the `rtl_macro_placer` command.
3. OpenROAD's `mpl` package supplies the `rtl_macro_placer` command and the underlying hierarchical simulated annealing placement algorithm.
4. After placement, `flow/scripts/placement_blockages.tcl` creates blockages around macros.

## Tcl files

### `flow/scripts/macro_place.tcl`

Role:
- Main macro placement stage script in the flow.
- Loads the design checkpoint and SDC for the current floorplan stage.
- Sources the macro placement utility script.
- Runs `report_design_area`, writes out a post-placement database, and exports macro placement to a Tcl file.

What it does:
- `source $::env(SCRIPTS_DIR)/load.tcl`
- `erase_non_stage_variables floorplan`
- `load_design 2_1_floorplan.odb 2_1_floorplan.sdc`
- `source_step_tcl PRE MACRO_PLACE`
- `source $::env(SCRIPTS_DIR)/macro_place_util.tcl`
- `source_step_tcl POST MACRO_PLACE`
- `report_design_area`
- `orfs_write_db $::env(RESULTS_DIR)/2_2_floorplan_macro.odb`
- `write_macro_placement $::env(RESULTS_DIR)/2_2_floorplan_macro.tcl`

Why it matters:
- It is the entry point for macro placement in the flow.
- It guarantees the design is loaded and the placement stage is recorded as a distinct stage.

### `flow/scripts/macro_place_util.tcl`

Role:
- Configures and launches the macro placement engine.
- Sets default report paths and runtime options.
- Applies macro wrappers if configured.
- Builds the Tcl command arguments for the OpenROAD placer.
- Invokes `rtl_macro_placer`.
- Creates macro blockages after placement.

Key behavior:
- Checks for macros using `[find_macros]` and skips placement if none are found.
- Sets default environment variables:
  - `RTLMP_RPT_DIR` → report directory for the macro placer
  - `RTLMP_RPT_FILE` → report filename (default `partition.txt`)
  - `RTLMP_BLOCKAGE_FILE` → macro blockage file path
- Supports macro wrappers via the optional `MACRO_WRAPPERS` script.
- Reads `MACRO_PLACE_HALO` and optionally `MACRO_BLOCKAGE_HALO` to configure placement halos and blockages.
- Sources an optional external placement configuration file via `MACRO_PLACEMENT_TCL`.
- Builds placement arguments from `RTLMP_*` environment variables and defaults, including:
  - `-max_num_level`
  - `-max_num_inst`
  - `-min_num_inst`
  - `-max_num_macro`
  - `-min_num_macro`
  - `-halo_width`
  - `-halo_height`
  - `-min_ar`
  - `-area_weight`
  - `-wirelength_weight`
  - `-outline_weight`
  - `-boundary_weight`
  - `-notch_weight`
  - `-report_directory`
  - `-fence_lx`, `-fence_ly`, `-fence_ux`, `-fence_uy`
  - `-target_util`
- If `RTLMP_ARGS` is defined, it overrides constructed args.
- Runs `rtl_macro_placer {*}$all_args`.
- Calls `source $::env(SCRIPTS_DIR)/placement_blockages.tcl` and `block_channels $blockage_width`.

Why it matters:
- This script is the flow-level translator between environment variables and the OpenROAD macro placer command.
- It is where macro-specific flow policies are defined.

### `flow/scripts/placement_blockages.tcl`

Role:
- Converts placed macro bounding boxes into soft blockages.
- Expands macro shapes by the configured channel width.
- Clips the blockage region to the core area.
- Creates `dbBlockage` objects for the resulting expanded macro regions.

Key behavior:
- Iterates all instances in the block and selects macros.
- Makes a shape set from macro bounding boxes.
- Bloats the shapes by `channel_width_in_microns * units`.
- Intersects with core area.
- Creates soft blockages from resulting rectangles.

Why it matters:
- It prevents subsequent placement/routing passes from moving standard cells into macro channel regions.
- It is executed immediately after `rtl_macro_placer`.

## C++ files

The actual macro placement algorithm is implemented in OpenROAD's `mpl` package under:
- `dashboard/gui/resources/tools/OpenROAD/src/mpl/include/mpl/`
- `dashboard/gui/resources/tools/OpenROAD/src/mpl/src/`

### `dashboard/gui/resources/tools/OpenROAD/src/mpl/src/MakeMacroPlacer.cpp`

Role:
- Integration entry point for the OpenROAD Tcl extension.
- Declares `extern int Mpl_Init(Tcl_Interp* interp)`.
- Implements `mpl::initMacroPlacer` which calls `Mpl_Init` and loads Tcl initialization scripts.

Why it matters:
- This file wires the macro placer into the OpenROAD Tcl interpreter.
- It is the bridge between Tcl and the C++ macro placement library.

### `dashboard/gui/resources/tools/OpenROAD/src/mpl/include/mpl/rtl_mp.h`

Role:
- Declares the `mpl::MacroPlacer` class API.
- Defines methods for:
  - macro placement execution (`place`)
  - placing an individual macro (`placeMacro`)
  - detecting overlaps
  - guidance regions and halos
  - debug controls

Why it matters:
- It documents the runtime interface used by the macro placement engine.
- It shows the parameters exposed to the flow layer.

### `dashboard/gui/resources/tools/OpenROAD/src/mpl/src/rtl_mp.cpp`

Role:
- Implements the `MacroPlacer` C++ wrapper.
- Receives flow-level parameters and forwards them to the hierarchical placer.
- Configures cluster sizes, halo, fences, weights, and report settings.
- Calls `hier_rtlmp_->run()`.
- Includes utility functions for placing individual macros and checking overlap.

Why it matters:
- This is the runtime implementation invoked by `rtl_macro_placer`.
- It maps Tcl-level arguments to algorithm-level settings.

### `dashboard/gui/resources/tools/OpenROAD/src/mpl/src/hier_rtlmp.cpp`

Role:
- Implements the hierarchical macro placement algorithm.
- Provides setters for all placement and clustering weights.
- Implements `run()` and the top-level placement workflow.
- Coordinates algorithm stages:
  1. Multilevel autoclustering
  2. Coarse shaping
  3. Fine shaping via cluster placement
  4. Hierarchical macro placement
  5. Boundary pushing
  6. Orientation improvement

Why it matters:
- This file contains the core algorithm flow for macro placement.
- It is the primary organizer of the algorithmic stages used by the placer.

### `dashboard/gui/resources/tools/OpenROAD/src/mpl/src/SimulatedAnnealingCore.cpp`

Role:
- Implements simulated annealing for cluster and macro placement.
- Provides sequence pair initialization and optimization.
- Computes cost functions such as area, wirelength, guidance, fence, and outline penalties.
- Supports random perturbations and move operators used by the placer.

Why it matters:
- This is the key optimization engine for macro placement quality.
- It controls how the placer searches for improved macro arrangements.

### `dashboard/gui/resources/tools/OpenROAD/src/mpl/src/clusterEngine.cpp`

Role:
- Implements hierarchical clustering of macros and standard cells.
- Converts logical hierarchy into a physical clustering tree.
- Computes macro with halo area, floorplan shape, and macro placement feasibility.
- Breaks mixed clusters and creates cluster structures for the SA optimizer.

Why it matters:
- It prepares the problem for hierarchical placement.
- It determines clustering and shape information before the placement stage.

## File relationships

- `macro_place.tcl` is the user-visible flow stage that triggers macro placement.
- `macro_place_util.tcl` performs macro-specific configuration and invokes the low-level `rtl_macro_placer` command.
- `placement_blockages.tcl` creates post-placement soft blockages for the macro channel regions.
- `MakeMacroPlacer.cpp` registers the macro placer with the OpenROAD Tcl interpreter.
- `rtl_mp.cpp` receives the command parameters and executes the placer.
- `hier_rtlmp.cpp` orchestrates the high-level macro placement algorithm.
- `SimulatedAnnealingCore.cpp` performs the optimization search used by the placer.
- `clusterEngine.cpp` creates the cluster hierarchy and validates the macro placement problem.

## Notes

- The Tcl files are flow-level orchestration scripts.
- The C++ files are algorithmic implementation and runtime integration.
- The `rtl_macro_placer` command is exposed through the OpenROAD `mpl` extension and is the key runtime entry point.

If you need a next-level trace, the exact command is launched from `macro_place_util.tcl` with `rtl_macro_placer {*}$all_args` and the implementation is rooted in `dashboard/gui/resources/tools/OpenROAD/src/mpl/src/rtl_mp.cpp`.
