# TAHPP Hybrid Macro Placement Integration

## What was added

This update introduces a new hybrid macro placement entry point for OpenROAD, based on the existing `mpl` macro placer.
It is designed to start the implementation of a timing-aware hybrid macro placement flow for Ariane133 on Nangate45.

## Key files

- `dashboard/gui/resources/tools/OpenROAD/src/mpl/include/mpl/rtl_mp.h`
  - Added a new `MacroPlacer::placeHybrid()` API.
- `dashboard/gui/resources/tools/OpenROAD/src/mpl/src/rtl_mp.cpp`
  - Implemented the new `placeHybrid()` algorithm.
  - Added helper cost functions and a local macro move pass.
- `dashboard/gui/resources/tools/OpenROAD/src/mpl/src/mpl.i`
  - Exported the new `rtl_macro_placer_hybrid_cmd` SWIG wrapper.
- `dashboard/gui/resources/tools/OpenROAD/src/mpl/src/mpl.tcl`
  - Added a new Tcl command `rtl_macro_placer_hybrid`.
- `flow/scripts/macro_place_util.tcl`
  - Added flow-level support for selecting the hybrid placer and an optional virtual global placement warm start.
- `flow/scripts/macro_place.tcl`
  - Added a macro placement WNS/timing report stage after macro placement.
- `flow/scripts/variables.yaml`
  - Documented the new hybrid macro placement environment variables.

## How to use

The new hybrid path is enabled by setting an environment variable in the flow:

```tcl
set ::env(RTLMP_HYBRID) 1
```

Optional tuning variables:

- `RTLMP_VIRTUAL_PLACE` - enable a lightweight `global_placement` warm start before macro placement.
- `RTLMP_TIMING_WT` - maps to `-timing_weight`
- `RTLMP_HYBRID_ITERATIONS` - maps to `-hybrid_iterations`
- `RTLMP_CANDIDATE_RADIUS` - maps to `-candidate_radius`

Example:

```tcl
set ::env(RTLMP_HYBRID) 1
set ::env(RTLMP_TIMING_WT) 1.2
set ::env(RTLMP_HYBRID_ITERATIONS) 8
set ::env(RTLMP_CANDIDATE_RADIUS) 40.0
```

## What the new algorithm does

- Calls the existing RTL macro placer for a baseline placement.
- Uses STA-driven timing pressure on nets connected to each macro to choose the highest-risk macros.
- Computes a timing-weighted candidate cost using HPWL plus negative-slack pressure.
- Performs a limited local search on top macros to improve timing while respecting macro overlap and block boundaries.

## Where the files belong

These files must be kept in the OpenROAD `mpl` plugin area, not in the wrapper flow scripts:

- OpenROAD plugin source and interface: `dashboard/gui/resources/tools/OpenROAD/src/mpl`
- Flow integration and environment control: `flow/scripts/macro_place_util.tcl`

## Notes

This is an initial implementation to make the hybrid algorithm runnable and reproducible.
The current cost model is a lightweight proxy suitable for early algorithm development.
Further work can refine the timing model, expand the population search, and add GPU-accelerated candidate evaluation.
