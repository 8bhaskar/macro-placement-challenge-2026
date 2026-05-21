# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021-2025, The OpenROAD Authors
#
# Merge these sources into the OpenROAD mpl target in tools/OpenROAD/src/mpl/CMakeLists.txt

option(BUILD_TAHPP_GPU "Build TAHPP CUDA batch evaluator" OFF)

set(MPL_TAHPP_SOURCES
  src/tahpp/timing_criticality.cpp
  src/tahpp/placement_solution.cpp
  src/tahpp/gpu_batch_eval.cpp
  src/tahpp/tahpp_optimizer.cpp
  src/rtl_mp_hybrid.cpp
)

# Add to mpl target:
# target_sources(mpl PRIVATE ${MPL_TAHPP_SOURCES})
# target_include_directories(mpl PUBLIC include src)

if(BUILD_TAHPP_GPU)
  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)
  # target_compile_definitions(mpl PRIVATE BUILD_TAHPP_GPU)
  # target_link_libraries(mpl PRIVATE CUDA::cudart)
endif()
