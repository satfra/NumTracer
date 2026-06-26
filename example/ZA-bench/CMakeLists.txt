# ZA-bench — standalone benchmark comparing the NumTracer full-QCD AAqbq ZA kernel against the
# FormTracer analogue (generation time, emitted size, compile time, runtime ns/eval), and a
# faithfulness gate that the global-collection rework is output-preserving vs the pre-rework path.
#
# This consumes already-generated kernels under generated/{form,nt,nt_before}/ (populate them with
# run_bench.sh / gen_*.sh for live regeneration). It uses ONLY the in-tree NumTracer headers + a
# standalone shim — no DiFfRG / Kokkos / CUDA. Configure & build:
#   cmake -S . -B build && cmake --build build && ctest --test-dir build --output-on-failure
# Full live-regeneration + size/compile/runtime report:  ./run_bench.sh
cmake_minimum_required(VERSION 3.20)
project(ZA_bench CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
add_compile_options(-O3 -march=native -ftemplate-depth=4000)

# in-tree NumTracer headers (the kernels #include numtracer/...) + the shim + generated kernels
set(NT_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/../../numtracer/include" CACHE PATH "NumTracer include dir")
set(BENCH_INCLUDES
  "${CMAKE_CURRENT_SOURCE_DIR}/bench"
  "${CMAKE_CURRENT_SOURCE_DIR}/generated"
  "${CMAKE_CURRENT_SOURCE_DIR}/bench/shiminc"
  "${NT_INCLUDE}")

enable_testing()

# Build a compare binary per NT variant whose kernel.hh exists.
function(add_variant variant hdr cls)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/generated/${variant}/kernel.hh")
    add_executable(compare_za_${variant} bench/compare_za_full.cpp)
    target_include_directories(compare_za_${variant} PRIVATE ${BENCH_INCLUDES})
    target_compile_definitions(compare_za_${variant} PRIVATE
      NT_HDR="${hdr}" NT_CLASS=${cls} NT_NAME="${variant}")
    add_test(NAME za_${variant} COMMAND compare_za_${variant})
  endif()
endfunction()

add_variant(nt        "nt/kernel.hh"        ZA_nt_kernel)
add_variant(nt_before "nt_before/kernel.hh" ZA_nt_before_kernel)

# convenience: full live-regeneration + size/compile/runtime report
add_custom_target(bench COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_bench.sh
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} USES_TERMINAL)
