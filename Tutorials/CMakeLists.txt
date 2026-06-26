# The NumTracer tutorial programs — the runnable versions of the code in
# numtracer/documentation/tutorials/*.md, one subdirectory per tutorial page.
#
# This is a STANDALONE CMake project, deliberately not wired into numtracer/'s
# own build: it consumes NumTracer exactly like an external user would, via
#   find_package(NumTracer REQUIRED HINTS ~/.local/share/NumTracer)
# (install with `cmake --install numtracer/build`). When no installed NumTracer
# is found it falls back to the in-tree library, so the tutorials also build
# straight from a repo checkout:
#   cmake -S . -B build && cmake --build build       # then ./build/<tutorial>
#   ctest --test-dir build                            # runs every tutorial
#
# Every program depends ONLY on NumTracer::NumTracer — no test fixtures, no
# Wolfram — and self-checks against an oracle or closed form (exit non-zero on
# mismatch), so `ctest` is a meaningful gate.
cmake_minimum_required(VERSION 3.20)
project(NumTracerTutorials CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(NumTracer QUIET HINTS $ENV{HOME}/.local/share/NumTracer)
if(NumTracer_FOUND)
  message(STATUS "Using installed NumTracer: ${NumTracer_DIR}")
elseif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../numtracer/CMakeLists.txt)
  message(STATUS "NumTracer not installed — falling back to the in-tree library")
  set(NUMTRACER_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  set(NUMTRACER_INSTALL_MATHEMATICA OFF CACHE BOOL "" FORCE)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../numtracer numtracer EXCLUDE_FROM_ALL)
else()
  message(FATAL_ERROR "NumTracer not found: install it (cmake --install numtracer/build) "
                      "or build inside the repository checkout")
endif()

# The docs build the tutorials with -O2 and without -march=native: this project
# mimics an ordinary external consumer, not the tuned in-tree test build.
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-O2;-Wall;-Wextra>")

enable_testing()

function(nt_tutorial name src)
  add_executable(${name} ${src})
  target_link_libraries(${name} PRIVATE NumTracer::NumTracer ${ARGN})
  add_test(NAME ${name} COMMAND ${name})
endfunction()

# 01 — colour factors fold to a number
nt_tutorial(colour_factors 01-color-contraction/colour_factors.cpp)

# 02 — the Dirac trace, hand-coded: raw 4x4 matrix product, then the token API
nt_tutorial(trace_raw    02-dirac-traces/trace_raw.cpp)
nt_tutorial(trace_tokens 02-dirac-traces/trace_tokens.cpp)

# 03 — Lorentz networks reduce to scalar products
nt_tutorial(lorentz_networks 03-lorentz-networks/lorentz_networks.cpp)

# 04 — hand-code a full diagram, validated two ways (numeric / closed form)
nt_tutorial(full_diagram 04-full-diagram/full_diagram.cpp)

# 05 — generating kernels from Mathematica: a walkthrough, not a compiled program.
#      See 05-generating-kernels/README.md and tests/gen/gen_zq_numeric.wls. The
#      optional target below regenerates the committed numeric Zq kernel (needs
#      wolframscript + FunKit); the committed header is what a consumer build uses.
find_program(NT_WOLFRAMSCRIPT NAMES wolframscript wolfram math)
if(NT_WOLFRAMSCRIPT AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../numtracer/tests/gen/gen_zq_numeric.wls)
  add_custom_target(tutorial_codegen
    COMMAND ${NT_WOLFRAMSCRIPT} -script
            ${CMAKE_CURRENT_SOURCE_DIR}/../numtracer/tests/gen/gen_zq_numeric.wls
    COMMENT "Regenerating tests/gen/Zq_num_kernel.hh via the Mathematica codegen"
    VERBATIM)
endif()
