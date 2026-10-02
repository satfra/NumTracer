# The NumTracer tutorial programs — the runnable versions of the code in
# numtracer/documentation/tutorials/step-*.md, one subdirectory per tutorial
# step (step-NN-<slug>/, matching step-NN.md 1:1).
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
# The pages quote these files with Sphinx `literalinclude`, pulling named
# regions delimited by `// @snip begin: <name>` / `// @snip end: <name>`
# comments — so the documentation can never drift from the code that runs.
#
# THREE TIERS OF RUNNABILITY (see documentation/tutorials/index.md):
#   A  C++20 compiler only ............ nt_tutorial();     always a ctest test
#   B  + a Wolfram kernel and FunKit .. nt_tutorial_wls(); a ctest test when Wolfram is found
#   C  + TensorBases and DiFfRG ....... not built here; the pages give the commands
#
# Tier B needs FunKit because MakeNTKernel delegates all scalar/boilerplate C++ emission to
# FunKit's COEN. It does NOT need TensorBases, FORM or DiFfRG: the .wls scripts here build their
# DSL networks by hand.
#
# Every tier-A/B program self-checks against an oracle or closed form (exit
# non-zero on mismatch), so `ctest` is a meaningful gate.
cmake_minimum_required(VERSION 3.20)
project(NumTracerTutorials CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(NumTracer QUIET HINTS $ENV{HOME}/.local/share/NumTracer)
if(NumTracer_FOUND)
  message(STATUS "Using installed NumTracer: ${NumTracer_DIR}")
  set(NT_SOURCE_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/..")
elseif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../numtracer/CMakeLists.txt)
  message(STATUS "NumTracer not installed — falling back to the in-tree library")
  set(NUMTRACER_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  set(NUMTRACER_INSTALL_MATHEMATICA OFF CACHE BOOL "" FORCE)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../numtracer numtracer EXCLUDE_FROM_ALL)
  set(NT_SOURCE_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/..")
else()
  message(FATAL_ERROR "NumTracer not found: install it (cmake --install numtracer/build) "
                      "or build inside the repository checkout")
endif()

# The docs build the tutorials with -O2 and without -march=native: this project
# mimics an ordinary external consumer, not the tuned in-tree test build.
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-O2;-Wall;-Wextra>")

enable_testing()

# ---------------------------------------------------------------------------
# Tier A — a self-checking C++ program, built and tested unconditionally.
# ---------------------------------------------------------------------------
function(nt_tutorial name src)
  add_executable(${name} ${src})
  target_link_libraries(${name} PRIVATE NumTracer::NumTracer ${ARGN})
  add_test(NAME ${name} COMMAND ${name})
endfunction()

# ---------------------------------------------------------------------------
# Tier B — a Mathematica codegen driver plus a C++ checker for what it emitted.
#
#   nt_tutorial_wls(<name> <script.wls> <checker.cpp> <emitted-header>...)
#
# The .wls scripts build their DSL network BY HAND, so they need a Wolfram kernel
# and FunKit but no TensorBases and no FORM. Each writes its kernel headers
# into ${CMAKE_CURRENT_BINARY_DIR}/gen/, which the checker includes; the checker
# validates the generated kernel against a closed form or the runtime engine.
#
# Generation is a build-time custom command, so `cmake --build . --target <name>`
# regenerates only when the script changes. Without wolframscript the target is
# skipped entirely (with a status message) and no test is registered.
# ---------------------------------------------------------------------------
find_program(NT_WOLFRAMSCRIPT NAMES wolframscript wolfram math)
if(NT_WOLFRAMSCRIPT)
  message(STATUS "Wolfram kernel found (${NT_WOLFRAMSCRIPT}) — tier-B codegen tutorials enabled")
else()
  message(STATUS "No wolframscript/wolfram/math on PATH — tier-B codegen tutorials skipped")
endif()

add_custom_target(tutorials_codegen COMMENT "Running every tier-B Mathematica codegen tutorial")

function(nt_tutorial_wls name script checker)
  if(NOT NT_WOLFRAMSCRIPT)
    return()
  endif()
  set(_gen ${CMAKE_CURRENT_BINARY_DIR}/gen/${name})
  # Created at CONFIGURE time: it is the custom command's WORKING_DIRECTORY, so the build-time
  # `cd` into it happens before any COMMAND could have made it.
  file(MAKE_DIRECTORY ${_gen})
  set(_outs "")
  foreach(h IN LISTS ARGN)
    list(APPEND _outs ${_gen}/${h})
  endforeach()

  add_custom_command(
    OUTPUT ${_outs}
    # Delete the previous outputs FIRST. wolframscript exits 0 even on a syntax error in the
    # script, so without this a broken .wls leaves the stale headers in place and the build goes
    # green against the previous run's kernels — a false pass of exactly the kind the codegen
    # tests warn about. With the outputs removed, a script that fails to write them fails the build.
    COMMAND ${CMAKE_COMMAND} -E rm -f ${_outs}
    # NT_TUTORIAL_OUT tells the script where to write; NUMTRACER_DIR lets it
    # Get[] the front-end straight out of the checkout when it is not installed.
    COMMAND ${CMAKE_COMMAND} -E env
            NT_TUTORIAL_OUT=${_gen}
            NUMTRACER_DIR=${NT_SOURCE_ROOT}/numtracer
            ${NT_WOLFRAMSCRIPT} -script ${CMAKE_CURRENT_SOURCE_DIR}/${script}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${script}
    WORKING_DIRECTORY ${_gen}
    COMMENT "Generating kernels for ${name} (${script})"
    VERBATIM)

  add_custom_target(${name}_gen DEPENDS ${_outs})
  add_dependencies(tutorials_codegen ${name}_gen)

  add_executable(${name} ${checker} ${_outs})
  target_include_directories(${name} PRIVATE ${_gen})
  target_link_libraries(${name} PRIVATE NumTracer::NumTracer)
  # A generated kernel emits every declared parameter whether or not this particular flow uses it
  # (the signature is the flow's contract, not a minimal set), so an unused `p` in constant() is
  # normal and not a defect. Real consumers do the same — see QCD_Nf2's flows/CMakeLists.txt.
  target_compile_options(${name} PRIVATE -Wno-unused-parameter -Wno-unused-variable)
  add_test(NAME ${name} COMMAND ${name})
endfunction()

# ===========================================================================
# I. The engine, by hand — tier A
# ===========================================================================
# step-01 — hello, tensor network: contract vectors through a metric and a projector (no physics)
nt_tutorial(hello_network step-01-hello-network/hello_network.cpp)

# step-02 — colour factors fold to a number
nt_tutorial(colour_factors step-02-colour-factors/colour_factors.cpp)

# step-03 — the Dirac trace, hand-coded: raw 4x4 matrix product, then the token API
nt_tutorial(trace_raw    step-03-dirac-traces/trace_raw.cpp)
nt_tutorial(trace_tokens step-03-dirac-traces/trace_tokens.cpp)

# step-04 — Lorentz networks reduce to scalar products
nt_tutorial(lorentz_networks step-04-lorentz-networks/lorentz_networks.cpp)

# step-05 — hand-code a full diagram, validated two ways, then lower it to a kernel
nt_tutorial(full_diagram step-05-full-diagram/full_diagram.cpp)

# ===========================================================================
# IV. Advanced sectors — tier A
# ===========================================================================
# step-16 — finite temperature: the electric/magnetic split and a broken-O(4) trace
nt_tutorial(finite_temperature step-16-finite-temperature/finite_temperature.cpp)

# step-17 — per-flavour and per-component dressings: the group-diagonal SU(N) fold
nt_tutorial(dressed_flavour step-17-dressed-flavour/dressed_flavour.cpp)

# step-18 — dressed propagator numerators: one trace collected into a DPoly, not 2^D diagrams
nt_tutorial(dressed_numerators step-18-dressed-numerators/dressed_numerators.cpp)

# ===========================================================================
# II. The front-end — tier B (needs wolframscript + FunKit)
# ===========================================================================
# step-06 — the smallest end-to-end generation: a DSL network in, a C++ kernel out
nt_tutorial_wls(first_kernel
                step-06-first-kernel/first_kernel.wls
                step-06-first-kernel/check_first_kernel.cpp
                first_kernel.hh first_kernels.hh)

# step-07 — frames: the general 3-point frame must reduce to the symmetric-point one
nt_tutorial_wls(frames
                step-07-frames/frames.wls
                step-07-frames/check_frames.cpp
                sp3_kernel.hh sp3_kernels.hh gen3_kernel.hh gen3_kernels.hh)

# step-08 — one micro-kernel per DSL head, each checked against its closed form
nt_tutorial_wls(dsl_heads
                step-08-dsl-heads/dsl_heads.wls
                step-08-dsl-heads/check_dsl_heads.cpp
                H_metric_kernel.hh H_metric_kernels.hh
                H_transproj_kernel.hh H_transproj_kernels.hh
                H_longproj_kernel.hh H_longproj_kernels.hh
                H_diracid_kernel.hh H_diracid_kernels.hh
                H_dirac2_kernel.hh H_dirac2_kernels.hh
                H_dirac4_kernel.hh H_dirac4_kernels.hh
                H_sunT_kernel.hh H_sunT_kernels.hh
                H_sunf_kernel.hh H_sunf_kernels.hh
                H_dress_kernel.hh H_dress_kernels.hh)

# step-09 — the same network emitted three ways; every variant must agree to the last bit
nt_tutorial_wls(options
                step-09-codegen-options/options.wls
                step-09-codegen-options/check_options.cpp
                Opt_default_kernel.hh Opt_default_kernels.hh
                Opt_const_kernel.hh Opt_const_kernels.hh
                Opt_ns_kernel.hh Opt_ns_kernels.hh)

# ===========================================================================
# IV. Advanced sectors — tier B
# ===========================================================================
# step-19 — fixed-component gammas (the feature four-quark Fierz bases need)
nt_tutorial_wls(fixed_components
                step-19-four-fermi/fixed_components.wls
                step-19-four-fermi/check_fixed_components.cpp
                Fx_t_kernel.hh Fx_t_kernels.hh
                Fx_x_kernel.hh Fx_x_kernels.hh
                Fx_tt_kernel.hh Fx_tt_kernels.hh
                Fx_4_kernel.hh Fx_4_kernels.hh
                Fx_4ref_kernel.hh Fx_4ref_kernels.hh
                Fx_g_kernel.hh Fx_g_kernels.hh)
