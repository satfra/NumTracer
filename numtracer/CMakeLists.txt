cmake_minimum_required(VERSION 3.20)
project(NumTracer VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Default to a per-user prefix (same convention as DiFfRG: consumers say
# `find_package(NumTracer REQUIRED HINTS ~/.local/share/NumTracer)`), but never
# override an explicit -DCMAKE_INSTALL_PREFIX.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local/share/NumTracer"
      CACHE PATH "NumTracer install prefix" FORCE)
endif()

# Deliberately NO fast-math flags. `-fno-math-errno` measured ~8% faster on the A4 kernel
# (752 vs 814 ns) and is itself NaN-safe (only skips errno on libm calls), but `-ffast-math`
# / `-ffinite-math-only` assume no NaN/Inf and silently mask them — unacceptable for fRG
# kernels where a NaN signals a real problem. Keep the build NaN-honest; the win is small and
# also speeds the FORM reference, so it does not close the NumTracer-vs-FORM gap anyway.
# CXX-only: once the GPU tests enable_language(CUDA), these must not leak into nvcc
# (-march=native is not an nvcc flag; CUDA targets get their own flags below).
# Directory-scoped, so nothing of it is exported to consumers of the installed package.
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-O3;-march=native;-Wall;-Wextra>")

# Opt-in sanitizers for validation builds (off by default; never affects the shipped
# package). Accepts a semicolon list of sanitizer names, e.g.
#   cmake -S . -B build-asan -DNUMTRACER_SANITIZE="ADDRESS;UNDEFINED"
# Recognised: ADDRESS, UNDEFINED, THREAD, MEMORY (case-insensitive). When enabled we drop
# to -O1 and keep frame pointers so the reports carry clean, inlining-free stack traces;
# the flags go to both compile and link lines (sanitizer runtimes must be linked in).
set(NUMTRACER_SANITIZE "OFF" CACHE STRING
    "Sanitizers for CXX builds: OFF, or a list of ADDRESS;UNDEFINED;THREAD;MEMORY")
if(NOT NUMTRACER_SANITIZE STREQUAL "OFF")
  set(_nt_san_list "")
  foreach(_s IN LISTS NUMTRACER_SANITIZE)
    string(TOLOWER "${_s}" _s_lc)
    if(_s_lc STREQUAL "address")
      list(APPEND _nt_san_list "address")
    elseif(_s_lc STREQUAL "undefined")
      list(APPEND _nt_san_list "undefined")
    elseif(_s_lc STREQUAL "thread")
      list(APPEND _nt_san_list "thread")
    elseif(_s_lc STREQUAL "memory")
      list(APPEND _nt_san_list "memory")
    else()
      message(FATAL_ERROR "NUMTRACER_SANITIZE: unknown sanitizer '${_s}'")
    endif()
  endforeach()
  string(REPLACE ";" "," _nt_san "${_nt_san_list}")
  message(STATUS "NumTracer: sanitizers enabled (-fsanitize=${_nt_san}, -O1)")
  add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-O1;-g;-fno-omit-frame-pointer;-fsanitize=${_nt_san}>")
  add_link_options("-fsanitize=${_nt_san}")
endif()

# ##############################################################################
# The library: header-only interface targets (the installable product)
# ##############################################################################

add_library(NumTracer INTERFACE)
add_library(NumTracer::NumTracer ALIAS NumTracer)
target_compile_features(NumTracer INTERFACE cxx_std_20)
target_include_directories(NumTracer INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)

# Opt-in companion target carrying the raised constexpr budgets. Long contraction
# chains exhaust GCC's default constexpr op/loop limits ("constexpr evaluation
# exceeded ... operations"); consumers hitting that link NumTracer::constexpr_budget
# in addition to NumTracer::NumTracer. Kept off the main target so small consumers
# (and compilers without these flags) are not penalised. Guarded per compiler and
# per language so the flags never reach nvcc or a C compile.
add_library(numtracer_constexpr_budget INTERFACE)
set_target_properties(numtracer_constexpr_budget PROPERTIES EXPORT_NAME constexpr_budget)
add_library(NumTracer::constexpr_budget ALIAS numtracer_constexpr_budget)
target_compile_options(numtracer_constexpr_budget INTERFACE
  "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:GNU>>:-fconstexpr-ops-limit=2000000000;-fconstexpr-loop-limit=100000000;-fconstexpr-depth=2000000>"
  "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:Clang>>:-fconstexpr-steps=2000000000;-fconstexpr-depth=2000000>")

# ##############################################################################
# Install & export (headers + CMake package config; arch-independent)
# ##############################################################################

include(CMakePackageConfigHelpers)

install(DIRECTORY include/ DESTINATION include)
install(TARGETS NumTracer numtracer_constexpr_budget EXPORT NumTracerTargets)
install(EXPORT NumTracerTargets NAMESPACE NumTracer::
        DESTINATION lib/cmake/NumTracer)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/NumTracerConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/NumTracerConfig.cmake
  INSTALL_DESTINATION lib/cmake/NumTracer)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/NumTracerConfigVersion.cmake
  COMPATIBILITY SameMajorVersion ARCH_INDEPENDENT)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/NumTracerConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/NumTracerConfigVersion.cmake
        DESTINATION lib/cmake/NumTracer)

# ##############################################################################
# Mathematica package install (to $UserBaseDirectory/Applications, FunKit-style)
# ##############################################################################
#
# The Wolfram front-end (NumTrace / MakeNTKernel) installs beside FunKit so that
# `Needs["NumTracer`"]` works from anywhere. The Applications dir is queried from
# wolframscript at configure time; if no Wolfram installation is found the install
# is skipped with a warning (the C++ library install is unaffected).
option(NUMTRACER_INSTALL_MATHEMATICA
       "Install the Wolfram package to $UserBaseDirectory/Applications" ON)

if(NUMTRACER_INSTALL_MATHEMATICA)
  find_program(WolframLanguage_wolframscript_EXE NAMES wolframscript wolfram math)
  if(DEFINED ENV{WOLFRAM_SCRIPT})
    set(WolframLanguage_wolframscript_EXE $ENV{WOLFRAM_SCRIPT})
  endif()
  if(WolframLanguage_wolframscript_EXE)
    execute_process(
      COMMAND ${WolframLanguage_wolframscript_EXE} -script
              ${CMAKE_CURRENT_SOURCE_DIR}/cmake/get_wolfram_app_dir.m
      OUTPUT_VARIABLE WOLFRAM_APP_DIR
      OUTPUT_STRIP_TRAILING_WHITESPACE
      RESULT_VARIABLE _wolfram_rc)
    string(REPLACE "\n" "" WOLFRAM_APP_DIR "${WOLFRAM_APP_DIR}")
    if(_wolfram_rc EQUAL 0 AND WOLFRAM_APP_DIR)
      message(STATUS "Wolfram Language application directory: ${WOLFRAM_APP_DIR}")
      # Record the header install location for the package: the "Generate" backend
      # compiles the emitted generator with g++ against these headers.
      configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/NumTracerPaths.m.in
                     ${CMAKE_CURRENT_BINARY_DIR}/NumTracerPaths.m @ONLY)
      install(FILES mathematica/NumTracer.m
                    mathematica/DSL.m
                    mathematica/Frames.m
                    mathematica/Codegen.m
                    mathematica/FunKitAdapter.m
                    mathematica/DiFfRG_compat.m
                    mathematica/PacletInfo.m
                    ${CMAKE_CURRENT_BINARY_DIR}/NumTracerPaths.m
              DESTINATION ${WOLFRAM_APP_DIR}/NumTracer)
    else()
      message(WARNING "wolframscript found but querying $UserBaseDirectory/Applications "
                      "failed — skipping the Mathematica package install")
    endif()
  else()
    message(WARNING "wolframscript not found — skipping the Mathematica package install "
                    "(set the WOLFRAM_SCRIPT environment variable, or configure with "
                    "-DNUMTRACER_INSTALL_MATHEMATICA=OFF to silence this warning)")
  endif()
endif()

# ##############################################################################
# Tests & benchmarks (top-level builds only; off under add_subdirectory)
# ##############################################################################

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(NT_IS_TOP_LEVEL ON)
else()
  set(NT_IS_TOP_LEVEL OFF)
endif()
option(NUMTRACER_BUILD_TESTS "Build the NumTracer tests and benchmarks" ${NT_IS_TOP_LEVEL})

if(NUMTRACER_BUILD_TESTS)
  # On-the-fly color-contraction.
  add_executable(sun_contract tests/sun_contract.cpp)
  target_link_libraries(sun_contract PRIVATE NumTracer::NumTracer)

  # Numeric (matrix-product) contraction backend (numeric/*): MPoly-entry gamma trace via 4x4
  # matrix products + projector index-sum, vs chiral_gamma_trace and direct numeric contraction;
  # validates the monomial-vs-non-monomial projector-atom cancellation. The Part-B correctness gate.
  add_executable(test_numeric_contract tests/test_numeric_contract.cpp)
  target_link_libraries(test_numeric_contract PRIVATE NumTracer::NumTracer)

  # Projector identities for the numeric backend: transverse P_T and longitudinal P_L factors
  # (nprojT / nprojL) checked against the 4x4 oracle and the algebraic identities (P_T+P_L=δ,
  # idempotency, P_T·P_L=0, traces 3/1, transversality/longitudinality).
  add_executable(test_projector_identity tests/test_projector_identity.cpp)
  target_link_libraries(test_projector_identity PRIVATE NumTracer::NumTracer)

  # Projector-algebra fusion (ndetail::fuse_projectors): same-momentum projectors sharing a dummy
  # index collapse pre-expansion — idempotency P·P→P / tr P, orthogonality L⟂{T,E,M}, E⟂M → 0;
  # structural folds (all symmetric slot variants, chains, traces, no-op negatives) + value-preserving.
  add_executable(test_projector_fusion tests/test_projector_fusion.cpp)
  target_link_libraries(test_projector_fusion PRIVATE NumTracer::NumTracer)

  # Dressing-collection layer (numeric/dpoly.hpp + numeric_value_dressed): DPoly arithmetic, the
  # dressing-free reduction (DPoly == MPoly path, same emitted program), and a dressed quark line
  # collected-vs-distributed. The symbolic-dressing-collection correctness gate.
  add_executable(test_dpoly tests/test_dpoly.cpp)
  target_link_libraries(test_dpoly PRIVATE NumTracer::NumTracer)

  # The build-time numeric colour/flavour net (network/sun_net.hpp): fundamental T^a, two groups,
  # vs analytic group factors.
  add_executable(test_col_numeric tests/test_col_numeric.cpp)
  target_link_libraries(test_col_numeric PRIVATE NumTracer::NumTracer)

  # Packed renv layout of the dense momentum builders: a GAPPED mask (0b1001) where the k-th
  # present component lives at renv[Base+k] != renv[Base+component]. Guards the packing the
  # contiguous-mask kernels/tutorials don't exercise.
  add_executable(test_dense_pack tests/test_dense_pack.cpp)
  target_link_libraries(test_dense_pack PRIVATE NumTracer::NumTracer)

  # Phase-2 example 2: pion <- quark loop -> pion (Yukawa) — fundamental flavour x Dirac x colour.
  add_executable(compare_pion_quark tests/refshim/compare_pion_quark.cpp)
  target_include_directories(compare_pion_quark PRIVATE tests/gen)
  target_link_libraries(compare_pion_quark PRIVATE NumTracer::NumTracer NumTracer::constexpr_budget)

  # Phase-2 example 3: gluon <- quark loop -> gluon — fundamental colour via the numeric path.
  add_executable(compare_glu_quark tests/refshim/compare_glu_quark.cpp)
  target_include_directories(compare_glu_quark PRIVATE tests/gen)
  target_link_libraries(compare_glu_quark PRIVATE NumTracer::NumTracer NumTracer::constexpr_budget)

  add_executable(compare_zaqbq1 tests/refshim/compare_zaqbq1.cpp)
  target_include_directories(compare_zaqbq1 PRIVATE tests/gen tests/refshim)
  target_link_libraries(compare_zaqbq1 PRIVATE NumTracer::NumTracer NumTracer::constexpr_budget)

  # Symbolic dressing collection: collected Zq (DressingCollection -> True) vs the distributed Zq.
  add_executable(compare_zq_collect tests/refshim/compare_zq_collect.cpp)
  target_include_directories(compare_zq_collect PRIVATE tests/gen tests/refshim)
  target_link_libraries(compare_zq_collect PRIVATE NumTracer::NumTracer NumTracer::constexpr_budget)

  # NUMERIC (matrix-product) backend validation: ZA quark loop vs Dense oracle; ZA3 quark triangle
  # vs the FORM oracle (the headline — generates in seconds at <0.4 GB where the inv path OOMs).
  add_executable(compare_za_num tests/refshim/compare_za_num.cpp)
  target_include_directories(compare_za_num PRIVATE tests/gen tests/refshim)
  target_link_libraries(compare_za_num PRIVATE NumTracer::NumTracer NumTracer::constexpr_budget)

  add_executable(compare_za3_num tests/refshim/compare_za3_num.cpp)
  target_include_directories(compare_za3_num PRIVATE tests/gen tests/refshim)
  target_link_libraries(compare_za3_num PRIVATE NumTracer::NumTracer NumTracer::constexpr_budget)

  # pure-YM ghost-gluon vertex: numeric vs dense oracle. Guards the disconnected-Lorentz-component
  # (product-not-sum) codegen bug — the only flow with such a diagram.
  add_executable(compare_zacbc_num tests/refshim/compare_zacbc_num.cpp)
  target_include_directories(compare_zacbc_num PRIVATE tests/gen tests/refshim)
  target_link_libraries(compare_zacbc_num PRIVATE NumTracer::NumTracer NumTracer::constexpr_budget)

  # ZA3 with the FULL quark-gluon vertex {1,4,7} tensor basis (the dense struct-7 case). Gate for the
  # struct-7 codegen fix: with ZAqbq4=ZAqbq7=0 the 1/4/7 numeric kernel must reduce EXACTLY to the
  # FORM-validated struct-1 numeric kernel (ZA3_num_kernel), and differ from it with the real 4/7
  # dressings. (The Dense oracle is impractical at this size; struct-1 is the trusted reference.)
  add_executable(compare_za3_147_num tests/refshim/compare_za3_147_num.cpp)
  target_include_directories(compare_za3_147_num PRIVATE tests/gen tests/refshim)
  target_link_libraries(compare_za3_147_num PRIVATE NumTracer::NumTracer NumTracer::constexpr_budget)

  # ZA3 1/4/7 numeric vs the copied FormTracer oracle (independent physical check of structures 4/7).
  add_executable(compare_za3_147_form tests/refshim/compare_za3_147_form.cpp)
  target_include_directories(compare_za3_147_form PRIVATE tests/gen tests/refshim)
  target_link_libraries(compare_za3_147_form PRIVATE NumTracer::NumTracer NumTracer::constexpr_budget)

  # ns/eval timing of the 1/4/7 numeric kernel vs the FORM oracle (build on demand: not a pass/fail
  # test). `cmake --build build --target bench_147 && ./build/bench_147` reproduces the README numbers.
  add_executable(bench_147 EXCLUDE_FROM_ALL tests/refshim/bench_147.cpp)
  target_include_directories(bench_147 PRIVATE tests/gen tests/refshim)
  target_link_libraries(bench_147 PRIVATE NumTracer::NumTracer NumTracer::constexpr_budget)

  # Dense-backend stack-overflow reproducer + ns/eval timing on the largest committed dense kernel
  # (ZA3 1/4/7). Cross-checks the dense oracle against the numeric backend. Build on demand (not a
  # default pass/fail test — it is slow and exercises the heap-backed dense fold):
  #   cmake -S . -B build -DNUMTRACER_BUILD_ZA3_147_DENSE=ON
  #   cmake --build build --target bench_za3_147_dense && ./build/bench_za3_147_dense
  option(NUMTRACER_BUILD_ZA3_147_DENSE "Build the dense stack-overflow reproducer + dense timing harnesses" OFF)
  if(NUMTRACER_BUILD_ZA3_147_DENSE)
    # Stack-overflow reproducer on the largest dense kernel (one eval is tens of seconds).
    add_executable(bench_za3_147_dense tests/refshim/bench_za3_147_dense.cpp)
    target_include_directories(bench_za3_147_dense PRIVATE tests/gen tests/refshim)
    target_link_libraries(bench_za3_147_dense PRIVATE NumTracer::NumTracer NumTracer::constexpr_budget)
    # Fast dense timing on the smaller ZA kernel (~ms/eval) — the tool for measuring dense-fold
    # changes: `NT_DENSE_NP=20000 ./build/bench_za_dense`.
    add_executable(bench_za_dense tests/refshim/bench_za_dense.cpp)
    target_include_directories(bench_za_dense PRIVATE tests/gen tests/refshim)
    target_link_libraries(bench_za_dense PRIVATE NumTracer::NumTracer NumTracer::constexpr_budget)
  endif()

  add_executable(compare_za4_num tests/refshim/compare_za4_num.cpp)
  target_include_directories(compare_za4_num PRIVATE tests/gen tests/refshim)
  target_link_libraries(compare_za4_num PRIVATE NumTracer::NumTracer NumTracer::constexpr_budget)

  # ZA4 (four-gluon box, quark loop) with the FULL {1,4,7} vertex — the heaviest 1/4/7 target, made
  # tractable by the σ^{μν} commutator fold (~14 min generation vs the multi-hour unfolded contraction).
  # Its numeric kernel is ~15 MB and compiles heavily at -O3, so it is NOT committed and this test is
  # OFF by default. Generate the kernel with
  #   wolfram -script tests/gen/gen_qcd_za4_147_numeric.wls
  # (writes tests/gen/ZA4_147_num_kernel{,s}.hh), then configure -DNUMTRACER_BUILD_ZA4_147=ON. Validation
  # is by struct-reduction: with ZAqbq4=ZAqbq7=0 the 1/4/7 kernel must reduce to the FORM-validated
  # ZA4_num_kernel (measured max rel err 8.3e-12), and differ from it with the real 4/7 dressings.
  option(NUMTRACER_BUILD_ZA4_147 "Build the ZA4 1/4/7 struct-reduction test (needs the locally-generated ~15MB kernel)" OFF)
  if(NUMTRACER_BUILD_ZA4_147 AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/gen/ZA4_147_num_kernel.hh")
    add_executable(compare_za4_147_num tests/refshim/compare_za4_147_num.cpp)
    target_include_directories(compare_za4_147_num PRIVATE tests/gen tests/refshim)
    target_link_libraries(compare_za4_147_num PRIVATE NumTracer::NumTracer NumTracer::constexpr_budget)
    add_test(NAME flow_za4_147_num COMMAND compare_za4_147_num)
  endif()

  add_executable(compare_zq_num tests/refshim/compare_zq_num.cpp)
  target_include_directories(compare_zq_num PRIVATE tests/gen tests/refshim)
  target_link_libraries(compare_zq_num PRIVATE NumTracer::NumTracer NumTracer::constexpr_budget)

  # Per-flavour dressings: dressing the u-quark and d-quark DIFFERENTLY (broken flavour symmetry).
  # Two explicit flavour loops, each with its own runtime dressing (Gu/Gd, GuDot/GdDot) and NO
  # SU(Nf) flavour group; numeric kernel vs dense oracle + the per-flavour independence/breaking
  # checks (kernels generated by tests/gen/gen_flavour_split.wls).
  add_executable(compare_flavour_split tests/refshim/compare_flavour_split.cpp)
  target_include_directories(compare_flavour_split PRIVATE tests/gen tests/refshim)
  target_link_libraries(compare_flavour_split PRIVATE NumTracer::NumTracer NumTracer::constexpr_budget)

  # Per-component (group-diagonal) dressings: dress u/d differently WITHIN the SU(2) flavour group
  # (no per-flavour diagram split). Validated by exact equivalence to flow_flavour_split (the
  # separate-diagram kernel), so no dense oracle is needed (kernels by gen/gen_flavour_ingroup.wls).
  add_executable(compare_flavour_ingroup tests/refshim/compare_flavour_ingroup.cpp)
  target_include_directories(compare_flavour_ingroup PRIVATE tests/gen tests/refshim)
  target_link_libraries(compare_flavour_ingroup PRIVATE NumTracer::NumTracer NumTracer::constexpr_budget)

  # Per-component ADJOINT dressing (gluon condensation): an SU(3)-adjoint line dressed component-by-
  # component (Z_A^a, 8 components) folded through the group. Validated by the collapse property
  # (all components equal == the colour-blind kernel) + per-component breaking — no dense oracle
  # needed (kernels by gen/gen_gluon_condensate.wls).
  add_executable(compare_gluon_condensate tests/refshim/compare_gluon_condensate.cpp)
  target_include_directories(compare_gluon_condensate PRIVATE tests/gen tests/refshim)
  target_link_libraries(compare_gluon_condensate PRIVATE NumTracer::NumTracer NumTracer::constexpr_budget)

  # Longitudinal-projector codegen gate: numeric backend vs dense oracle for a hand-built net
  # mixing ntLongProj + ntTransProj (kernels generated by tests/gen/gen_proj_numeric.wls).
  add_executable(compare_proj_num tests/compare_proj_num.cpp)
  target_include_directories(compare_proj_num PRIVATE tests/gen)
  target_link_libraries(compare_proj_num PRIVATE NumTracer::NumTracer)

  # Finite-temperature codegen gate: numeric backend vs dense oracle for a hand-built finite-T net
  # exercising ntElectricProj/ntMagneticProj/ntTransProj + ntSPS + ntVec[q,0] over ftFrameQuark
  # (kernels generated by tests/gen/gen_ftproj_numeric.wls).
  add_executable(compare_ftproj_num tests/compare_ftproj_num.cpp)
  target_include_directories(compare_ftproj_num PRIVATE tests/gen)
  target_link_libraries(compare_ftproj_num PRIVATE NumTracer::NumTracer)

  # Disconnected Dirac x Lorentz codegen gate: numeric backend vs dense oracle for a hand-built net
  # whose diagrams factorize into a closed Dirac trace, a constant colour trace, and an index-disjoint
  # pure-Lorentz scalar that MULTIPLY (kernels generated by tests/gen/gen_discdirac_numeric.wls). The
  # regression guard for MakeNTKernel::disconnectmix (the hSigL sigma-quark vertex): the numeric backend
  # factors the disconnected Lorentz scalar as its own trace and multiplies it in, matching dense.
  add_executable(compare_discdirac_num tests/compare_discdirac_num.cpp)
  target_include_directories(compare_discdirac_num PRIVATE tests/gen)
  target_link_libraries(compare_discdirac_num PRIVATE NumTracer::NumTracer)

  # σL δ-closed meson-loop regression: a quark line threading two pion γ5 vertices closed by a BARE
  # spinor δ external projector (the σL scalar component) — tr(δ·γ5·S·γ5·S) ∝ l1·p1. Guards the
  # orderDiracFacs undirected-walk fix (BUG_meson_sector_dressed_collapse.md): the old directed walk
  # collapsed this δ-closed loop to a single γ5 → tr(γ5)=0. Numeric (gamma-trace theorem) vs the Dense
  # DTensor oracle, which contracts entry-for-entry and is immune to the walk (kernels generated by
  # tests/gen/gen_sigl_dirac.wls).
  add_executable(compare_sigl_dirac tests/refshim/compare_sigl_dirac.cpp)
  target_include_directories(compare_sigl_dirac PRIVATE tests/gen)
  target_link_libraries(compare_sigl_dirac PRIVATE NumTracer::NumTracer NumTracer::constexpr_budget)

  # Finite-temperature quark self-energy (hand-coded, dependency-free): the numeric engine contracts
  # the closed Dirac trace tr[(γ·p) γ^μ (γ·q) γ^ν] against the finite-T gluon line
  # ZAE·P_E(l) + ZAM·P_M(l) in the heat-bath frame (p0 = πT), checked against the closed-form trace
  # identity oracle. (The codegen DENSE oracle overflows the stack at -O3 for a 4-γ × projector trace;
  # this validates the same physics through the numeric path without that fragility.)
  add_executable(test_finiteT_quark tests/test_finiteT_quark.cpp)
  target_link_libraries(test_finiteT_quark PRIVATE NumTracer::NumTracer)

  enable_testing()
  add_test(NAME numeric_contract COMMAND test_numeric_contract)
  add_test(NAME dpoly COMMAND test_dpoly)
  add_test(NAME flow_zq_collect COMMAND compare_zq_collect)
  add_test(NAME projector_identity COMMAND test_projector_identity)
  add_test(NAME projector_fusion COMMAND test_projector_fusion)
  add_test(NAME proj_num COMMAND compare_proj_num)
  add_test(NAME ftproj_num COMMAND compare_ftproj_num)
  add_test(NAME discdirac_num COMMAND compare_discdirac_num)
  add_test(NAME flow_sigl_dirac COMMAND compare_sigl_dirac)
  add_test(NAME finiteT_quark COMMAND test_finiteT_quark)
  add_test(NAME col_numeric COMMAND test_col_numeric)
  add_test(NAME dense_pack COMMAND test_dense_pack)
  add_test(NAME flow_zaqbq1 COMMAND compare_zaqbq1)
  add_test(NAME flow_pion_quark COMMAND compare_pion_quark)
  add_test(NAME flow_glu_quark COMMAND compare_glu_quark)
  add_test(NAME flow_za_num COMMAND compare_za_num)
  add_test(NAME flow_za3_num COMMAND compare_za3_num)
  add_test(NAME flow_zacbc_num COMMAND compare_zacbc_num)
  add_test(NAME flow_za3_147_num COMMAND compare_za3_147_num)
  add_test(NAME flow_za3_147_form COMMAND compare_za3_147_form)
  add_test(NAME flow_za4_num COMMAND compare_za4_num)
  add_test(NAME flow_zq_num COMMAND compare_zq_num)
  add_test(NAME flow_flavour_split COMMAND compare_flavour_split)
  add_test(NAME flow_flavour_ingroup COMMAND compare_flavour_ingroup)
  add_test(NAME flow_gluon_condensate COMMAND compare_gluon_condensate)
  add_test(NAME color_su3 COMMAND sun_contract 3)
endif()

# ##############################################################################
# GPU integration tests (CUDA + GSL; OFF by default so the core stays dependency-free)
# ##############################################################################
#
# Full 4D loop integrals of the generated ZA3/ZA4 kernels over a grid of external momenta,
# DiFfRG-style two-phase map/reduce on the GPU (include/numtracer/cuda/integrator.cuh).
# Needs the __host__ __device__ kernel variants (tests/gen/ZA{3,4}_gpu_kernel*.hh, regenerate
# with tests/gen/gen_qcd_za{3,4}_gpu.wls), the CUDA toolkit and GSL (quadrature nodes/weights).
# Configure with -DNUMTRACER_GPU_TESTS=ON; if the default host compiler is too new for nvcc,
# add -DCMAKE_CUDA_HOST_COMPILER=g++-14. clang CUDA also works (-DCMAKE_CUDA_COMPILER=clang++).
option(NUMTRACER_GPU_TESTS "Build the CUDA GPU loop-integral integration tests" OFF)
# Same tests expressed in Kokkos (CUDA backend) instead of raw CUDA + cub, for a framework
# overhead comparison. Uses the DiFfRG-bundled Kokkos install (override NT_KOKKOS_ROOT).
option(NUMTRACER_KOKKOS_TESTS "Build the Kokkos twins of the GPU integration tests" OFF)

if(NUMTRACER_GPU_TESTS OR NUMTRACER_KOKKOS_TESTS)
  find_package(GSL REQUIRED)
  set(NT_CONSTEXPR_FLAGS -fconstexpr-ops-limit=2000000000 -fconstexpr-loop-limit=100000000
                         -fconstexpr-depth=2000000)
  # CPU same-quadrature reference (the committed HOST kernels), shared by the CUDA and Kokkos
  # test mains and always compiled as plain C++ by the configured CXX compiler.
  add_library(nt_gpu_ref STATIC tests/gpu/ref_za3.cpp tests/gpu/ref_za4.cpp)
  # tests/gpu BEFORE tests/gen: the generated kernels' `#include "shim.hpp"` must resolve to
  # the host/device shim, not tests/refshim's host-only one (which is not on this path at all).
  target_include_directories(nt_gpu_ref PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/tests/gpu
                                               ${CMAKE_CURRENT_SOURCE_DIR}/tests/gen)
  target_link_libraries(nt_gpu_ref PUBLIC NumTracer::NumTracer GSL::gsl)
  target_compile_options(nt_gpu_ref PRIVATE $<$<CXX_COMPILER_ID:GNU>:${NT_CONSTEXPR_FLAGS}>)
endif()

if(NUMTRACER_GPU_TESTS)
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES native)
  endif()
  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 20)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)

  function(nt_gpu_test name)
    add_executable(${name} ${ARGN})
    target_link_libraries(${name} PRIVATE nt_gpu_ref)
    target_compile_options(${name} PRIVATE
      # CUDA TU via nvcc: relaxed constexpr (the kernels' compile-time SU(3) colour nets call
      # constexpr host functions from __device__ code) + the budgets for the gcc host pass
      # (-Xcompiler: the gcc host pass of every __host__ __device__ function evaluates the same
      #  constexpr nets, so it needs the budget too)
      "$<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<CUDA_COMPILER_ID:NVIDIA>>:--expt-relaxed-constexpr;-O3;-lineinfo;-Xcompiler=-O3,-fconstexpr-ops-limit=2000000000,-fconstexpr-loop-limit=100000000,-fconstexpr-depth=2000000>"
      # CUDA TU via clang: it evaluates constexpr itself, so raise its step limits
      "$<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<CUDA_COMPILER_ID:Clang>>:-O3;-fconstexpr-steps=2000000000;-fconstexpr-depth=2000000>")
    add_test(NAME ${name} COMMAND ${name})
    set_tests_properties(${name} PROPERTIES LABELS gpu)
  endfunction()

  nt_gpu_test(integrate_za3_gpu tests/gpu/integrate_za3_gpu.cu)
  nt_gpu_test(integrate_za4_gpu tests/gpu/integrate_za4_gpu.cu)
endif()

if(NUMTRACER_KOKKOS_TESTS)
  set(NT_KOKKOS_ROOT "$ENV{HOME}/.local/share/DiFfRG/bundled" CACHE PATH
      "Kokkos install prefix (default: the DiFfRG bundled install)")
  find_package(Kokkos REQUIRED PATHS ${NT_KOKKOS_ROOT} NO_DEFAULT_PATH)

  function(nt_kokkos_test name src)
    add_executable(${name} ${src})
    target_link_libraries(${name} PRIVATE nt_gpu_ref Kokkos::kokkos)
    # the Kokkos TU still needs the gcc constexpr budgets (nvcc_wrapper forwards unknown
    # flags to the g++-14 host pass); Kokkos::kokkos itself injects --expt-relaxed-constexpr,
    # --extended-lambda and the arch flags.
    target_compile_options(${name} PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${NT_CONSTEXPR_FLAGS}>)
    # route this target's C++ compiles through kokkos_launch_compiler -> nvcc_wrapper(g++-14)
    kokkos_compilation(TARGET ${name})
    add_test(NAME ${name} COMMAND ${name})
    set_tests_properties(${name} PROPERTIES LABELS "gpu;kokkos")
  endfunction()

  nt_kokkos_test(integrate_za3_kokkos tests/gpu/integrate_za3_kokkos.cpp)
  nt_kokkos_test(integrate_za4_kokkos tests/gpu/integrate_za4_kokkos.cpp)
endif()

# ##############################################################################
# Documentation
# ##############################################################################
#
# A Sphinx site bundling the C++ API (Doxygen wrapped by doxysphinx) and the
# narrative pages. It is orchestrated by documentation/build.sh, which needs
# Doxygen plus a Python 3.9–3.12 environment with the packages in
# documentation/requirements.txt. OFF by default so an ordinary build never
# needs the docs toolchain; we also degrade gracefully when it is missing.
option(NumTracer_DOCUMENTATION "Build the documentation site" OFF)

if(NumTracer_DOCUMENTATION)
  find_package(Doxygen)
  find_program(SPHINX_BUILD_EXECUTABLE sphinx-build)
  find_program(DOXYSPHINX_EXECUTABLE doxysphinx)

  if(DOXYGEN_FOUND AND SPHINX_BUILD_EXECUTABLE AND DOXYSPHINX_EXECUTABLE)
    message(STATUS "Documentation target configured (Doxygen + doxysphinx + Sphinx)")
    add_custom_target(
      documentation
      COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/documentation/build.sh
              ${CMAKE_CURRENT_BINARY_DIR}/documentation
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/documentation
      COMMENT "Building NumTracer documentation (C++ API + narrative)"
      VERBATIM)
  else()
    # Report exactly which tools are missing so the user knows what to supply.
    set(_doc_missing "")
    if(NOT DOXYGEN_FOUND)
      list(APPEND _doc_missing "doxygen")
    endif()
    if(NOT SPHINX_BUILD_EXECUTABLE)
      list(APPEND _doc_missing "sphinx-build")
    endif()
    if(NOT DOXYSPHINX_EXECUTABLE)
      list(APPEND _doc_missing "doxysphinx")
    endif()
    string(REPLACE ";" ", " _doc_missing_str "${_doc_missing}")
    message(WARNING
      "The 'documentation' target was NOT created: missing ${_doc_missing_str}.\n"
      "  * doxygen — install via your package manager (e.g. 'pacman -S doxygen').\n"
      "  * sphinx-build / doxysphinx — come from the Python toolchain. doxysphinx requires\n"
      "    Python 3.9-3.12 (NOT 3.13+), so create such an environment and install the deps:\n"
      "        python3.12 -m venv ${CMAKE_CURRENT_SOURCE_DIR}/documentation/.venv-docs\n"
      "        source ${CMAKE_CURRENT_SOURCE_DIR}/documentation/.venv-docs/bin/activate\n"
      "        pip install -r ${CMAKE_CURRENT_SOURCE_DIR}/documentation/requirements.txt\n"
      "    then re-run cmake with that environment ACTIVE so the tools are on PATH.\n"
      "  Or skip cmake entirely and build the docs directly:\n"
      "        bash ${CMAKE_CURRENT_SOURCE_DIR}/documentation/build.sh\n"
      "  Configure with -DNumTracer_DOCUMENTATION=OFF to silence this warning.")
  endif()
endif()
