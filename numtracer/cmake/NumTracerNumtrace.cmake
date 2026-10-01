# NumTracerNumtrace.cmake — offline kernel generation as a build target.
#
#   numtracer_add_numtrace(<flows_target> <flow_dir>)
#
# Emitting a flow from Wolfram writes three things: the generator sources (into a gen/ directory that
# is normally gitignored), a probe source when the flow is complex, and a per-flow
# <flow_dir>/<Name>/numtrace.json manifest. This function turns the manifests into a `numtrace` target
# that compiles and runs whatever still needs generating, and makes <flows_target> depend on it — so
# the flows library can never be built against stale or missing kernels.
#
# The switch is the manifest's "generated" field, not a file timestamp. It has to be: kernels.hh is
# committed while gen/ is not, so a fresh clone has the outputs and none of the inputs, and a
# timestamp rule would either demand sources that do not exist or silently accept stale ones. A flow
# at 1 contributes nothing here — if every flow is at 1, `numtrace` is an empty target and does
# nothing, which is the normal state of a checked-out tree.
#
# Flows are generated ONE AT A TIME, and each gets the whole job budget. The generator binary sizes
# its own thread pool from hardware_concurrency(), so it is not one job however make counts it: left
# to itself, `make -j8` would start eight generators of 32 workers each, and a phase-A worker
# transiently holds ~1 GB on the dense 4-point flows. Two halves fix that:
#
#   * flow i+1's generator TARGET depends on flow i's completion stamp, so nothing of flow i+1 —
#     not even its compiles — starts until flow i has run, probed and marked. Within a flow the
#     unit TUs still compile across all -j slots, and then the run has the machine to itself.
#   * the driver resolves that same -jN (out of MAKEFLAGS, which GNU make exports to every recipe)
#     and hands it to the generator as NT_GEN_MAXW / NT_GEN_MAXW_B, capped by the flow's own
#     "maxw"/"maxw_b" if the notebook set SetNumTracerThreads around it.
#
# So `make -j8 numtrace` means eight compiler processes, then one generator with eight workers, per
# flow, in sorted order. A consequence to know about: because the edge sits on the generator target,
# asking for one flow (`make numtrace_<ns>`) also builds every pending flow ahead of it in that
# order. Expressing "mutually exclusive if both are scheduled" needs Ninja job pools, which the
# Makefile generator has no equivalent of.
#
# Per flow the chain is
#
#   gen_<ns>_num.cpp + _u*.cpp  --compile+link-->  numtrace_gen_<ns>
#                                      --run-->   <Name>/kernels.hh
#                                   --compile-->  numtrace_probe_<ns>      (complex flows only)
#                                      --run-->   <Name>/numtrace_verdict.hh
#
# and the last step in the chain also flips the switch, so a flow is marked done exactly once and
# only after everything it owes has succeeded. The verdict header carries the imaginary-part probe's
# three-way answer; kernel.hh holds all three bodies under an #if and the preprocessor picks one.

if(COMMAND numtracer_add_numtrace)
  return()
endif()

# Worker count for build tools that do not report their -j in the environment — Ninja above all.
# Under `make` this is never consulted: the driver reads the real -jN out of MAKEFLAGS at build
# time, which is the only place it exists (a configure-time value would be wrong the moment you
# changed -j without reconfiguring).
if(NOT DEFINED NUMTRACE_JOBS)
  cmake_host_system_information(RESULT _nt_cores QUERY NUMBER_OF_PHYSICAL_CORES)
  if(NOT _nt_cores OR _nt_cores LESS 1)
    set(_nt_cores 4)
  endif()
  set(NUMTRACE_JOBS "${_nt_cores}" CACHE STRING
      "numtrace: generator worker threads when the build tool does not report -jN (e.g. Ninja)")
  unset(_nt_cores)
endif()

function(numtracer_add_numtrace flows_target flow_dir)
  set(_driver "${NumTracer_NUMTRACE_DIR}/NumTracerNumtraceRun.cmake")
  if(NOT EXISTS "${_driver}")
    message(FATAL_ERROR "numtrace: driver script not found at ${_driver}")
  endif()

  file(GLOB _manifests "${flow_dir}/*/numtrace.json")

  # Sorted, because the flows are now run in this order and a glob-order chain would shuffle a long
  # sequential build from machine to machine.
  list(SORT _manifests)

  # A re-emit from Wolfram rewrites the manifests (switch back to 0). Registering them here makes the
  # next build reconfigure by itself, so the new work is picked up without anyone remembering to
  # re-run cmake. The glob is registered too, so an entirely new flow directory also triggers it.
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${_manifests})

  # How many flows this call will generate, for the "[k/N] tracing <flow>" banner. Known only after
  # a pass over the switches, and the banner is baked into the command line, so it needs its own
  # cheap pre-pass rather than being counted as we go.
  set(_total 0)
  foreach(_manifest IN LISTS _manifests)
    file(READ "${_manifest}" _json)
    string(JSON _generated GET "${_json}" "generated")
    if(NOT _generated)
      math(EXPR _total "${_total} + 1")
    endif()
  endforeach()
  set(_idx 0)

  set(_outputs "")
  set(_pending "")

  foreach(_manifest IN LISTS _manifests)
    get_filename_component(_flowdir "${_manifest}" DIRECTORY)
    file(READ "${_manifest}" _json)

    string(JSON _generated GET "${_json}" "generated")
    if(_generated)
      continue()
    endif()

    string(JSON _name      GET "${_json}" "name")
    string(JSON _ns        GET "${_json}" "namespace")
    string(JSON _gendir    GET "${_json}" "gen_dir")
    string(JSON _generator GET "${_json}" "generator")
    string(JSON _decorator GET "${_json}" "decorator")
    string(JSON _mainopt   GET "${_json}" "main_opt")
    string(JSON _complex   GET "${_json}" "complex")
    string(JSON _kernels   GET "${_json}" "kernels")

    math(EXPR _idx "${_idx} + 1")

    # Per-flow worker caps, from whatever SetNumTracerThreads was in force when the flow was emitted
    # (0 = unset = take the build's -jN). Absent from every manifest written before these fields
    # existed, and string(JSON GET) is a hard error on a missing key — hence the ERROR_VARIABLE form.
    string(JSON _maxw  ERROR_VARIABLE _jerr GET "${_json}" "maxw")
    if(_jerr OR NOT _maxw)
      set(_maxw 0)
    endif()
    string(JSON _maxwb ERROR_VARIABLE _jerr GET "${_json}" "maxw_b")
    if(_jerr OR NOT _maxwb)
      set(_maxwb 0)
    endif()

    # Does this flow's kernel target device code? Enables gen.hpp's size-gated `__noinline__`, which
    # the generator reads as NT_GEN_DEVICE. Same ERROR_VARIABLE form and the same reason as the caps
    # above: absent from every manifest written before the field existed, and a missing key is a hard
    # error for string(JSON GET). Absent => false => the old all-inline emission, so old manifests
    # keep building exactly as they did.
    string(JSON _device ERROR_VARIABLE _jerr GET "${_json}" "device")
    if(_jerr OR NOT _device)
      set(_device 0)
    endif()

    # Manifest paths are normally relative to the flow directory, but tolerate absolute ones: a path
    # that could not be made relative at emit time (e.g. a not-yet-created file whose symlinked parent
    # resolved differently) is written out absolute, and must not be re-prefixed.
    if(NOT IS_ABSOLUTE "${_gendir}")
      set(_gendir "${_flowdir}/${_gendir}")
    endif()
    get_filename_component(_gendir "${_gendir}" ABSOLUTE)

    # unit sources
    string(JSON _nunits LENGTH "${_json}" "units")
    set(_sources "${_gendir}/${_generator}")
    math(EXPR _last "${_nunits} - 1")
    if(_nunits GREATER 0)
      foreach(_i RANGE ${_last})
        string(JSON _unit GET "${_json}" "units" ${_i})
        list(APPEND _sources "${_gendir}/${_unit}")
      endforeach()
    endif()

    # A committed switch of 0 with no sources on disk is unbuildable, and the resulting CMake error
    # would otherwise be an obscure "cannot find source file". Say what actually has to happen.
    foreach(_src IN LISTS _sources)
      if(NOT EXISTS "${_src}")
        message(FATAL_ERROR
          "numtrace: flow '${_name}' is marked ungenerated but its generator source is missing:\n"
          "  ${_src}\n"
          "Re-run the flow notebook to re-emit it, or set \"generated\": 1 in\n"
          "  ${_manifest}\n"
          "if the committed kernels are known to be current.")
      endif()
    endforeach()

    # ---- the generator -------------------------------------------------------------------------
    set(_gentgt "numtrace_gen_${_ns}")
    add_executable(${_gentgt} EXCLUDE_FROM_ALL ${_sources})

    # THE SERIALIZATION. A target-level edge onto the previous flow's completion target, so this
    # flow's unit TUs cannot start compiling until that one has run, probed and marked. Both the
    # Makefile and Ninja generators honour it, and it applies to EXCLUDE_FROM_ALL targets. The tail
    # is a GLOBAL property, not a local: a project with several flow directories calls this function
    # once per directory, and the chain has to continue across those calls rather than restart.
    get_property(_prev GLOBAL PROPERTY NUMTRACE_CHAIN_TAIL)
    if(_prev)
      add_dependencies(${_gentgt} ${_prev})
    endif()

    target_link_libraries(${_gentgt} PRIVATE NumTracer::NumTracer)
    target_compile_options(${_gentgt} PRIVATE -ftemplate-depth=4000)
    set_target_properties(${_gentgt} PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)

    # The main TU carries the contraction loop and must stay optimised or the RUN slows ~10x; the
    # net-builder units are pure straight-line codegen where -O0 is both faster to compile and
    # runtime-neutral. -fno-exceptions -fno-rtti on the units drops the exception-cleanup landing pads
    # for their thousands of destructible braced-init temporaries (~90% of the -O0 unit compile); the
    # unit code never throws. Same split and flags the Wolfram generator compile uses.
    set_source_files_properties("${_gendir}/${_generator}"
      TARGET_DIRECTORY ${_gentgt} PROPERTIES COMPILE_OPTIONS "${_mainopt}")
    foreach(_src IN LISTS _sources)
      if(NOT _src STREQUAL "${_gendir}/${_generator}")
        set_source_files_properties("${_src}"
          TARGET_DIRECTORY ${_gentgt} PROPERTIES COMPILE_OPTIONS "-O0;-fno-exceptions;-fno-rtti")
      endif()
    endforeach()

    if(IS_ABSOLUTE "${_kernels}")
      set(_kernels_hh "${_kernels}")
    else()
      set(_kernels_hh "${_flowdir}/${_kernels}")
    endif()
    set(_run_cmd
      COMMAND ${CMAKE_COMMAND}
        -DMODE=run
        "-DGEN=$<TARGET_FILE:${_gentgt}>"
        "-DOUT=${_kernels_hh}"
        "-DNS=${_ns}"
        "-DFLOW=${_name}"
        "-DDECOR=${_decorator}"
        "-DJOBS=${NUMTRACE_JOBS}"
        "-DMAXW=${_maxw}"
        "-DMAXWB=${_maxwb}"
        "-DDEVICE=${_device}"
        "-DIDX=${_idx}"
        "-DTOTAL=${_total}"
        -P "${_driver}")

    # The traces header, for every flow. USES_TERMINAL: with the flows serialized, the build spends
    # nearly all its wall clock inside one generator, so its NT_GEN_PROFILE phase lines have to
    # stream rather than surface at the end. On Ninja it also lands the step in the depth-1 console
    # pool, reinforcing the chain.
    add_custom_command(
      OUTPUT "${_kernels_hh}"
      ${_run_cmd}
      DEPENDS ${_gentgt}
      COMMENT "NumTracer: [${_idx}/${_total}] tracing ${_name}"
      USES_TERMINAL
      VERBATIM)

    set(_needs_probe FALSE)
    if(_complex AND _json MATCHES "\"probe\"[ \t\r\n]*:")
      set(_needs_probe TRUE)
    endif()

    if(_needs_probe)
      # then the probe that reads it; the switch flips at the end of the chain.
      string(JSON _probe   GET "${_json}" "probe")
      string(JSON _macro   GET "${_json}" "verdict_macro")
      string(JSON _verdict GET "${_json}" "verdict")

      set(_probetgt "numtrace_probe_${_ns}")
      add_executable(${_probetgt} EXCLUDE_FROM_ALL "${_gendir}/${_probe}")
      target_include_directories(${_probetgt} PRIVATE "${_flowdir}")
      target_compile_options(${_probetgt} PRIVATE -O1 -w)
      set_target_properties(${_probetgt} PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)
      # the probe #includes the traces header, so it cannot be compiled before the generator has run
      set_source_files_properties("${_gendir}/${_probe}"
        TARGET_DIRECTORY ${_probetgt} PROPERTIES OBJECT_DEPENDS "${_kernels_hh}")

      if(IS_ABSOLUTE "${_verdict}")
        set(_verdict_hh "${_verdict}")
      else()
        set(_verdict_hh "${_flowdir}/${_verdict}")
      endif()
      add_custom_command(
        OUTPUT "${_verdict_hh}"
        COMMAND ${CMAKE_COMMAND}
          -DMODE=probe
          "-DPROBE=$<TARGET_FILE:${_probetgt}>"
          "-DOUT=${_verdict_hh}"
          "-DMACRO=${_macro}"
          "-DFLOW=${_name}"
          -P "${_driver}"
        DEPENDS ${_probetgt} "${_kernels_hh}"
        COMMENT "NumTracer: probing ${_name}"
        VERBATIM)

      set(_last "${_verdict_hh}")
      list(APPEND _outputs "${_kernels_hh}" "${_verdict_hh}")
    else()
      set(_last "${_kernels_hh}")
      list(APPEND _outputs "${_kernels_hh}")
    endif()

    # The switch flip gets its OWN output, a stamp in the build tree, rather than riding along as a
    # trailing COMMAND of the rule that produces kernels.hh. Make deletes a rule's output when the
    # recipe fails, so a failure in the flip would otherwise delete the committed traces header that
    # the run step had already written. The stamp is expendable; the header is not.
    set(_stampdir "${CMAKE_CURRENT_BINARY_DIR}/numtrace")
    set(_stamp "${_stampdir}/${_ns}.marked")
    file(MAKE_DIRECTORY "${_stampdir}")
    add_custom_command(
      OUTPUT "${_stamp}"
      COMMAND ${CMAKE_COMMAND} -DMODE=mark "-DMANIFEST=${_manifest}" -P "${_driver}"
      COMMAND ${CMAKE_COMMAND} -E touch "${_stamp}"
      DEPENDS "${_last}"
      COMMENT "NumTracer: ${_name} done"
      VERBATIM)
    list(APPEND _outputs "${_stamp}")

    # The flow's completion as a nameable target: it is what the NEXT flow's generator hangs off,
    # and it doubles as a hand-runnable `make numtrace_<ns>` (which, per the header note, drags in
    # every pending flow ahead of it).
    add_custom_target(numtrace_${_ns} DEPENDS "${_stamp}")
    set_property(GLOBAL PROPERTY NUMTRACE_CHAIN_TAIL numtrace_${_ns})

    list(APPEND _pending "${_name}")
  endforeach()

  # `numtrace` is always defined, so `make numtrace` and the dependency below work whether or not
  # anything is pending. With every switch at 1 it is an empty target. A project with more than one
  # flow directory calls this function once per directory, so the second call hangs its work off a
  # per-directory sub-target rather than redefining `numtrace`.
  if(TARGET numtrace)
    string(MAKE_C_IDENTIFIER "numtrace_${flows_target}" _sub)
    add_custom_target(${_sub} DEPENDS ${_outputs})
    add_dependencies(numtrace ${_sub})
  else()
    add_custom_target(numtrace DEPENDS ${_outputs})
  endif()
  add_dependencies(${flows_target} numtrace)

  if(_pending)
    list(LENGTH _pending _n)
    string(REPLACE ";" " " _pretty "${_pending}")
    message(STATUS "numtrace: ${_n} flow(s) to generate, one at a time in this order: ${_pretty}")
  elseif(_manifests)
    message(STATUS "numtrace: all flows already generated")
  else()
    # No manifests at all: either the flows predate offline generation (they were emitted with
    # "Offline" -> False, so their kernels are already committed) or nothing has been emitted yet.
    # Either way there is nothing for this target to do.
    message(STATUS "numtrace: no numtrace manifests under ${flow_dir} — nothing to generate")
  endif()
endfunction()
