# NumTracerNumtraceRun.cmake — build-step driver for the `numtrace` target (cmake -P).
#
# Two modes, selected with -DMODE=:
#
#   run   run a flow's generator binary and commit its stdout as the straight-line traces header,
#         with its thread pool sized to the build's own job budget (see "workers" below).
#         -DGEN=<binary> -DOUT=<kernels.hh> -DNS=<namespace> -DDECOR=<decorator>
#         -DFLOW=<name> -DJOBS=<n> -DMAXW=<n|0> -DMAXWB=<n|0> -DIDX=<k> -DTOTAL=<n>
#
#   probe run a flow's imaginary-part probe, which writes the verdict header, and report the verdict
#         in words. The probe's raw measurement line is captured rather than dumped into the build log.
#         -DPROBE=<binary> -DOUT=<numtrace_verdict.hh> -DMACRO=<name> -DFLOW=<name>
#
#   mark  flip a flow's numtrace.json switch to "generated": 1.
#         -DMANIFEST=<numtrace.json>
#
# `run` mirrors the safety checks the Wolfram path used to do inline: the generator writes to a temp
# file, which is validated (exit code 0 AND non-trivial size) before replacing the committed header —
# a crashed generator must never truncate a header that is checked into git.
#
# workers. The generator sizes its pool from hardware_concurrency() and only ever lets NT_GEN_MAXW
# lower it (NT_GEN_MAXW_B can also raise it), so the whole thread count is ours to set here. It is
# resolved fresh per flow, in this order:
#
#   1. an NT_GEN_MAXW / NT_GEN_MAXW_B already in the environment wins outright, so
#      `NT_GEN_MAXW=2 make -j8 numtrace` still works for debugging;
#   2. otherwise the build's own -jN, taken from MAKEFLAGS, which GNU make exports into every recipe
#      (measured on make 4.4.1: " -j8 --jobserver-auth=fifo:/tmp/GMfifo666907"). This is the only
#      place the number exists — it is chosen at build time, not at configure time;
#   3. failing that -DJOBS=, the NUMTRACE_JOBS cache value, for tools that export nothing (Ninja);
#   4. capped by -DMAXW=/-DMAXWB= when the flow's manifest carries one, i.e. when the notebook had
#      SetNumTracerThreads in force while emitting it. A cap only ever lowers: a manifest asking for
#      4 does not get 4 out of a `-j2` build.

cmake_minimum_required(VERSION 3.19)

if(MODE STREQUAL "run")
  foreach(v GEN OUT NS DECOR)
    if(NOT DEFINED ${v})
      message(FATAL_ERROR "numtrace run: -D${v}= is required")
    endif()
  endforeach()

  set(_tmp "${OUT}.tmp")
  set(_args -n "${NS}" -d "${DECOR}")

  # ---- worker count (see the header note) -------------------------------------------------------
  set(_jobs "${JOBS}")
  if(DEFINED ENV{MAKEFLAGS})
    set(_mf "$ENV{MAKEFLAGS}")
    # --jobs= is tried first so that --jobserver-auth=fifo:/tmp/GMfifo665822 cannot be misread; -j
    # then needs a digit right after it, which "-jobserver" does not have. A bare -j (unlimited)
    # matches neither and falls through to JOBS.
    if(_mf MATCHES "--jobs[= ]+([0-9]+)")
      set(_jobs "${CMAKE_MATCH_1}")
    elseif(_mf MATCHES "-j *([0-9]+)")
      set(_jobs "${CMAKE_MATCH_1}")
    endif()
  endif()
  if(NOT _jobs OR _jobs LESS 1)
    set(_jobs 1)
  endif()

  set(_w "${_jobs}")
  if(MAXW AND MAXW GREATER 0 AND MAXW LESS _w)
    set(_w "${MAXW}")
  endif()
  set(_wb "${_jobs}")
  if(MAXWB AND MAXWB GREATER 0 AND MAXWB LESS _wb)
    set(_wb "${MAXWB}")
  endif()

  # set(ENV{}) here is inherited by execute_process's child, so no `cmake -E env` wrapper is needed.
  if(DEFINED ENV{NT_GEN_MAXW})
    set(_w "$ENV{NT_GEN_MAXW} (from the environment)")
  else()
    set(ENV{NT_GEN_MAXW} "${_w}")
  endif()
  if(DEFINED ENV{NT_GEN_MAXW_B})
    set(_wb "$ENV{NT_GEN_MAXW_B} (from the environment)")
  else()
    set(ENV{NT_GEN_MAXW_B} "${_wb}")
  endif()

  # ---- device target ----------------------------------------------------------------------------
  # The offline twin of CodegenBuild.m's ntDeviceEnvPrefix. It enables gen.hpp's size-gated `__noinline__`,
  # which is device-only (the host has no register cliff and its emission stays byte-identical). This
  # is the whole reason the manifest carries a "device" field: nothing of the emitting Wolfram
  # kernel's environment survives the trip into a `cmake -P` build step, so a fact not written down
  # is lost -- and this one was, on every offline-generated flow.
  set(_dev "host")
  if(DEVICE)
    if(DEFINED ENV{NT_GEN_DEVICE})
      set(_dev "device (NT_GEN_DEVICE from the environment)")
    else()
      set(ENV{NT_GEN_DEVICE} "1")
      set(_dev "device")
    endif()
  endif()

  # With the flows serialized this banner is the only thing saying which one the build is sitting on
  # and how much of the machine it was given.
  set(_which "")
  if(TOTAL AND TOTAL GREATER 0)
    set(_which "[${IDX}/${TOTAL}] ")
  endif()
  if(NOT DEFINED FLOW OR FLOW STREQUAL "")
    set(FLOW "${NS}")
  endif()
  message(STATUS "NumTracer: ${_which}tracing ${FLOW} (W=${_w}, WB=${_wb}, ${_dev})")

  # ECHO_ERROR_VARIABLE tees the generator's stderr to the console as it arrives — the
  # NT_GEN_PROFILE phase lines — while still capturing it for the failure message below.
  execute_process(
    COMMAND "${GEN}" ${_args}
    OUTPUT_FILE "${_tmp}"
    ERROR_VARIABLE _err
    ECHO_ERROR_VARIABLE
    RESULT_VARIABLE _rc)

  set(_sz 0)
  if(EXISTS "${_tmp}")
    file(SIZE "${_tmp}" _sz)
  endif()

  if(NOT _rc EQUAL 0 OR _sz LESS 64)
    file(REMOVE "${_tmp}")
    message(FATAL_ERROR
      "numtrace: generator for '${NS}' failed (exit ${_rc}, ${_sz} bytes).\n"
      "The committed ${OUT} was left intact.\n${_err}")
  endif()

  file(RENAME "${_tmp}" "${OUT}")
  message(STATUS "NumTracer: wrote ${OUT} (${_sz} bytes)")

elseif(MODE STREQUAL "probe")
  foreach(v PROBE OUT MACRO FLOW)
    if(NOT DEFINED ${v})
      message(FATAL_ERROR "numtrace probe: -D${v}= is required")
    endif()
  endforeach()

  execute_process(
    COMMAND "${PROBE}" -o "${OUT}" -m "${MACRO}"
    OUTPUT_VARIABLE _out
    ERROR_VARIABLE _err
    RESULT_VARIABLE _rc)

  # Any failure is fatal: the verdict picks a preprocessor branch in a committed header, so a
  # probe that did not actually run must never be mistaken for "it came out complex".
  if(NOT _rc EQUAL 0 OR NOT EXISTS "${OUT}")
    message(FATAL_ERROR
      "numtrace: imaginary-part probe for '${FLOW}' failed (exit ${_rc}).\n${_err}${_out}")
  endif()

  # stdout is "mim mdiff mre mrim mrdiff mrrep nrep ok verdict"; report the verdict in words and keep
  # the relative imaginary part, which is the number worth seeing when a flow unexpectedly stays
  # complex.
  #
  # Indices are taken from the END, not the front. This block used to test `_nf EQUAL 7` against a
  # layout that had since grown to 9 fields, so EVERY offline build silently fell through to the bare
  # "probed" message and no verdict was ever reported — the failure mode of a hard-coded count is that
  # it goes quiet rather than complaining. verdict/ok are appended last and mrim sits in a stable
  # prefix, so this survives the next field being added too.
  # REGEX on runs of whitespace, not string(REPLACE " " ";"): the latter turns a double space into an
  # empty list element, which keeps LENGTH plausible while shifting every index past it.
  string(STRIP "${_out}" _out)
  string(REGEX REPLACE "[ \t]+" ";" _fields "${_out}")
  list(LENGTH _fields _nf)
  if(_nf GREATER_EQUAL 9)
    math(EXPR _iv "${_nf} - 1")
    math(EXPR _ipts "${_nf} - 2")
    list(GET _fields 3 _relim)
    list(GET _fields ${_ipts} _npts)
    list(GET _fields ${_iv} _v)
    if(_v EQUAL 2)
      set(_word "Pure (imaginary coefficients dropped)")
    elseif(_v EQUAL 1)
      set(_word "RePart (re/im split)")
    else()
      set(_word "complex")
    endif()
    if(_v EQUAL 0)
      # Not a STATUS line: a DiFfRG flow consumer accumulates into a `double`, so a complex kernel
      # body does not merely waste a lowered body — it fails to compile, far downstream, as a concept
      # error naming the integrator rather than this flow. Say so here, where the cause is known.
      message(WARNING
        "NumTracer: ${FLOW} probed ${_word} (rel|Im|=${_relim} over ${_npts} pts). A real-valued "
        "consumer cannot instantiate the complex kernel body. Either give this flow a complex "
        "integrator, or regenerate it with \"RealOutput\" -> True to take Re[...] deliberately.")
    else()
      message(STATUS "NumTracer: ${FLOW} verdict = ${_word}, rel|Im|=${_relim} over ${_npts} pts")
    endif()
  else()
    message(STATUS "NumTracer: ${FLOW} probed (unrecognised probe output: ${_nf} field(s))")
  endif()

elseif(MODE STREQUAL "mark")
  if(NOT DEFINED MANIFEST)
    message(FATAL_ERROR "numtrace mark: -DMANIFEST= is required")
  endif()
  file(READ "${MANIFEST}" _json)
  string(JSON _json SET "${_json}" "generated" 1)
  file(WRITE "${MANIFEST}" "${_json}")
  message(STATUS "NumTracer: ${MANIFEST} -> generated: 1")

else()
  message(FATAL_ERROR "numtrace: unknown -DMODE=${MODE} (expected 'run', 'probe' or 'mark')")
endif()
