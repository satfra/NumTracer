# Run a tutorial program and compare its stdout with the expected-output file the documentation
# quotes, so a page can never show output the program no longer prints.
#   cmake -DPROG=<exe> -DEXPECTED=<file> -P compare_output.cmake
execute_process(COMMAND ${PROG} OUTPUT_VARIABLE out RESULT_VARIABLE rc)
file(READ ${EXPECTED} want)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "${PROG} exited with ${rc}")
endif()
if(NOT out STREQUAL want)
  message(FATAL_ERROR "output of ${PROG} differs from ${EXPECTED} (the docs quote that file).\n"
                      "--- got ---\n${out}--- expected ---\n${want}\n"
                      "If the change is intended, regenerate the file: ${PROG} > ${EXPECTED}")
endif()
