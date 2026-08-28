# Driver of the hex-ecs regression tests, invoked by CTest through "cmake -P".
# Every test gets a fresh directory in the build tree, where the driver
#
#   1. lets hex-ecs write its sample input file,
#   2. optionally adds the projectile-exclusive knots needed by channel tests,
#   3. runs hex-ecs with the options of this test,
#   4. checks that the expected cross section files have been written.
#
# Required variables (passed via -D):
#   HEX_ECS   - path to the hex-ecs executable
#   WORKDIR   - directory to run in; erased and recreated first
#   OUTPUTS   - files that the solver has to produce, relative to WORKDIR
# Optional variables:
#   ARGS      - hex-ecs command line options, semicolon-separated
#   LAUNCHER  - MPI launcher command, semicolon-separated
#   CHANNELS  - if true, add the projectile-exclusive knots to the input file

foreach(var HEX_ECS WORKDIR OUTPUTS)
    if("${${var}}" STREQUAL "")
        message(FATAL_ERROR "run_test.cmake: ${var} has not been set")
    endif()
endforeach()

set(INPUT_FILE "example.inp")

## --------------------------------------------------------------------------------- ##

# Print the tail of a log file, so that a failure reported by CTest carries the
# reason with it and not just a path to look at.
function(tail_of_log path)

    if(NOT EXISTS "${path}")
        return()
    endif()

    file(READ "${path}" text)
    string(LENGTH "${text}" length)

    if(length GREATER 4000)
        math(EXPR offset "${length} - 4000")
        string(SUBSTRING "${text}" ${offset} -1 text)
        set(text "[...]\n${text}")
    endif()

    if(NOT text STREQUAL "")
        message("--- ${path} ---\n${text}")
    endif()

endfunction()

# Run one hex-ecs invocation in the working directory, keeping its output in
# logs beside the results. Aborts the test when the program fails.
function(run_hex_ecs stem)

    string(REPLACE ";" " " pretty "${ARGN}")
    message(STATUS "${pretty}")

    execute_process(
        COMMAND ${ARGN}
        WORKING_DIRECTORY "${WORKDIR}"
        OUTPUT_FILE "${WORKDIR}/${stem}.log"
        ERROR_FILE  "${WORKDIR}/${stem}.err"
        RESULT_VARIABLE result
    )

    if(NOT result EQUAL 0)
        tail_of_log("${WORKDIR}/${stem}.err")
        tail_of_log("${WORKDIR}/${stem}.log")
        message(FATAL_ERROR "${pretty}\nfailed with exit code ${result}")
    endif()

endfunction()

## --------------------------------------------------------------------------------- ##

# start from scratch; the tests reuse cached radial integrals otherwise
file(REMOVE_RECURSE "${WORKDIR}")
file(MAKE_DIRECTORY "${WORKDIR}")

run_hex_ecs(example "${HEX_ECS}" --example)

if(NOT EXISTS "${WORKDIR}/${INPUT_FILE}")
    message(FATAL_ERROR "hex-ecs --example did not write ${INPUT_FILE}")
endif()

# The sample input has no knots of its own for the projectile, which makes it an
# inner-region-only problem. Channel tests add such knots to solve the full one.
if(CHANNELS)

    set(ANCHOR "# b) Real knots that are exclusive to the projectile, if any. (Start from zero.)")

    file(READ "${WORKDIR}/${INPUT_FILE}" input)

    # the anchor is a plain comment line, so compare it literally rather than as a regex
    string(FIND "${input}" "${ANCHOR}\n" anchor_pos)

    if(anchor_pos EQUAL -1)
        message(FATAL_ERROR "Cannot add projectile knots: ${INPUT_FILE} has no line\n  ${ANCHOR}")
    endif()

    string(REPLACE "${ANCHOR}\n" "${ANCHOR}\n  L    0  100  101\n" input "${input}")

    file(WRITE "${WORKDIR}/${INPUT_FILE}" "${input}")

endif()

run_hex_ecs(run ${LAUNCHER} "${HEX_ECS}" --input "${INPUT_FILE}" ${ARGS})

## --------------------------------------------------------------------------------- ##

foreach(output IN LISTS OUTPUTS)

    if(NOT EXISTS "${WORKDIR}/${output}")
        tail_of_log("${WORKDIR}/run.log")
        message(FATAL_ERROR "The solver did not produce ${output}")
    endif()

    file(READ "${WORKDIR}/${output}" head LIMIT 1)

    if(head STREQUAL "")
        message(FATAL_ERROR "The solver produced an empty ${output}")
    endif()

endforeach()

message(STATUS "Test passed")
