# Driver of the hex-ecs exchange-folding tests, invoked by CTest through "cmake -P".
#
# The option --fold-exchange drops the angular blocks with l1 > l2 from the solved basis
# and reconstructs them from their mirror counterparts. That is an exact identity, so the
# calculation has to give the same answer as one that solves for every block. Every test
# therefore runs the same problem twice, once with the option and once without it, and
# compares the partial cross sections that come out.
#
# The problem is deliberately tiny. What is being tested is that the two runs solve the
# same equations, which does not need a physically converged calculation, so the driver
# shrinks the sample input to something that takes seconds rather than hours. It does have
# to have a non-zero total angular momentum: for L = 0 and natural parity every angular
# state has l1 = l2, and the option would have nothing to fold.
#
# Required variables (passed via -D):
#   HEX_ECS   - path to the hex-ecs executable
#   WORKDIR   - directory to run in; erased and recreated first
#   L, PI, NL - total angular momentum, parity and number of angular levels
# Optional variables:
#   ARGS      - hex-ecs command line options common to both runs, semicolon-separated
#   LAUNCHER  - MPI launcher command, semicolon-separated
#   CHANNELS  - if true, add the projectile-exclusive knots to the input file
#   DIGITS    - significant digits of the cross sections that have to agree (default 4)
#   ULPS      - tolerance in units of the last of those digits (default 5)

foreach(var HEX_ECS WORKDIR L PI NL)
    if("${${var}}" STREQUAL "")
        message(FATAL_ERROR "run_fold_test.cmake: ${var} has not been set")
    endif()
endforeach()

set(INPUT_FILE "example.inp")

if("${DIGITS}" STREQUAL "")
    set(DIGITS 4)
endif()

if("${ULPS}" STREQUAL "")
    set(ULPS 5)
endif()

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

# Run one hex-ecs invocation in the given directory, keeping its output in logs beside
# the results. Aborts the test when the program fails.
function(run_hex_ecs dir stem)

    string(REPLACE ";" " " pretty "${ARGN}")
    message(STATUS "${pretty}")

    execute_process(
        COMMAND ${ARGN}
        WORKING_DIRECTORY "${dir}"
        OUTPUT_FILE "${dir}/${stem}.log"
        ERROR_FILE  "${dir}/${stem}.err"
        RESULT_VARIABLE result
    )

    if(NOT result EQUAL 0)
        tail_of_log("${dir}/${stem}.err")
        tail_of_log("${dir}/${stem}.log")
        message(FATAL_ERROR "${pretty}\nfailed with exit code ${result}")
    endif()

endfunction()

## --------------------------------------------------------------------------------- ##

# Replace the data lines that follow the comment line 'anchor' and are terminated by the
# line " -1" with the given text. This is how the knot sequences and the energy list of
# the sample input file are rewritten.
function(replace_section text_var anchor lines)

    set(text "${${text_var}}")

    string(FIND "${text}" "${anchor}" start)

    if(start EQUAL -1)
        message(FATAL_ERROR "Cannot patch ${INPUT_FILE}: it has no line\n  ${anchor}")
    endif()

    # stop just before the end of the anchor line, so that the terminator of a section
    # that has no data lines at all is found where it is and not further down the file
    string(LENGTH "${anchor}" anchor_length)
    math(EXPR start "${start} + ${anchor_length}")

    string(SUBSTRING "${text}" ${start} -1 tail)
    string(FIND "${tail}" "\n -1\n" stop)

    if(stop EQUAL -1)
        message(FATAL_ERROR "Cannot patch ${INPUT_FILE}: the section below\n  ${anchor}\nis not terminated by -1")
    endif()

    math(EXPR stop "${start} + ${stop} + 1")

    string(SUBSTRING "${text}" 0 ${start} head)
    string(SUBSTRING "${text}" ${stop} -1 tail)

    set(${text_var} "${head}\n${lines}${tail}" PARENT_SCOPE)

endfunction()

# Replace the single data line that follows the comment line 'anchor'.
function(replace_line text_var anchor line)

    set(text "${${text_var}}")

    string(FIND "${text}" "${anchor}" start)

    if(start EQUAL -1)
        message(FATAL_ERROR "Cannot patch ${INPUT_FILE}: it has no line\n  ${anchor}")
    endif()

    string(LENGTH "${anchor}" anchor_length)
    math(EXPR start "${start} + ${anchor_length} + 1")

    string(SUBSTRING "${text}" ${start} -1 tail)
    string(FIND "${tail}" "\n" stop)

    math(EXPR stop "${start} + ${stop} + 1")

    string(SUBSTRING "${text}" 0 ${start} head)
    string(SUBSTRING "${text}" ${stop} -1 tail)

    set(${text_var} "${head}${line}\n${tail}" PARENT_SCOPE)

endfunction()

## --------------------------------------------------------------------------------- ##

# Split a decimal number into a sign, a string of significant digits and a decimal
# exponent, so that the value is <sign> 0.<digits> * 10^<exponent>. Everything is done on
# the strings, because CMake can only do integer arithmetic and the numbers have to be
# compared with a tolerance.
function(split_number text sign_var digits_var exponent_var)

    string(STRIP "${text}" text)

    # sign
    set(sign "+")
    if(text MATCHES "^-")
        set(sign "-")
    endif()
    string(REGEX REPLACE "^[-+]" "" text "${text}")

    # explicit exponent, if any
    set(exponent 0)
    if(text MATCHES "[eE]")
        string(REGEX REPLACE "^[^eE]*[eE]([-+]?[0-9]+)$" "\\1" exponent "${text}")
        string(REGEX REPLACE "[eE].*$" "" text "${text}")
        string(REGEX REPLACE "^\\+" "" exponent "${exponent}")
    endif()

    if(NOT text MATCHES "^[0-9]*\\.?[0-9]*$" OR text STREQUAL "" OR text STREQUAL ".")
        message(FATAL_ERROR "Cannot read \"${text}\" as a number.")
    endif()

    # digits before and after the decimal point
    set(fraction "")
    if(text MATCHES "\\.")
        string(REGEX REPLACE "^([0-9]*)\\..*$" "\\1" whole "${text}")
        string(REGEX REPLACE "^[0-9]*\\.(.*)$" "\\1" fraction "${text}")
    else()
        set(whole "${text}")
    endif()

    # the value is 0.<whole><fraction> * 10^(exponent + number of digits in <whole>)
    string(LENGTH "${whole}" length)
    math(EXPR exponent "${exponent} + ${length}")
    set(digits "${whole}${fraction}")

    # normalize, so that the first digit is significant
    while(digits MATCHES "^0")
        string(SUBSTRING "${digits}" 1 -1 digits)
        math(EXPR exponent "${exponent} - 1")
    endwhile()

    # zero has no significant digits and no meaningful exponent
    if(digits STREQUAL "")
        set(digits "0")
        set(exponent 0)
    endif()

    set(${sign_var} "${sign}" PARENT_SCOPE)
    set(${digits_var} "${digits}" PARENT_SCOPE)
    set(${exponent_var} "${exponent}" PARENT_SCOPE)

endfunction()

# Cut or pad a string of digits to the requested length and return it as an integer.
function(digits_to_integer digits length result_var)

    string(LENGTH "${digits}" have)

    while(have LESS length)
        set(digits "${digits}0")
        math(EXPR have "${have} + 1")
    endwhile()

    string(SUBSTRING "${digits}" 0 ${length} digits)

    # a leading zero would make CMake read the number as octal
    string(REGEX REPLACE "^0+" "" digits "${digits}")

    if(digits STREQUAL "")
        set(digits "0")
    endif()

    set(${result_var} "${digits}" PARENT_SCOPE)

endfunction()

# Whether two decimal numbers agree in their first DIGITS significant digits, up to a
# tolerance of ULPS units of the last of them.
function(numbers_agree a b result_var)

    split_number("${a}" sign_a digits_a exponent_a)
    split_number("${b}" sign_b digits_b exponent_b)

    set(${result_var} FALSE PARENT_SCOPE)

    if(NOT sign_a STREQUAL sign_b)
        return()
    endif()

    # bring both to a common exponent; they can differ by one when the values sit on a
    # rounding boundary, but anything more is a real disagreement
    if(exponent_a LESS exponent_b)
        set(spread_a 1)
        set(spread_b 0)
        math(EXPR gap "${exponent_b} - ${exponent_a}")
    else()
        set(spread_a 0)
        set(spread_b 1)
        math(EXPR gap "${exponent_a} - ${exponent_b}")
    endif()

    if(gap GREATER 1)
        return()
    endif()

    if(gap EQUAL 0)
        set(spread_a 0)
        set(spread_b 0)
        set(tolerance ${ULPS})
    else()
        math(EXPR tolerance "10 * ${ULPS}")
    endif()

    math(EXPR length "${DIGITS} + ${gap}")

    string(REPEAT "0" ${spread_a} pad_a)
    string(REPEAT "0" ${spread_b} pad_b)

    digits_to_integer("${pad_a}${digits_a}" ${length} number_a)
    digits_to_integer("${pad_b}${digits_b}" ${length} number_b)

    if(number_a LESS number_b)
        math(EXPR difference "${number_b} - ${number_a}")
    else()
        math(EXPR difference "${number_a} - ${number_b}")
    endif()

    if(difference LESS_EQUAL tolerance)
        set(${result_var} TRUE PARENT_SCOPE)
    endif()

endfunction()

# Compare the numbers in the data (non-comment) lines of two cross section files.
function(compare_results reference folded)

    foreach(path "${reference}" "${folded}")
        if(NOT EXISTS "${path}")
            message(FATAL_ERROR "The solver did not produce ${path}")
        endif()
    endforeach()

    file(STRINGS "${reference}" lines_reference REGEX "^[^#]")
    file(STRINGS "${folded}"    lines_folded    REGEX "^[^#]")

    list(LENGTH lines_reference count)
    list(LENGTH lines_folded count_folded)

    if(count EQUAL 0)
        message(FATAL_ERROR "There is no data in ${reference}")
    endif()

    # A partial wave that cannot be reached from the initial state gives files that hold
    # nothing but the energies. Comparing those would prove nothing at all.
    list(GET lines_reference 0 first)
    string(REGEX MATCHALL "[^ \t]+" fields "${first}")
    list(LENGTH fields fields)

    if(fields LESS 2)
        message(FATAL_ERROR
            "${reference} contains no cross sections, only the energies. This problem "
            "has no transition allowed by its total quantum numbers, so the test would "
            "prove nothing")
    endif()

    if(NOT count EQUAL count_folded)
        message(FATAL_ERROR "${reference} has ${count} data lines, but ${folded} has ${count_folded}")
    endif()

    math(EXPR last "${count} - 1")

    foreach(i RANGE ${last})

        list(GET lines_reference ${i} line_reference)
        list(GET lines_folded    ${i} line_folded)

        string(REGEX MATCHALL "[^ \t]+" fields_reference "${line_reference}")
        string(REGEX MATCHALL "[^ \t]+" fields_folded    "${line_folded}")

        list(LENGTH fields_reference fields)
        list(LENGTH fields_folded fields_other)

        if(NOT fields EQUAL fields_other)
            message(FATAL_ERROR
                "Different number of columns on line ${i}:\n"
                "  reference: ${line_reference}\n"
                "  folded:    ${line_folded}")
        endif()

        math(EXPR last_field "${fields} - 1")

        foreach(j RANGE ${last_field})

            list(GET fields_reference ${j} a)
            list(GET fields_folded    ${j} b)

            numbers_agree("${a}" "${b}" agree)

            if(NOT agree)
                message(FATAL_ERROR
                    "The folded calculation gave a different result in column ${j} of "
                    "line ${i} of ${reference}:\n"
                    "  reference: ${a}\n"
                    "  folded:    ${b}\n"
                    "(compared to ${DIGITS} significant digits with a tolerance of ${ULPS} units of the last one)")
            endif()

        endforeach()

    endforeach()

    message(STATUS "${count} line(s) of ${reference} reproduced by the folded calculation")

endfunction()

## --------------------------------------------------------------------------------- ##

# start from scratch; the tests reuse cached radial integrals otherwise
file(REMOVE_RECURSE "${WORKDIR}")
file(MAKE_DIRECTORY "${WORKDIR}")

run_hex_ecs("${WORKDIR}" example "${HEX_ECS}" --example)

if(NOT EXISTS "${WORKDIR}/${INPUT_FILE}")
    message(FATAL_ERROR "hex-ecs --example did not write ${INPUT_FILE}")
endif()

file(READ "${WORKDIR}/${INPUT_FILE}" input)

# A grid that is far too coarse for the physics but keeps the runs down to seconds. Both
# runs discretize the problem in exactly the same way, so this does not weaken the test.
replace_section(input
    "# a) Real knots of the basis that is common to atomic and projectile electron."
    "  L  0.0  0.0   4\n  G  0.1  6.0  0.3  1.2\n  L    7   20   14\n")

replace_section(input
    "# c) Complex region knots. (Start from zero.)"
    "  L    0   10   11\n")

# Without knots of its own for the projectile this is an inner-region-only problem. The
# channel variant adds them, which brings in the asymptotic channels: those are stored
# per escaping electron, so mirroring a solution segment has to swap them.
if(CHANNELS)
    replace_section(input
        "# b) Real knots that are exclusive to the projectile, if any. (Start from zero.)"
        "  L    0   12   13\n")

    # keep a few closed channels as well, so that the two groups differ in size
    replace_line(input
        "# Maximal energy (Ry) of states included in the asymptotic (outer) region."
        "  -0.1")
endif()

replace_line(input
    "# L   S   Pi  nL  limit exchange"
    "  ${L}   *   ${PI}   ${NL}   -1    1")

replace_section(input
    "#     E[xplicit] <Sample-1> <Sample-2> ... <Sample-last> -1"
    "  L  -0.35  -0.35  1\n")

file(WRITE "${WORKDIR}/${INPUT_FILE}" "${input}")

## --------------------------------------------------------------------------------- ##

foreach(run reference folded)

    file(MAKE_DIRECTORY "${WORKDIR}/${run}")
    file(COPY "${WORKDIR}/${INPUT_FILE}" DESTINATION "${WORKDIR}/${run}")

    set(options ${ARGS})
    if(run STREQUAL "folded")
        list(APPEND options --fold-exchange)
    endif()

    run_hex_ecs("${WORKDIR}/${run}" run ${LAUNCHER} "${HEX_ECS}" --input "${INPUT_FILE}" ${options})

endforeach()

## --------------------------------------------------------------------------------- ##

# Make sure the option had something to do. If the angular basis were not reduced, the two
# runs would agree for the trivial reason that they are the same calculation.
file(STRINGS "${WORKDIR}/folded/run.log" folding REGEX "exchange symmetry folded in")

if(NOT folding)
    tail_of_log("${WORKDIR}/folded/run.log")
    message(FATAL_ERROR "The run with --fold-exchange did not fold the angular basis")
endif()

string(REGEX MATCH "([0-9]+) of ([0-9]+) blocks" ignored "${folding}")
set(solved ${CMAKE_MATCH_1})
set(all ${CMAKE_MATCH_2})

if(NOT solved LESS all)
    message(FATAL_ERROR
        "The angular basis of this problem has no blocks to fold "
        "(${solved} of ${all} blocks are solved for), so the test would prove nothing")
endif()

message(STATUS "The folded run solved for ${solved} of ${all} angular blocks")

## --------------------------------------------------------------------------------- ##

foreach(spin 0 1)
    set(output "ics-L${L}-S${spin}-Pi${PI}.dat")
    compare_results("${WORKDIR}/reference/${output}" "${WORKDIR}/folded/${output}")
endforeach()

message(STATUS "Test passed")
