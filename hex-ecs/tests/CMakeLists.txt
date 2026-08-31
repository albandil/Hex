## --------------------------------------------------------------------------------- ##
##                                                                                   ##
##                       / /   / /    __    \ \  / /                                 ##
##                      / /__ / /   / _ \    \ \/ /                                  ##
##                     /  ___  /   | |/_/    / /\ \                                  ##
##                    / /   / /    \_\      / /  \ \                                 ##
##                                                                                   ##
##                                                                                   ##
##  Copyright (c) 2019, Jakub Benda, Charles University in Prague                    ##
##                                                                                   ##
## MIT License:                                                                      ##
##                                                                                   ##
##  Permission is hereby granted, free of charge, to any person obtaining a          ##
## copy of this software and associated documentation files (the "Software"),        ##
## to deal in the Software without restriction, including without limitation         ##
## the rights to use, copy, modify, merge, publish, distribute, sublicense,          ##
## and/or sell copies of the Software, and to permit persons to whom the             ##
## Software is furnished to do so, subject to the following conditions:              ##
##                                                                                   ##
##  The above copyright notice and this permission notice shall be included          ##
## in all copies or substantial portions of the Software.                            ##
##                                                                                   ##
##  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS          ##
## OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       ##
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE       ##
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, ##
## WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF         ##
## OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.  ##
##                                                                                   ##
## --------------------------------------------------------------------------------- ##

# Regression tests of hex-ecs, of two kinds.
#
# Those added by hex_ecs_add_test solve the sample problem written by
# "hex-ecs --example" using a different combination of preconditioners, matrix
# factorizers and parallelization modes, and check that the integral cross
# sections come out. They are full physical calculations, not unit tests, so
# expect minutes to hours per test rather than seconds.
#
# Those added by hex_ecs_add_fold_test check the exchange symmetry folding of the
# angular basis (the --fold-exchange option) by solving a much smaller problem twice,
# with and without the option, and comparing the results. Those run in seconds.
#
# Run them with
#     ctest --test-dir <build directory> --output-on-failure
# Tests that need a library the build was configured without are not created;
# the configuration summary below lists what has been left out.
#
# The distributed tests are started by the launcher that FindMPI has detected.
# Anything the launcher needs beyond the number of processes belongs to the cache
# variable MPIEXEC_PREFLAGS; with Open MPI the tests that combine several ranks
# with several threads want at least "--bind-to none" there. Note the two dashes:
# Open MPI 5 takes "-bind-to none" for a request to run a program called "none".

# Four hours. The KPA test of the sample problem needs close to two on a six-core
# machine, and the slowest energy of that run converges erratically, so the limit
# is meant as a bound on a hung test and not as a schedule.
set(HEX_ECS_TEST_TIMEOUT 14400 CACHE STRING "Time limit of a single hex-ecs test [s]")

# both spin channels of the lowest partial wave of the sample problem
set(HEX_ECS_TEST_OUTPUTS "ics-L0-S0-Pi0.dat" "ics-L0-S1-Pi0.dat")

set(hex-ecs_TESTS_SKIPPED "")

## --------------------------------------------------------------------------------- ##

# hex_ecs_add_test(<name>
#     [ARGS <hex-ecs options>...]    options that distinguish this test
#     [REQUIRES <feature>...]        one or more of umfpack, mumps, pardiso, opencl
#     [CHANNELS]                     extend the input file by projectile-only knots
#     [RANKS <n>]                    run under the MPI launcher with <n> processes
#     [THREADS <n>]                  value of OMP_NUM_THREADS for the test
# )
function(hex_ecs_add_test name)

    cmake_parse_arguments(TEST "CHANNELS" "RANKS;THREADS" "ARGS;REQUIRES" ${ARGN})

    if(TEST_UNPARSED_ARGUMENTS)
        message(FATAL_ERROR "hex_ecs_add_test(${name}): unknown arguments ${TEST_UNPARSED_ARGUMENTS}")
    endif()

    # an MPI test is only meaningful when there is something to launch it with
    if(TEST_RANKS)
        list(APPEND TEST_REQUIRES mpi)
    endif()

    # skip silently what this build cannot do; the missing features are summarized
    # once at the end of the configuration
    set(missing "")

    foreach(feature IN LISTS TEST_REQUIRES)
        if("${feature}" STREQUAL "umfpack" AND NOT WITH_UMFPACK)
            list(APPEND missing umfpack)
        elseif("${feature}" STREQUAL "mumps" AND NOT WITH_MUMPS)
            list(APPEND missing mumps)
        elseif("${feature}" STREQUAL "pardiso" AND NOT WITH_PARDISO AND NOT WITH_MKL)
            list(APPEND missing pardiso)
        elseif("${feature}" STREQUAL "opencl" AND NOT WITH_OPENCL)
            list(APPEND missing opencl)
        elseif("${feature}" STREQUAL "mpi" AND NOT MPIEXEC_EXECUTABLE)
            list(APPEND missing mpi)
        endif()
    endforeach()

    if(missing)
        set(hex-ecs_TESTS_SKIPPED ${hex-ecs_TESTS_SKIPPED} ${missing} PARENT_SCOPE)
        return()
    endif()

    ## ----------------------------------------------------------------------------- ##

    set(launcher "")

    if(TEST_RANKS)
        set(launcher
            "${MPIEXEC_EXECUTABLE}"
            ${MPIEXEC_NUMPROC_FLAG} ${TEST_RANKS}
            ${MPIEXEC_PREFLAGS}
        )
        list(APPEND TEST_ARGS --mpi)
    endif()

    add_test(
        NAME "hex-ecs.${name}"
        COMMAND ${CMAKE_COMMAND}
                -D "HEX_ECS=$<TARGET_FILE:hex-ecs>"
                -D "WORKDIR=${CMAKE_CURRENT_BINARY_DIR}/${name}"
                -D "OUTPUTS=${HEX_ECS_TEST_OUTPUTS}"
                -D "ARGS=${TEST_ARGS}"
                -D "LAUNCHER=${launcher}"
                -D "CHANNELS=${TEST_CHANNELS}"
                -P "${CMAKE_CURRENT_LIST_DIR}/cmake/run_test.cmake"
    )

    # every test saturates the machine with OpenMP threads or MPI ranks, so there is
    # nothing to gain from running several of them at the same time
    set_tests_properties("hex-ecs.${name}" PROPERTIES
        TIMEOUT     ${HEX_ECS_TEST_TIMEOUT}
        RUN_SERIAL  ON
    )

    if(TEST_THREADS)
        set_tests_properties("hex-ecs.${name}" PROPERTIES
            ENVIRONMENT "OMP_NUM_THREADS=${TEST_THREADS}"
        )
    endif()

endfunction()

## --------------------------------------------------------------------------------- ##

# inner region only, varying the preconditioner and the factorizer

hex_ecs_add_test(01-inner
    REQUIRES umfpack)

hex_ecs_add_test(02-inner-ILU
    ARGS --preconditioner ILU
    REQUIRES umfpack)

hex_ecs_add_test(03-inner-ILU-mumps
    ARGS --preconditioner ILU --lu mumps
    REQUIRES mumps
    RANKS 1)

hex_ecs_add_test(04-inner-ILU-pardiso
    ARGS --preconditioner ILU --lu pardiso
    REQUIRES pardiso)

hex_ecs_add_test(05-inner-ILU-lightweight
    ARGS --preconditioner ILU --lightweight-full
    REQUIRES umfpack)

hex_ecs_add_test(06-inner-KPA
    ARGS --preconditioner KPA)

hex_ecs_add_test(07-inner-GPU
    ARGS --preconditioner GPU
    REQUIRES opencl)

## --------------------------------------------------------------------------------- ##

# full problem including the projectile channels

hex_ecs_add_test(08-channels
    CHANNELS
    REQUIRES umfpack)

hex_ecs_add_test(09-channels-mumps
    ARGS --lu mumps
    CHANNELS
    REQUIRES mumps
    RANKS 1)

hex_ecs_add_test(10-channels-pardiso
    ARGS --lu pardiso
    CHANNELS
    REQUIRES pardiso)

hex_ecs_add_test(11-channels-lightweight
    ARGS --lightweight-full
    CHANNELS
    REQUIRES umfpack)

## --------------------------------------------------------------------------------- ##

# domain decomposition

hex_ecs_add_test(12-inner-DOM-2x2-ILU
    ARGS --preconditioner DOM --dom-xpanels 2 --dom-ypanels 2
         --dom-preconditioner ILU --tolerance 0.1
    REQUIRES umfpack)

hex_ecs_add_test(13-inner-DOM-2x2-KPA
    ARGS --preconditioner DOM --dom-xpanels 2 --dom-ypanels 2
         --dom-preconditioner KPA --tolerance 0.1)

hex_ecs_add_test(14-inner-DOM-1x2-ILU
    ARGS --preconditioner DOM --dom-xpanels 1 --dom-ypanels 2
         --dom-preconditioner ILU --tolerance 0.1
    REQUIRES umfpack)

hex_ecs_add_test(15-inner-DOM-2x1-ILU
    ARGS --preconditioner DOM --dom-xpanels 2 --dom-ypanels 1
         --dom-preconditioner ILU --tolerance 0.1
    REQUIRES umfpack)

hex_ecs_add_test(16-inner-DOM-2x2-KPA-lightweight
    ARGS --preconditioner DOM --dom-xpanels 2 --dom-ypanels 2
         --dom-preconditioner KPA --lightweight-full --tolerance 0.1)

## --------------------------------------------------------------------------------- ##

# out-of-core storage of the matrices

hex_ecs_add_test(17-inner-ILU-out-of-core
    ARGS --preconditioner ILU --out-of-core
    REQUIRES umfpack)

hex_ecs_add_test(18-inner-ILU-out-of-core-lightweight
    ARGS --preconditioner ILU --out-of-core --lightweight-full
    REQUIRES umfpack)

hex_ecs_add_test(19-channels-lightweight-analytic-eigenstates
    ARGS --lightweight-full --analytic-eigenstates
    CHANNELS
    REQUIRES umfpack)

## --------------------------------------------------------------------------------- ##

# distributed runs

hex_ecs_add_test(20-inner-ILU-mpi
    ARGS --shared-scratch --preconditioner ILU
    REQUIRES umfpack
    RANKS 2 THREADS 2)

hex_ecs_add_test(21-inner-ILU-mpi-lightweight
    ARGS --shared-scratch --preconditioner ILU --lightweight-full
    REQUIRES umfpack
    RANKS 2 THREADS 2)

hex_ecs_add_test(22-inner-ILU-mpi-distributed-mumps
    ARGS --groupsize 2 --shared-scratch --preconditioner ILU --lu mumps
    REQUIRES mumps
    RANKS 4 THREADS 1)

hex_ecs_add_test(23-inner-ILU-mpi-out-of-core
    ARGS --shared-scratch --preconditioner ILU --out-of-core
    REQUIRES umfpack
    RANKS 2 THREADS 2)

hex_ecs_add_test(24-inner-ILU-mpi-distributed-mumps-out-of-core
    ARGS --groupsize 2 --shared-scratch --preconditioner ILU --out-of-core --lu mumps
    REQUIRES mumps
    RANKS 4 THREADS 1)

## --------------------------------------------------------------------------------- ##

# coupled preconditioner

hex_ecs_add_test(25-inner-coupled
    ARGS --groupsize 4 --shared-scratch --preconditioner coupled
         --lightweight-full --lu mumps
    REQUIRES mumps
    RANKS 4 THREADS 1)

hex_ecs_add_test(26-inner-coupled-channels
    ARGS --groupsize 4 --shared-scratch --preconditioner coupled
         --lightweight-full --couple-channels --lu mumps --mumps-out-of-core
    REQUIRES mumps
    RANKS 4 THREADS 1)

hex_ecs_add_test(27-channels-coupled
    ARGS --groupsize 4 --shared-scratch --preconditioner coupled
         --lightweight-full --lu mumps
    CHANNELS
    REQUIRES mumps
    RANKS 4 THREADS 1)

hex_ecs_add_test(28-channels-coupled-channels
    ARGS --groupsize 4 --shared-scratch --preconditioner coupled
         --couple-channels --lightweight-full --lu mumps
    CHANNELS
    REQUIRES mumps
    RANKS 4 THREADS 1)

## --------------------------------------------------------------------------------- ##

# hex_ecs_add_fold_test(<name>
#     [L <number>]                   total angular momentum (default 1)
#     [PI <number>]                  total parity (default 0)
#     [NL <number>]                  number of angular levels (default 1)
#     [ARGS <hex-ecs options>...]    options given to both runs of this test
#     [REQUIRES <feature>...]        one or more of umfpack, mumps, pardiso, opencl
#     [CHANNELS]                     extend the input file by projectile-only knots
#     [RANKS <n>]                    run under the MPI launcher with <n> processes
#     [THREADS <n>]                  value of OMP_NUM_THREADS for the test
# )
#
# Solves a small problem twice, with and without --fold-exchange, and checks that the
# partial cross sections agree. Unlike the tests above these are not physical calculations
# and run in seconds; see tests/cmake/run_fold_test.cmake for what is being compared.
function(hex_ecs_add_fold_test name)

    cmake_parse_arguments(TEST "CHANNELS" "L;PI;NL;RANKS;THREADS" "ARGS;REQUIRES" ${ARGN})

    if(TEST_UNPARSED_ARGUMENTS)
        message(FATAL_ERROR "hex_ecs_add_fold_test(${name}): unknown arguments ${TEST_UNPARSED_ARGUMENTS}")
    endif()

    # the default is the smallest basis in which some blocks can be folded at all
    # (dereference before comparing, so that an explicit zero is not overridden)
    if("${TEST_L}" STREQUAL "")
        set(TEST_L 1)
    endif()
    if("${TEST_PI}" STREQUAL "")
        set(TEST_PI 0)
    endif()
    if("${TEST_NL}" STREQUAL "")
        set(TEST_NL 1)
    endif()

    # an MPI test is only meaningful when there is something to launch it with
    if(TEST_RANKS)
        list(APPEND TEST_REQUIRES mpi)
    endif()

    # skip silently what this build cannot do; the missing features are summarized
    # once at the end of the configuration
    set(missing "")

    foreach(feature IN LISTS TEST_REQUIRES)
        if("${feature}" STREQUAL "umfpack" AND NOT WITH_UMFPACK)
            list(APPEND missing umfpack)
        elseif("${feature}" STREQUAL "mumps" AND NOT WITH_MUMPS)
            list(APPEND missing mumps)
        elseif("${feature}" STREQUAL "pardiso" AND NOT WITH_PARDISO AND NOT WITH_MKL)
            list(APPEND missing pardiso)
        elseif("${feature}" STREQUAL "opencl" AND NOT WITH_OPENCL)
            list(APPEND missing opencl)
        elseif("${feature}" STREQUAL "mpi" AND NOT MPIEXEC_EXECUTABLE)
            list(APPEND missing mpi)
        endif()
    endforeach()

    if(missing)
        set(hex-ecs_TESTS_SKIPPED ${hex-ecs_TESTS_SKIPPED} ${missing} PARENT_SCOPE)
        return()
    endif()

    ## ----------------------------------------------------------------------------- ##

    set(launcher "")

    if(TEST_RANKS)
        set(launcher
            "${MPIEXEC_EXECUTABLE}"
            ${MPIEXEC_NUMPROC_FLAG} ${TEST_RANKS}
            ${MPIEXEC_PREFLAGS}
        )
        list(APPEND TEST_ARGS --mpi)
    endif()

    add_test(
        NAME "hex-ecs.${name}"
        COMMAND ${CMAKE_COMMAND}
                -D "HEX_ECS=$<TARGET_FILE:hex-ecs>"
                -D "WORKDIR=${CMAKE_CURRENT_BINARY_DIR}/${name}"
                -D "ARGS=${TEST_ARGS}"
                -D "LAUNCHER=${launcher}"
                -D "CHANNELS=${TEST_CHANNELS}"
                -D "L=${TEST_L}"
                -D "PI=${TEST_PI}"
                -D "NL=${TEST_NL}"
                -P "${CMAKE_CURRENT_LIST_DIR}/cmake/run_fold_test.cmake"
    )

    set_tests_properties("hex-ecs.${name}" PROPERTIES
        TIMEOUT     ${HEX_ECS_TEST_TIMEOUT}
        RUN_SERIAL  ON
    )

    if(TEST_THREADS)
        set_tests_properties("hex-ecs.${name}" PROPERTIES
            ENVIRONMENT "OMP_NUM_THREADS=${TEST_THREADS}"
        )
    endif()

endfunction()

## --------------------------------------------------------------------------------- ##

# folding of the exchange symmetry into the angular basis
#
# The sample problem these use is small enough for a dense factorization, so most of them
# ask for "--lu lapack" and run in any build. That keeps the coverage of the option itself
# independent of which sparse factorizer happens to be available.

hex_ecs_add_fold_test(29-fold-inner-ILU
    ARGS --preconditioner ILU --lu lapack)

hex_ecs_add_fold_test(30-fold-inner-KPA
    ARGS --preconditioner KPA)

# The columns of an angular block row are dealt out to the members of an MPI group. There
# are more of them than there are solved blocks once the basis is folded, so the two do not
# line up as they do otherwise.
hex_ecs_add_fold_test(31-fold-inner-ILU-mpi-distributed
    ARGS --groupsize 2 --shared-scratch --preconditioner ILU --lu lapack
    RANKS 2 THREADS 1)

# three angular levels, so that blocks are folded across several of them
hex_ecs_add_fold_test(32-fold-inner-ILU-deep
    L 3 NL 2
    ARGS --preconditioner ILU --lu lapack)

hex_ecs_add_fold_test(33-fold-inner-ILU-lightweight
    ARGS --preconditioner ILU --lu lapack --lightweight-full)

hex_ecs_add_fold_test(34-fold-inner-ILU-out-of-core
    ARGS --preconditioner ILU --lu lapack --out-of-core)

# mirroring a solution segment has to swap the two groups of asymptotic channels, which
# only the tests with the projectile knots exercise
hex_ecs_add_fold_test(35-fold-channels-ILU
    CHANNELS
    ARGS --preconditioner ILU --lu lapack)

hex_ecs_add_fold_test(36-fold-channels-ILU-lightweight
    CHANNELS
    ARGS --preconditioner ILU --lu lapack --lightweight-full)

# the same with the factorizer that a production run would use
hex_ecs_add_fold_test(37-fold-channels-umfpack
    CHANNELS
    ARGS --preconditioner ILU
    REQUIRES umfpack)

hex_ecs_add_fold_test(38-fold-inner-ILU-mpi
    ARGS --shared-scratch --preconditioner ILU --lu lapack
    RANKS 2 THREADS 2)

## --------------------------------------------------------------------------------- ##

if(hex-ecs_TESTS_SKIPPED)
    list(REMOVE_DUPLICATES hex-ecs_TESTS_SKIPPED)
    list(SORT hex-ecs_TESTS_SKIPPED)
    string(REPLACE ";" ", " skipped "${hex-ecs_TESTS_SKIPPED}")
    message(STATUS "Some hex-ecs tests were left out, they need: ${skipped}")
endif()
