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

# Tests of the C and Fortran bindings of libhex-db. Both programs fill a small
# database with a few T-matrices and print the complete cross section computed
# from them; the expected numbers are checked by CTest.
#
# The Fortran binding is only tested when a Fortran compiler is available. It is
# the sole use of Fortran in the project, so the language is enabled here and not
# in the top-level CMakeLists.txt.
#
# The third test is a regression test of the ionization quantities; it reports its
# own verdict through the exit code.

# the three cross sections that the sample T-matrices give, one per line
set(HEX_DB_EXPECTED "16\\.5930[^\n]*\n[^\n]*14\\.0468[^\n]*\n[^\n]*12\\.3309")

## --------------------------------------------------------------------------------- ##

add_executable(test-hexdb-lib-c "${CMAKE_CURRENT_LIST_DIR}/test-hexdb-lib-c.c")

# keep the test programs out of the directory with the installed ones
set_target_properties(test-hexdb-lib-c PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

target_link_libraries(test-hexdb-lib-c PRIVATE libhex-db)

# the programs write their database and SQL batch into the current directory
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/c")

add_test(
    NAME "hex-db.binding-c"
    COMMAND test-hexdb-lib-c
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/c"
)

set_tests_properties("hex-db.binding-c" PROPERTIES
    PASS_REGULAR_EXPRESSION "${HEX_DB_EXPECTED}"
)

## --------------------------------------------------------------------------------- ##

include(CheckLanguage)
check_language(Fortran)

if(CMAKE_Fortran_COMPILER)

    enable_language(Fortran)

    add_executable(test-hexdb-lib-f "${CMAKE_CURRENT_LIST_DIR}/test-hexdb-lib-f.f90")

    set_target_properties(test-hexdb-lib-f PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    )

    target_link_libraries(test-hexdb-lib-f PRIVATE libhex-db)

    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran")

    add_test(
        NAME "hex-db.binding-fortran"
        COMMAND test-hexdb-lib-f
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran"
    )

    set_tests_properties("hex-db.binding-fortran" PROPERTIES
        PASS_REGULAR_EXPRESSION "${HEX_DB_EXPECTED}"
    )

else()

    message(STATUS "The hex-db Fortran binding test was left out, no Fortran compiler")

endif()

## --------------------------------------------------------------------------------- ##

add_executable(test-hexdb-ionization "${CMAKE_CURRENT_LIST_DIR}/test-hexdb-ionization.cpp")

set_target_properties(test-hexdb-ionization PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

target_link_libraries(test-hexdb-ionization PRIVATE libhex-common libhex-db)

target_include_directories(test-hexdb-ionization PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/../../common"
    "${CMAKE_CURRENT_LIST_DIR}/../../libs"
)

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/ionization")

add_test(
    NAME "hex-db.ionization"
    COMMAND test-hexdb-ionization
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/ionization"
)

## --------------------------------------------------------------------------------- ##
