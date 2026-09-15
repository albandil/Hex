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

# Finds UMFPACK, the sparse LU factorization of SuiteSparse.
#
# Looked for in this order:
#   1. UMFPACK_LIBRARIES, if the user has set it -- taken as given, not verified
#   2. UMFPACKConfig.cmake, exported by SuiteSparse 6 and later
#   3. umfpack.pc, shipped by some distributions
#   4. a plain search for umfpack.h and libumfpack
#
# Results: UMFPACK_FOUND, UMFPACK_INCLUDE_DIRS, UMFPACK_LIBRARIES, UMFPACK::UMFPACK

include("${CMAKE_CURRENT_LIST_DIR}/HexFindHelper.cmake")

if(UMFPACK_LIBRARIES)

    hex_find_result(UMFPACK
        REQUIRED_VARS UMFPACK_LIBRARIES
        TARGET        UMFPACK::UMFPACK
        INCLUDE_DIRS  ${UMFPACK_INCLUDE_DIRS}
        LIBRARIES     ${UMFPACK_LIBRARIES}
    )

    return()

endif()

# the config package of SuiteSparse; CONFIG mode, so this does not recurse into
# the present module
find_package(UMFPACK CONFIG QUIET)

if(UMFPACK_FOUND AND TARGET SuiteSparse::UMFPACK)

    hex_find_result(UMFPACK
        REQUIRED_VARS UMFPACK_LIBRARY
        TARGET        UMFPACK::UMFPACK
        LIBRARIES     SuiteSparse::UMFPACK
        VERSION       "${UMFPACK_VERSION}"
    )

    return()

endif()

find_package(PkgConfig QUIET)

if(PKG_CONFIG_FOUND)
    pkg_check_modules(PC_UMFPACK QUIET umfpack)
endif()

find_path(UMFPACK_INCLUDE_DIR
    NAMES umfpack.h
    HINTS ${PC_UMFPACK_INCLUDE_DIRS}
    PATH_SUFFIXES suitesparse SuiteSparse ufsparse
)

find_library(UMFPACK_LIBRARY
    NAMES umfpack
    HINTS ${PC_UMFPACK_LIBRARY_DIRS}
)

# UMFPACK needs AMD and the SuiteSparse configuration object; shared builds carry
# them as dependencies, static ones do not
find_library(UMFPACK_AMD_LIBRARY NAMES amd HINTS ${PC_UMFPACK_LIBRARY_DIRS})
find_library(UMFPACK_CONFIG_LIBRARY NAMES suitesparseconfig HINTS ${PC_UMFPACK_LIBRARY_DIRS})

set(_umfpack_libraries ${UMFPACK_LIBRARY})

foreach(_extra ${UMFPACK_AMD_LIBRARY} ${UMFPACK_CONFIG_LIBRARY})
    list(APPEND _umfpack_libraries "${_extra}")
endforeach()

hex_find_result(UMFPACK
    REQUIRED_VARS UMFPACK_LIBRARY UMFPACK_INCLUDE_DIR
    TARGET        UMFPACK::UMFPACK
    INCLUDE_DIRS  ${UMFPACK_INCLUDE_DIR}
    LIBRARIES     ${_umfpack_libraries}
    VERSION       "${PC_UMFPACK_VERSION}"
)

unset(_umfpack_libraries)

mark_as_advanced(UMFPACK_INCLUDE_DIR UMFPACK_LIBRARY UMFPACK_AMD_LIBRARY UMFPACK_CONFIG_LIBRARY)
