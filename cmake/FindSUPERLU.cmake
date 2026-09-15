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

# Finds SuperLU, the serial sparse LU factorization.
#
# Looked for in this order:
#   1. SUPERLU_LIBRARIES, if the user has set it -- taken as given, not verified
#   2. superluConfig.cmake, exported by SuperLU 5.3 and later
#   3. superlu.pc, shipped by most distributions
#   4. a plain search for slu_zdefs.h and libsuperlu
#
# Results: SUPERLU_FOUND, SUPERLU_INCLUDE_DIRS, SUPERLU_LIBRARIES, SuperLU::SuperLU

include("${CMAKE_CURRENT_LIST_DIR}/HexFindHelper.cmake")

if(SUPERLU_LIBRARIES)

    hex_find_result(SUPERLU
        REQUIRED_VARS SUPERLU_LIBRARIES
        TARGET        SuperLU::SuperLU
        INCLUDE_DIRS  ${SUPERLU_INCLUDE_DIRS}
        LIBRARIES     ${SUPERLU_LIBRARIES}
    )

    return()

endif()

find_package(superlu CONFIG QUIET)

if(superlu_FOUND AND TARGET superlu::superlu)

    set(superlu_FOUND superlu::superlu)

    hex_find_result(SUPERLU
        REQUIRED_VARS superlu_FOUND
        TARGET        SuperLU::SuperLU
        LIBRARIES     superlu::superlu
        VERSION       "${superlu_VERSION}"
    )

    return()

endif()

find_package(PkgConfig QUIET)

if(PKG_CONFIG_FOUND)
    pkg_check_modules(PC_SUPERLU QUIET superlu)
endif()

# hex-ecs includes slu_zdefs.h in FP64 and slu_cdefs.h in FP32; both sit next to
# each other, so looking for one of them is enough
find_path(SUPERLU_INCLUDE_DIR
    NAMES slu_zdefs.h
    HINTS ${PC_SUPERLU_INCLUDE_DIRS}
    PATH_SUFFIXES superlu SuperLU
)

find_library(SUPERLU_LIBRARY
    NAMES superlu
    HINTS ${PC_SUPERLU_LIBRARY_DIRS}
)

hex_find_result(SUPERLU
    REQUIRED_VARS SUPERLU_LIBRARY SUPERLU_INCLUDE_DIR
    TARGET        SuperLU::SuperLU
    INCLUDE_DIRS  ${SUPERLU_INCLUDE_DIR}
    LIBRARIES     ${SUPERLU_LIBRARY}
    VERSION       "${PC_SUPERLU_VERSION}"
)

mark_as_advanced(SUPERLU_INCLUDE_DIR SUPERLU_LIBRARY)
