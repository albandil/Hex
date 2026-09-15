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

# Finds SuperLU_DIST, the distributed sparse LU factorization.
#
# Looked for in this order:
#   1. SUPERLU_DIST_LIBRARIES, if the user has set it -- taken as given
#   2. superlu_distConfig.cmake, exported by recent releases
#   3. a plain search for superlu_zdefs.h and libsuperlu_dist
#
# The search covers the directories of the MPI library, because distributions that
# offer several MPI flavours install SuperLU_DIST next to the flavour it was built
# against instead of into the default library path.
#
# Results: SUPERLU_DIST_FOUND, SUPERLU_DIST_INCLUDE_DIRS, SUPERLU_DIST_LIBRARIES,
#          SuperLU::SuperLUDist

include("${CMAKE_CURRENT_LIST_DIR}/HexFindHelper.cmake")

if(SUPERLU_DIST_LIBRARIES)

    hex_find_result(SUPERLU_DIST
        REQUIRED_VARS SUPERLU_DIST_LIBRARIES
        TARGET        SuperLU::SuperLUDist
        INCLUDE_DIRS  ${SUPERLU_DIST_INCLUDE_DIRS}
        LIBRARIES     ${SUPERLU_DIST_LIBRARIES}
    )

    return()

endif()

find_package(superlu_dist CONFIG QUIET)

if(superlu_dist_FOUND AND TARGET superlu_dist::superlu_dist)

    set(superlu_dist_FOUND superlu_dist::superlu_dist)

    hex_find_result(SUPERLU_DIST
        REQUIRED_VARS superlu_dist_FOUND
        TARGET        SuperLU::SuperLUDist
        LIBRARIES     superlu_dist::superlu_dist
        VERSION       "${superlu_dist_VERSION}"
    )

    return()

endif()

hex_mpi_library_dirs(_mpi_dirs)

find_path(SUPERLU_DIST_INCLUDE_DIR
    NAMES superlu_zdefs.h
    PATH_SUFFIXES superlu_dist superlu-dist SuperLU_DIST
)

find_library(SUPERLU_DIST_LIBRARY
    NAMES superlu_dist superlu-dist
    HINTS ${_mpi_dirs}
)

hex_find_result(SUPERLU_DIST
    REQUIRED_VARS SUPERLU_DIST_LIBRARY SUPERLU_DIST_INCLUDE_DIR
    TARGET        SuperLU::SuperLUDist
    INCLUDE_DIRS  ${SUPERLU_DIST_INCLUDE_DIR}
    LIBRARIES     ${SUPERLU_DIST_LIBRARY}
)

unset(_mpi_dirs)

mark_as_advanced(SUPERLU_DIST_INCLUDE_DIR SUPERLU_DIST_LIBRARY)
