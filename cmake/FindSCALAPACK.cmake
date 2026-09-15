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

# Finds ScaLAPACK, the distributed dense linear algebra library.
#
# Looked for in this order:
#   1. SCALAPACK_LIBRARIES, if the user has set it -- taken as given
#   2. MKL::MKL_SCALAPACK, when the build already uses Intel MKL
#   3. scalapack.pc, shipped by some distributions
#   4. a plain search for libscalapack, including the directories of the MPI library
#
# hex-ecs declares the ScaLAPACK routines itself, so no header is needed -- the
# library alone is the whole dependency.
#
# Results: SCALAPACK_FOUND, SCALAPACK_INCLUDE_DIRS, SCALAPACK_LIBRARIES,
#          ScaLAPACK::ScaLAPACK

include("${CMAKE_CURRENT_LIST_DIR}/HexFindHelper.cmake")

if(SCALAPACK_LIBRARIES)

    hex_find_result(SCALAPACK
        REQUIRED_VARS SCALAPACK_LIBRARIES
        TARGET        ScaLAPACK::ScaLAPACK
        INCLUDE_DIRS  ${SCALAPACK_INCLUDE_DIRS}
        LIBRARIES     ${SCALAPACK_LIBRARIES}
    )

    return()

endif()

# MKL brings its own ScaLAPACK, and mixing it with the Netlib one does not work
if(TARGET MKL::MKL_SCALAPACK)

    set(SCALAPACK_MKL_TARGET MKL::MKL_SCALAPACK)

    hex_find_result(SCALAPACK
        REQUIRED_VARS SCALAPACK_MKL_TARGET
        TARGET        ScaLAPACK::ScaLAPACK
        LIBRARIES     MKL::MKL_SCALAPACK
    )

    return()

endif()

find_package(PkgConfig QUIET)

if(PKG_CONFIG_FOUND)
    pkg_check_modules(PC_SCALAPACK QUIET scalapack)
endif()

hex_mpi_library_dirs(_mpi_dirs)

find_library(SCALAPACK_LIBRARY
    NAMES scalapack scalapack-openmpi scalapack-mpich scalapack-mpich2 scalapack-lam
    HINTS ${PC_SCALAPACK_LIBRARY_DIRS} ${_mpi_dirs}
)

# builds that keep BLACS apart from ScaLAPACK
find_library(SCALAPACK_BLACS_LIBRARY
    NAMES blacs-openmpi blacs-mpich blacs
    HINTS ${PC_SCALAPACK_LIBRARY_DIRS} ${_mpi_dirs}
)

set(_scalapack_libraries ${SCALAPACK_LIBRARY})

if(SCALAPACK_BLACS_LIBRARY)
    list(APPEND _scalapack_libraries "${SCALAPACK_BLACS_LIBRARY}")
endif()

hex_find_result(SCALAPACK
    REQUIRED_VARS SCALAPACK_LIBRARY
    TARGET        ScaLAPACK::ScaLAPACK
    LIBRARIES     ${_scalapack_libraries}
    VERSION       "${PC_SCALAPACK_VERSION}"
)

unset(_mpi_dirs)
unset(_scalapack_libraries)

mark_as_advanced(SCALAPACK_LIBRARY SCALAPACK_BLACS_LIBRARY)
