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

# Finds Intel MKL, used by hex-ecs for its PARDISO sparse LU factorization.
#
# Looked for in this order:
#   1. MKL_LIBRARIES, if the user has set it -- taken as given, not verified
#   2. MKLConfig.cmake, which every oneAPI installation ships
#
# The config package is what builds the MKL link line, whose order and choice of
# threading and interface layer are easy to get wrong by hand. It is steered by the
# cache variables MKL_INTERFACE, MKL_THREADING, MKL_LINK and MKL_MPI, and it also
# exports MKL::MKL_SCALAPACK, which FindSCALAPACK picks up.
#
# The interface layer defaults to lp64 here, not to the ilp64 default of MKL
# itself: the PARDISO prototypes in hex-ecs/src/factorizers/lu-pardiso.h declare
# their integer arguments as "int". This is independent of _LONGINT, which only
# selects the integer type Hex uses towards UMFPACK and SuperLU.
#
# Results: MKL_FOUND, MKL_INCLUDE_DIRS, MKL_LIBRARIES, MKL::MKL

include("${CMAKE_CURRENT_LIST_DIR}/HexFindHelper.cmake")

if(MKL_LIBRARIES)

    hex_find_result(MKL
        REQUIRED_VARS MKL_LIBRARIES
        TARGET        MKL::MKL
        INCLUDE_DIRS  ${MKL_INCLUDE_DIRS}
        LIBRARIES     ${MKL_LIBRARIES}
    )

    return()

endif()

if(NOT DEFINED MKL_INTERFACE)
    set(MKL_INTERFACE "lp64" CACHE STRING "MKL integer interface layer (lp64 or ilp64)")
endif()

# setvars.sh of oneAPI exports MKLROOT, while find_package looks for MKL_ROOT
if(NOT MKL_ROOT AND NOT DEFINED ENV{MKL_ROOT} AND DEFINED ENV{MKLROOT})
    set(MKL_ROOT "$ENV{MKLROOT}")
endif()

find_package(MKL CONFIG QUIET)

if(MKL_FOUND AND TARGET MKL::MKL)
    message(STATUS "Found MKL: ${MKL_ROOT} (interface ${MKL_INTERFACE})")
    return()
endif()

hex_find_result(MKL
    REQUIRED_VARS MKL_ROOT
    TARGET        MKL::MKL
)
