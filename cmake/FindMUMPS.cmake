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

# Finds MUMPS, the distributed multifrontal sparse solver.
#
# Looked for in this order:
#   1. MUMPS_LIBRARIES, if the user has set it -- taken as given, not verified
#   2. MUMPSConfig.cmake, exported by recent releases and by the spack/conda builds
#   3. a plain search for the headers and the libraries
#
# MUMPS is the dependency that most often has to be pointed at by hand: the library
# is built in several flavours (sequential, with ParMETIS, with PT-Scotch) which
# distributions install side by side with decorated names, and the ordering
# libraries it needs are not recorded anywhere. Setting MUMPS_LIBRARIES to a
# complete link line, "-L..." entries included, bypasses the search entirely, e.g.
#
#   -D MUMPS_INCLUDE_DIRS=/usr/include/mumps
#   -D MUMPS_LIBRARIES="-L/usr/lib64/mpi/gcc/openmpi5/lib64;zmumps_ptscotch;mumps_common_ptscotch;ptesmumps;ptscotch;ptscotcherr"
#
# Results: MUMPS_FOUND, MUMPS_INCLUDE_DIRS, MUMPS_LIBRARIES, MUMPS::MUMPS

include("${CMAKE_CURRENT_LIST_DIR}/HexFindHelper.cmake")

if(MUMPS_LIBRARIES)

    hex_find_result(MUMPS
        REQUIRED_VARS MUMPS_LIBRARIES
        TARGET        MUMPS::MUMPS
        INCLUDE_DIRS  ${MUMPS_INCLUDE_DIRS}
        LIBRARIES     ${MUMPS_LIBRARIES}
    )

    return()

endif()

find_package(MUMPS CONFIG QUIET)

if(MUMPS_FOUND AND TARGET MUMPS::MUMPS)
    return()
endif()

hex_mpi_library_dirs(_mpi_dirs)

# hex-ecs uses the complex-valued interface: zmumps_c.h in FP64, cmumps_c.h in FP32
find_path(MUMPS_INCLUDE_DIR
    NAMES zmumps_c.h
    PATH_SUFFIXES mumps MUMPS mumps_seq openmpi-x86_64 mpich-x86_64
)

# the flavour suffixes of the distributions, plainest name first
set(_mumps_suffixes "" _ptscotch _parmetis _scotch _metis _seq)

set(_mumps_z_names "")
set(_mumps_common_names "")

foreach(_suffix IN LISTS _mumps_suffixes)
    list(APPEND _mumps_z_names      "zmumps${_suffix}")
    list(APPEND _mumps_common_names "mumps_common${_suffix}")
endforeach()

find_library(MUMPS_zmumps_LIBRARY NAMES ${_mumps_z_names}      HINTS ${_mpi_dirs})
find_library(MUMPS_common_LIBRARY NAMES ${_mumps_common_names} HINTS ${_mpi_dirs})

# the built-in ordering of MUMPS; present in most builds, not required by all
find_library(MUMPS_pord_LIBRARY NAMES pord pord_seq HINTS ${_mpi_dirs})

set(_mumps_libraries ${MUMPS_zmumps_LIBRARY} ${MUMPS_common_LIBRARY})

if(MUMPS_pord_LIBRARY)
    list(APPEND _mumps_libraries "${MUMPS_pord_LIBRARY}")
endif()

hex_find_result(MUMPS
    REQUIRED_VARS MUMPS_zmumps_LIBRARY MUMPS_common_LIBRARY MUMPS_INCLUDE_DIR
    TARGET        MUMPS::MUMPS
    INCLUDE_DIRS  ${MUMPS_INCLUDE_DIR}
    LIBRARIES     ${_mumps_libraries}
)

unset(_mpi_dirs)
unset(_mumps_suffixes)
unset(_mumps_z_names)
unset(_mumps_common_names)
unset(_mumps_libraries)

mark_as_advanced(MUMPS_INCLUDE_DIR MUMPS_zmumps_LIBRARY MUMPS_common_LIBRARY MUMPS_pord_LIBRARY)
