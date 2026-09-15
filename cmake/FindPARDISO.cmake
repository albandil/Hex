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

# Finds the stand-alone PARDISO library from pardiso-project.org.
#
# This is the separately licensed library, not the PARDISO that comes inside Intel
# MKL; for the latter use WITH_MKL. hex-ecs declares its entry points itself, so
# there is no header to find, and the library is distributed as a single shared
# object that the user usually has to point at:
#
#   -D PARDISO_LIBRARIES=/opt/pardiso/libpardiso700-GNU831-X86-64.so
#
# Results: PARDISO_FOUND, PARDISO_INCLUDE_DIRS, PARDISO_LIBRARIES, PARDISO::PARDISO

include("${CMAKE_CURRENT_LIST_DIR}/HexFindHelper.cmake")

if(PARDISO_LIBRARIES)

    hex_find_result(PARDISO
        REQUIRED_VARS PARDISO_LIBRARIES
        TARGET        PARDISO::PARDISO
        INCLUDE_DIRS  ${PARDISO_INCLUDE_DIRS}
        LIBRARIES     ${PARDISO_LIBRARIES}
    )

    return()

endif()

# the release archives carry the compiler and the architecture in the file name,
# so a bare "pardiso" is only one of the possibilities
find_library(PARDISO_LIBRARY
    NAMES pardiso pardiso700 pardiso600 pardiso500
)

hex_find_result(PARDISO
    REQUIRED_VARS PARDISO_LIBRARY
    TARGET        PARDISO::PARDISO
    LIBRARIES     ${PARDISO_LIBRARY}
)

mark_as_advanced(PARDISO_LIBRARY)
