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

# Finds GiNaC, the symbolic algebra library, and the CLN it is built on.
#
# Looked for in this order:
#   1. GINAC_LIBRARIES, if the user has set it -- taken as given, not verified
#   2. ginac-config.cmake, exported by GiNaC 1.8 and later
#   3. ginac.pc, which GiNaC itself installs
#   4. a plain search for ginac/ginac.h and libginac
#
# GiNaC exposes CLN in its own headers, so CLN is part of the result in every case.
#
# Results: GINAC_FOUND, GINAC_INCLUDE_DIRS, GINAC_LIBRARIES, GiNaC::GiNaC

include("${CMAKE_CURRENT_LIST_DIR}/HexFindHelper.cmake")

if(GINAC_LIBRARIES)

    hex_find_result(GINAC
        REQUIRED_VARS GINAC_LIBRARIES
        TARGET        GiNaC::GiNaC
        INCLUDE_DIRS  ${GINAC_INCLUDE_DIRS}
        LIBRARIES     ${GINAC_LIBRARIES}
    )

    return()

endif()

find_package(ginac CONFIG QUIET)

if(ginac_FOUND AND TARGET ginac::ginac)

    set(ginac_FOUND ginac::ginac)

    hex_find_result(GINAC
        REQUIRED_VARS ginac_FOUND
        TARGET        GiNaC::GiNaC
        LIBRARIES     ginac::ginac
        VERSION       "${ginac_VERSION}"
    )

    return()

endif()

find_package(PkgConfig QUIET)

if(PKG_CONFIG_FOUND)
    pkg_check_modules(PC_GINAC QUIET ginac)
endif()

find_path(GINAC_INCLUDE_DIR
    NAMES ginac/ginac.h
    HINTS ${PC_GINAC_INCLUDE_DIRS}
)

find_library(GINAC_LIBRARY
    NAMES ginac
    HINTS ${PC_GINAC_LIBRARY_DIRS}
)

# the headers of GiNaC include those of CLN, and its symbols resolve against it
find_package(CLN QUIET)

set(_ginac_includes ${GINAC_INCLUDE_DIR})
set(_ginac_libraries ${GINAC_LIBRARY})

if(CLN_FOUND)
    list(APPEND _ginac_includes ${CLN_INCLUDE_DIRS})
    list(APPEND _ginac_libraries ${CLN_LIBRARIES})
endif()

if(_ginac_includes)
    list(REMOVE_DUPLICATES _ginac_includes)
endif()

hex_find_result(GINAC
    REQUIRED_VARS GINAC_LIBRARY GINAC_INCLUDE_DIR
    TARGET        GiNaC::GiNaC
    INCLUDE_DIRS  ${_ginac_includes}
    LIBRARIES     ${_ginac_libraries}
    VERSION       "${PC_GINAC_VERSION}"
)

unset(_ginac_includes)
unset(_ginac_libraries)

mark_as_advanced(GINAC_INCLUDE_DIR GINAC_LIBRARY)
