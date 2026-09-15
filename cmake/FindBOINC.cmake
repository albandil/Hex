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

# Finds the BOINC client libraries, needed to build Hex as a BOINC application.
#
# Looked for in this order:
#   1. BOINC_LIBRARIES, if the user has set it -- taken as given, not verified
#   2. a plain search for boinc_api.h, libboinc_api and libboinc
#
# The two libraries have to be linked in this order: libboinc_api calls into
# libboinc, not the other way round.
#
# Results: BOINC_FOUND, BOINC_INCLUDE_DIRS, BOINC_LIBRARIES, BOINC::BOINC

include("${CMAKE_CURRENT_LIST_DIR}/HexFindHelper.cmake")

if(BOINC_LIBRARIES)

    hex_find_result(BOINC
        REQUIRED_VARS BOINC_LIBRARIES
        TARGET        BOINC::BOINC
        INCLUDE_DIRS  ${BOINC_INCLUDE_DIRS}
        LIBRARIES     ${BOINC_LIBRARIES}
    )

    return()

endif()

find_path(BOINC_INCLUDE_DIR
    NAMES boinc_api.h
    PATH_SUFFIXES boinc BOINC
)

find_library(BOINC_API_LIBRARY NAMES boinc_api)
find_library(BOINC_LIBRARY     NAMES boinc)

hex_find_result(BOINC
    REQUIRED_VARS BOINC_API_LIBRARY BOINC_LIBRARY BOINC_INCLUDE_DIR
    TARGET        BOINC::BOINC
    INCLUDE_DIRS  ${BOINC_INCLUDE_DIR}
    LIBRARIES     ${BOINC_API_LIBRARY} ${BOINC_LIBRARY}
)

mark_as_advanced(BOINC_INCLUDE_DIR BOINC_API_LIBRARY BOINC_LIBRARY)
