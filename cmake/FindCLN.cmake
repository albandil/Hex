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

# Finds CLN, the Class Library for Numbers.
#
# Looked for in this order:
#   1. CLN_LIBRARIES, if the user has set it -- taken as given, not verified
#   2. cln.pc, which CLN itself installs
#   3. a plain search for cln/cln.h and libcln
#
# Results: CLN_FOUND, CLN_INCLUDE_DIRS, CLN_LIBRARIES, CLN::CLN and cln::cln

include("${CMAKE_CURRENT_LIST_DIR}/HexFindHelper.cmake")

if(CLN_LIBRARIES)

    # taken as given
    set(_cln_includes ${CLN_INCLUDE_DIRS})
    set(_cln_libraries ${CLN_LIBRARIES})
    set(_cln_required CLN_LIBRARIES)
    set(_cln_version "")

else()

    find_package(PkgConfig QUIET)

    if(PKG_CONFIG_FOUND)
        pkg_check_modules(PC_CLN QUIET cln)
    endif()

    find_path(CLN_INCLUDE_DIR
        NAMES cln/cln.h
        HINTS ${PC_CLN_INCLUDE_DIRS}
    )

    find_library(CLN_LIBRARY
        NAMES cln
        HINTS ${PC_CLN_LIBRARY_DIRS}
    )

    set(_cln_includes ${CLN_INCLUDE_DIR})
    set(_cln_libraries ${CLN_LIBRARY})
    set(_cln_required CLN_LIBRARY CLN_INCLUDE_DIR)
    set(_cln_version "${PC_CLN_VERSION}")

endif()

hex_find_result(CLN
    REQUIRED_VARS ${_cln_required}
    TARGET        CLN::CLN
    INCLUDE_DIRS  ${_cln_includes}
    LIBRARIES     ${_cln_libraries}
    VERSION       "${_cln_version}"
)

# GiNaC's installed config package does find_package(CLN REQUIRED) and then links
# the target cln::cln, the name created by the FindCLN.cmake that GiNaC bundles but
# does not install. Provide that name as well, so that ginac-config.cmake resolves
# against this module instead of failing on an unknown target.
if(CLN_FOUND AND NOT TARGET cln::cln)

    add_library(cln::cln INTERFACE IMPORTED GLOBAL)

    if(CLN_INCLUDE_DIRS)
        set_target_properties(cln::cln PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${CLN_INCLUDE_DIRS}"
        )
    endif()

    if(CLN_LIBRARIES)
        set_target_properties(cln::cln PROPERTIES
            INTERFACE_LINK_LIBRARIES "${CLN_LIBRARIES}"
        )
    endif()

endif()

unset(_cln_includes)
unset(_cln_libraries)
unset(_cln_required)
unset(_cln_version)

mark_as_advanced(CLN_INCLUDE_DIR CLN_LIBRARY)
