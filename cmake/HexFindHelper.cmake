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

# Shared tail of the Find modules of Hex.
#
# hex_find_result(<PKG>
#     REQUIRED_VARS <var>...        variables that must be set for success
#     [TARGET       <target>]       imported target to create
#     [INCLUDE_DIRS <dir>...]       what to publish as <PKG>_INCLUDE_DIRS
#     [LIBRARIES    <lib>...]       what to publish as <PKG>_LIBRARIES
#     [VERSION      <version>]      version to report and to check a request against
# )
#
# Reports the outcome through find_package_handle_standard_args, which honours the
# REQUIRED and QUIET keywords of the find_package call and prints the usual
# "-- Found <PKG>: ..." line, then publishes the result variables and wraps them in
# an interface imported target.
#
# This is a macro on purpose: find_package_handle_standard_args reads and writes
# variables of the scope it is called from, and a macro does not add one.

include(FindPackageHandleStandardArgs)

macro(hex_find_result PKG)

    cmake_parse_arguments(HFR
        ""
        "TARGET;VERSION"
        "REQUIRED_VARS;INCLUDE_DIRS;LIBRARIES"
        ${ARGN}
    )

    if(HFR_VERSION)
        find_package_handle_standard_args(${PKG}
            REQUIRED_VARS ${HFR_REQUIRED_VARS}
            VERSION_VAR   HFR_VERSION
        )
    else()
        find_package_handle_standard_args(${PKG}
            REQUIRED_VARS ${HFR_REQUIRED_VARS}
        )
    endif()

    if(${PKG}_FOUND)

        set(${PKG}_INCLUDE_DIRS ${HFR_INCLUDE_DIRS})
        set(${PKG}_LIBRARIES    ${HFR_LIBRARIES})

        if(HFR_TARGET AND NOT TARGET ${HFR_TARGET})

            add_library(${HFR_TARGET} INTERFACE IMPORTED GLOBAL)

            if(${PKG}_INCLUDE_DIRS)
                set_target_properties(${HFR_TARGET} PROPERTIES
                    INTERFACE_INCLUDE_DIRECTORIES "${${PKG}_INCLUDE_DIRS}"
                )
            endif()

            if(${PKG}_LIBRARIES)
                set_target_properties(${HFR_TARGET} PROPERTIES
                    INTERFACE_LINK_LIBRARIES "${${PKG}_LIBRARIES}"
                )
            endif()

        endif()

    endif()

    unset(HFR_TARGET)
    unset(HFR_VERSION)
    unset(HFR_REQUIRED_VARS)
    unset(HFR_INCLUDE_DIRS)
    unset(HFR_LIBRARIES)
    unset(HFR_UNPARSED_ARGUMENTS)

endmacro()

# Directories where the MPI library itself lives. Distributions that keep several
# MPI flavours side by side (openSUSE, Fedora, Debian) install the MPI-dependent
# solvers -- ScaLAPACK, MUMPS, SuperLU_DIST -- next to it rather than in /usr/lib,
# where find_library would look.
function(hex_mpi_library_dirs out)

    set(dirs "")

    foreach(library IN LISTS MPI_C_LIBRARIES MPI_CXX_LIBRARIES)
        if(EXISTS "${library}")
            get_filename_component(directory "${library}" DIRECTORY)
            list(APPEND dirs "${directory}")
        endif()
    endforeach()

    if(dirs)
        list(REMOVE_DUPLICATES dirs)
    endif()

    set(${out} "${dirs}" PARENT_SCOPE)

endfunction()
