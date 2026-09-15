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

include(FeatureSummary)

# hex_dependency(<name>
#     ENABLED       <bool>           whether the dependency is available and wanted
#     [DEFINE       <macro>]         macro to define for the code that uses it
#     [DESCRIPTION  <text>]          one line for the configuration summary
#     [TARGETS      <target>...]     imported targets to link (from find_package)
#     [INCLUDE_DIRS <dir>...]        include directories to add (hand-specified)
#     [LIBRARIES    <lib>...]        libraries to link (hand-specified)
# )
#
# Declares the interface library Hex::<name>. The target is created whether or not
# the dependency is enabled, so that the subdirectories can link it unconditionally
# and stay free of if(WITH_...) blocks: a disabled dependency is simply an empty
# target that contributes nothing.
#
# Dependencies found by find_package are passed as TARGETS, which carry their own
# include directories, libraries and flags. The ones the user still has to point at
# by hand are passed as INCLUDE_DIRS and LIBRARIES; those show up in the summary,
# as there is no "-- Found ..." line to tell what was actually used.
function(hex_dependency name)

    # the plain ${ARGN} form, not PARSE_ARGV: the latter escapes the semicolons
    # inside each argument, which would turn a list of libraries into one element
    cmake_parse_arguments(HD
        ""
        "ENABLED;DEFINE;DESCRIPTION"
        "TARGETS;INCLUDE_DIRS;LIBRARIES"
        ${ARGN}
    )

    if(DEFINED HD_UNPARSED_ARGUMENTS)
        message(FATAL_ERROR "hex_dependency(${name}): unexpected argument(s) ${HD_UNPARSED_ARGUMENTS}")
    endif()

    add_library(hex-dep-${name} INTERFACE)
    add_library(Hex::${name} ALIAS hex-dep-${name})

    if(HD_ENABLED)

        if(HD_DEFINE)
            target_compile_definitions(hex-dep-${name} INTERFACE "${HD_DEFINE}")
        endif()

        # SYSTEM: these are somebody else's headers, do not warn about them
        if(HD_INCLUDE_DIRS)
            target_include_directories(hex-dep-${name} SYSTEM INTERFACE ${HD_INCLUDE_DIRS})
        endif()

        foreach(target IN LISTS HD_TARGETS)
            if(NOT TARGET "${target}")
                message(FATAL_ERROR "hex_dependency(${name}): no such target ${target}")
            endif()
        endforeach()

        if(HD_TARGETS OR HD_LIBRARIES)
            target_link_libraries(hex-dep-${name} INTERFACE ${HD_TARGETS} ${HD_LIBRARIES})
        endif()

    endif()

    # record the outcome for the summary printed at the end of the configuration
    set(summary "${HD_DESCRIPTION}")
    if(HD_ENABLED AND HD_LIBRARIES)
        string(REPLACE ";" " " libraries "${HD_LIBRARIES}")
        set(summary "${summary} -- using ${libraries}")
    endif()

    # The feature is registered under its target name, not under "${name}": a bare
    # name that happens to match a package someone called find_package on -- superlu
    # and ginac both do -- is taken for that package by FeatureSummary and dropped
    # from the feature lists. The target name is also what the subdirectories link.
    if(HD_ENABLED)
        add_feature_info("Hex::${name}" TRUE "${summary}")
    else()
        add_feature_info("Hex::${name}" FALSE "${summary}")
    endif()

endfunction()
