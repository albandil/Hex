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

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

set(libhex-common_SOURCES
    "hex-arrays.cpp"
    "hex-blas.cpp"
    "hex-born.cpp"
    "hex-hdffile.cpp"
    "hex-hydrogen.cpp"
    "hex-matrix.cpp"
    "hex-misc.cpp"
    "hex-special.cpp"
    "hex-spgrid.cpp"
    "hex-symbolic.cpp"
    "hex-vec3d.cpp"
    "hex-version.cpp"
    "hex-vtkfile.cpp"
)

if(HDF5_FOUND)
    set(libhex-common_SOURCES ${libhex-common_SOURCES} "hex-h5file.cpp")
endif()

add_library(libhex-common SHARED ${libhex-common_SOURCES})

set_target_properties(libhex-common PROPERTIES PREFIX "")

target_include_directories(libhex-common PUBLIC
    ../libs
)

## --------------------------------------------------------------------------------- ##

# The version and git commit shown in the "logo" banner (hex-version.cpp) are
# derived from git rather than hand-maintained: the version is the name of the
# nearest reachable tag (e.g. v2.03) and the commit is the short hash. Both fall
# back to "unknown" outside a git checkout, e.g. a source tarball.
#
# hex-version-git.h is (re)written by a custom target that runs on every build,
# not only on reconfiguration, so a rebuild always reflects the current commit;
# it is only rewritten when its content changes, so this does not force
# hex-version.cpp to recompile except when the version or commit actually moved.

find_package(Git QUIET)

set(HEX_VERSION_HEADER "${CMAKE_CURRENT_BINARY_DIR}/hex-version-git.h")

add_custom_target(hex-version-check
    COMMAND ${CMAKE_COMMAND}
        "-DHEX_SOURCE_DIR=${PROJECT_SOURCE_DIR}"
        "-DHEX_GIT_EXECUTABLE=${GIT_EXECUTABLE}"
        "-DHEX_OUTPUT_FILE=${HEX_VERSION_HEADER}"
        -P "${PROJECT_SOURCE_DIR}/cmake/GenerateHexVersion.cmake"
    BYPRODUCTS "${HEX_VERSION_HEADER}"
    COMMENT "Determining Hex version from git"
    VERBATIM
)

add_dependencies(libhex-common hex-version-check)

target_include_directories(libhex-common PRIVATE
    "${CMAKE_CURRENT_BINARY_DIR}"
)

# PUBLIC throughout: the installed headers of libhex-common include the headers of
# GSL, HDF5, CLN/GiNaC and png++, so its consumers need them too. The Hex::* targets
# are empty for the dependencies that are switched off, hence no conditionals here.
target_link_libraries(libhex-common PUBLIC
    Hex::config
    Hex::blas
    Hex::boinc
    Hex::cln
    Hex::ginac
    Hex::gsl
    Hex::hdf5
    Hex::lapack
    Hex::mpi
    Hex::openmp
    Hex::png
)

## --------------------------------------------------------------------------------- ##

if(BUILD_TESTING)
    add_subdirectory(test)
endif()

## --------------------------------------------------------------------------------- ##

install(TARGETS libhex-common
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)

install(DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/hex"
    FILES_MATCHING
    PATTERN "test" EXCLUDE
    PATTERN "*.h"
)

# bundled header-only libraries used by the public headers above
install(DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/../libs/png++"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/hex"
    FILES_MATCHING PATTERN "*.hpp"
)
