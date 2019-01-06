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

cmake_minimum_required(VERSION 3.13)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

set(libhex-common_SOURCES
    "hex-arrays.cpp"
    "hex-blas.cpp"
    "hex-born.cpp"
    "hex-hdffile.cpp"
    "hex-hydrogen.cpp"
    "hex-luft.cpp"
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

target_compile_options(libhex-common PUBLIC ${OpenMP_CXX_FLAGS} ${OPTFLAGS})

target_include_directories(libhex-common PUBLIC
    ../libs
    ${CLN_INCLUDE_DIRS}
    ${GSL_INCLUDE_DIRS}
    ${GINAC_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
    ${MUMPS_INCLUDE_DIRS}
    ${UMFPACK_INCLUDE_DIRS}
)

target_link_libraries(libhex-common
    ${CLN_LIBRARIES}
    ${GINAC_LIBRARIES}
    ${GSL_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${MUMPS_LIBRARIES}
    ${UMFPACK_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
)
