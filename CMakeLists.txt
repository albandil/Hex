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

project(Hex)

enable_language(Fortran)

## --------------------------------------------------------------------------------- ##

option(WITH_BOINC        "Build as a BOINC application"                           OFF)
option(WITH_CLN          "Build with CLN (Class Library for Numbers)"             OFF)
option(WITH_GINAC        "Build with GiNaC (symbolic algebra library)"            OFF)
option(WITH_MKL          "Build with MKL PARDISO"                                 OFF)
option(WITH_MUMPS        "Build with MUMPS (enables coupled preconditioner)"      OFF)
option(WITH_OPENCL       "Build with OpenCL support (enables GPU preconditioner)" OFF)
option(WITH_PARDISO      "Build with the stand-alone PARDISO library"             OFF)
option(WITH_PNG          "Build with libpng library (for debugging)"              OFF)
option(WITH_SCALAPACK    "Build with ScaLAPACK library"                           OFF)
option(WITH_SUPERLU      "Build with SuperLU library"                             OFF)
option(WITH_SUPERLU_DIST "Build with SuperLU_DIST library"                        OFF)
option(WITH_UMFPACK      "Build with UMFPACK library"                             OFF)

option(BUILD_DEVEL       "Build also unfinished programs"                         OFF)
option(SINGLE_PRECISION  "Use single precision for real numbers"                  OFF)

option(DISABLE_PARALLEL_PRECONDITION "Avoid nested parallelism in preconditioning" ON)

## --------------------------------------------------------------------------------- ##

set(BOINC_INCLUDE_DIRS     "" CACHE STRING "BOINC include directories")
set(BOINC_LIBRARIES        "" CACHE STRING "BOINC libraries (libboinc_api, libboinc)")

set(CLN_INCLUDE_DIRS       "" CACHE STRING "CLN include directories")
set(CLN_LIBRARIES          "" CACHE STRING "CLN libraries")

set(GINAC_INCLUDE_DIRS     "" CACHE STRING "GINAC include directories")
set(GINAC_LIBRARIES        "" CACHE STRING "GINAC libraries")

set(GSL_INCLUDE_DIRS       "" CACHE STRING "GNU Scientific Library include directories")
set(GSL_LIBRARIES          "" CACHE STRING "GNU Scientific Library libraries")

set(OpenCL_INCLUDE_DIRS    "" CACHE STRING "OpenCL include directories")
set(OpenCL_LIBRARIES       "" CACHE STRING "OpenCL libraries")

set(MKL_INCLUDE_DIRS       "" CACHE STRING "MKL include directories (not used)")
set(MKL_LIBRARIES          "" CACHE STRING "MKL libraries")

set(MUMPS_INCLUDE_DIRS     "" CACHE STRING "MUMPS include directories")
set(MUMPS_LIBRARIES        "" CACHE STRING "MUMPS libraries")

set(PARDISO_INCLUDE_DIRS   "" CACHE STRING "Stand-alone PARDISO include directories (not used)")
set(PARDISO_LIBRARIES      "" CACHE STRING "Stand-alone PARDISO libraries")

set(PNG_INCLUDE_DIRS       "" CACHE STRING "libpng include directories")
set(PNG_LIBRARIES          "" CACHE STRING "libpng libraries")

set(SQLite3_INCLUDE_DIRS   "" CACHE STRING "SQLite3 include directories")
set(SQLite3_LIBRARIES      "" CACHE STRING "SQLite3 libraries")

set(SCALAPACK_INCLUDE_DIRS "" CACHE STRING "SCALAPACK include directories (not used)")
set(SCALAPACK_LIBRARIES    "" CACHE STRING "SCALAPACK libraries")

set(UMFPACK_INCLUDE_DIRS   "" CACHE STRING "UMFPACK include directories")
set(UMFPACK_LIBRARIES      "" CACHE STRING "UMFPACK libraries")

## --------------------------------------------------------------------------------- ##

if(WITH_CLN AND NOT CLN_LIBRARIES)
    set(CLN_LIBRARIES cln)
endif()

if(WITH_GINAC AND NOT GINAC_LIBRARIES)
    set(GINAC_LIBRARIES ginac)
endif()

if(NOT GSL_LIBRARIES)
    set(GSL_LIBRARIES gsl gslcblas)
endif()

if(NOT SQLite3_LIBRARIES)
    set(SQLite3_LIBRARIES sqlite3)
endif()

## --------------------------------------------------------------------------------- ##

find_package(OpenMP REQUIRED)
find_package(BLAS   REQUIRED)
find_package(LAPACK REQUIRED)

find_package(MPI)
find_package(HDF5)

## --------------------------------------------------------------------------------- ##

# use long integers
add_definitions(-D_LONGINT)

# release requirements on IEEE conformance to enable fast math
if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
    set(OPTFLAGS ${OPTFLAGS} -funsafe-math-optimizations -fno-rounding-math -fno-signed-zeros)
endif()

# do not check complex overflows in FP64, where they do not occur
if(NOT SINGLE_PRECISION)
    set(OPTFLAGS ${OPTFLAGS} -fcx-limited-range)
endif()

if(WITH_BOINC)
    message(STATUS "Building for BOINC")
    message(STATUS "    includes:  ${BOINC_INCLUDE_DIRS}")
    message(STATUS "    libraries: ${BOINC_LIBRARIES}")
    add_definitions(-DWITH_CLN)
endif()

if(WITH_CLN)
    message(STATUS "Building with CLN")
    message(STATUS "    includes:  ${CLN_INCLUDE_DIRS}")
    message(STATUS "    libraries: ${CLN_LIBRARIES}")
    add_definitions(-DWITH_CLN)
endif()

if(WITH_GINAC)
    message(STATUS "Building with GiNaC")
    message(STATUS "    includes:  ${GINAC_INCLUDE_DIRS}")
    message(STATUS "    libraries: ${GINAC_LIBRARIES}")
    add_definitions(-DWITH_GINAC)
endif()

if(WITH_OPENCL)
    message(STATUS "Building with OpenCL")
    message(STATUS "    includes:  ${OpenCL_INCLUDE_DIRS}")
    message(STATUS "    libraries: ${OpenCL_LIBRARIES}")
    add_definitions(-DWITH_OPENCL)
endif()

if(WITH_PNG)
    if(NOT PNG_LIBRARIES)
        set(PNG_LIBRARIES libpng)
    endif()
    message(STATUS "Building with PNG")
    message(STATUS "    includes:  ${PNG_INCLUDE_DIRS}")
    message(STATUS "    libraries: ${PNG_LIBRARIES}")
    add_definitions(-DWITH_PNG)
endif()

if(WITH_MKL)
    message(STATUS "Building with MKL")
    message(STATUS "    includes:  ${MKL_INCLUDE_DIRS}")
    message(STATUS "    libraries: ${MKL_LIBRARIES}")
    add_definitions(-DWITH_MKL)
endif()

if(WITH_MUMPS)
    message(STATUS "Building with MUMPS")
    message(STATUS "    includes:  ${MUMPS_INCLUDE_DIRS}")
    message(STATUS "    libraries: ${MUMPS_LIBRARIES}")
    add_definitions(-DWITH_MUMPS)
endif()

if(WITH_SCALAPACK)
    message(STATUS "Building with ScaLAPACK")
    message(STATUS "    includes:  ${SCALAPACK_INCLUDE_DIRS}")
    message(STATUS "    libraries: ${SCALAPACK_LIBRARIES}")
    add_definitions(-DWITH_SCALAPACK)
endif()

if(WITH_UMFPACK)
    message(STATUS "Building with UMFPACK")
    message(STATUS "    includes:  ${UMFPACK_INCLUDE_DIRS}")
    message(STATUS "    libraries: ${UMFPACK_LIBRARIES}")
    add_definitions(-DWITH_UMFPACK)
endif()

if(MPI_FOUND)
    message(STATUS "Building with MPI")
    add_definitions(-DWITH_MPI)
    set(CMAKE_CXX_COMPILER "${MPI_CXX_COMPILER}")
endif()

## --------------------------------------------------------------------------------- ##

add_subdirectory(common)
add_subdirectory(utilities)

add_subdirectory(hex-db)
add_subdirectory(hex-dwba)
add_subdirectory(hex-ecs)

if(WITH_DEVEL)
    add_subdirectory(hex-pecs)
endif()

## --------------------------------------------------------------------------------- ##
