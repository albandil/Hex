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

cmake_minimum_required(VERSION 3.20)

project(Hex)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

include(HexDependency)

# hex-ecs uses std::filesystem
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

## --------------------------------------------------------------------------------- ##

option(WITH_BOINC        "Build as a BOINC application"                           OFF)
option(WITH_CLN          "Build with CLN (Class Library for Numbers)"             OFF)
option(WITH_GINAC        "Build with GiNaC (symbolic algebra library)"            OFF)
option(WITH_HDF5         "Build with HDF5"                                        ON)
option(WITH_MKL          "Build with MKL PARDISO"                                 OFF)
option(WITH_MPI          "Build with MPI"                                         ON)
option(WITH_MUMPS        "Build with MUMPS (enables coupled preconditioner)"      OFF)
option(WITH_OPENCL       "Build with OpenCL support (enables GPU preconditioner)" OFF)
option(WITH_PARDISO      "Build with the stand-alone PARDISO library"             OFF)
option(WITH_PNG          "Build with libpng library (for debugging)"              OFF)
option(WITH_SCALAPACK    "Build with ScaLAPACK library"                           OFF)
option(WITH_SUPERLU      "Build with SuperLU library"                             OFF)
option(WITH_SUPERLU_DIST "Build with SuperLU_DIST library"                        OFF)
option(WITH_UMFPACK      "Build with UMFPACK library"                             OFF)

option(BUILD_DEVEL       "Build also unfinished programs"                         OFF)
option(BUILD_TESTING     "Make the hex-ecs regression tests available to CTest"   ON)
option(SINGLE_PRECISION  "Use single precision for real numbers"                  OFF)

option(DISABLE_PARALLEL_PRECONDITION "Avoid nested parallelism in preconditioning" ON)

## --------------------------------------------------------------------------------- ##

# The dependencies below have no find module in CMake itself, so Hex brings its own
# in cmake/. Each of them accepts <PKG>_LIBRARIES (and <PKG>_INCLUDE_DIRS) as a
# complete, ready-made answer: setting those bypasses the search, which is what the
# more baroque installations of MUMPS and PARDISO need. Otherwise the usual hints
# apply -- CMAKE_PREFIX_PATH, <PKG>_ROOT, and the <PKG>_LIBRARY cache entries that
# the search itself creates.

## --------------------------------------------------------------------------------- ##

# unconditional dependencies
find_package(OpenMP  REQUIRED)
find_package(BLAS    REQUIRED)
find_package(LAPACK  REQUIRED)
find_package(GSL     REQUIRED)
find_package(SQLite3 REQUIRED)

# Optional dependencies. An enabled WITH_* option is taken as a hard requirement,
# so that a missing library is reported here rather than at link time. All of these
# honour the usual search hints -- CMAKE_PREFIX_PATH, <PKG>_ROOT, and the cache
# entries the searches create, such as OpenCL_LIBRARY or SUPERLU_INCLUDE_DIR.
if(WITH_MPI)
    find_package(MPI REQUIRED COMPONENTS C)
endif()

if(WITH_OPENCL)
    find_package(OpenCL REQUIRED)
endif()

if(WITH_PNG)
    find_package(PNG REQUIRED)
endif()

if(WITH_BOINC)
    find_package(BOINC REQUIRED)
endif()

if(WITH_CLN)
    find_package(CLN REQUIRED)
endif()

if(WITH_GINAC)
    find_package(GINAC REQUIRED)
endif()

if(WITH_MKL)
    find_package(MKL REQUIRED)
endif()

if(WITH_MUMPS)
    find_package(MUMPS REQUIRED)
endif()

if(WITH_PARDISO)
    find_package(PARDISO REQUIRED)
endif()

if(WITH_SUPERLU)
    find_package(SUPERLU REQUIRED)
endif()

if(WITH_SUPERLU_DIST)
    find_package(SUPERLU_DIST REQUIRED)
endif()

if(WITH_UMFPACK)
    find_package(UMFPACK REQUIRED)
endif()

# after MKL, which brings a ScaLAPACK of its own that must not be mixed with the
# Netlib one
if(WITH_SCALAPACK)
    find_package(SCALAPACK REQUIRED)
endif()

# HDF5 is enabled by default and degrades gracefully: the parts of the code that
# need it are gated on HDF5_FOUND, so a missing HDF5 is not an error.
if(WITH_HDF5)
    find_package(HDF5)
endif()

## --------------------------------------------------------------------------------- ##

# Hex::config holds what applies to every target of the project: the global
# preprocessor macros and the optimization flags. Everything links it, directly or
# through libhex-common, so this is the single place to change them.

add_library(hex-config INTERFACE)
add_library(Hex::config ALIAS hex-config)

target_compile_definitions(hex-config INTERFACE
    _LONGINT                                    # use long integers
    $<$<BOOL:${SINGLE_PRECISION}>:SINGLE>       # select the floating-point type
)

# These are GCC/Clang/icpx spellings, and they are meant for the C++ sources only:
# the project also builds a C and a Fortran test program, and gfortran does not
# accept them. Unlike a plain if(CMAKE_BUILD_TYPE STREQUAL "Release"), the
# generator expression is also correct under a multi-configuration generator.
if(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang|Intel|IntelLLVM)$")

    # release requirements on IEEE conformance to enable fast math
    target_compile_options(hex-config INTERFACE
        "$<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:-fno-rounding-math>"
        "$<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:-fno-signed-zeros>"
    )

    # do not check complex overflows in FP64, where they do not occur
    if(NOT SINGLE_PRECISION)
        target_compile_options(hex-config INTERFACE
            "$<$<COMPILE_LANGUAGE:CXX>:-fcx-limited-range>"
        )
    endif()

endif()

## --------------------------------------------------------------------------------- ##

# Every dependency becomes an interface library Hex::<name> carrying its include
# directories, its libraries and the macro that switches on the code using it. The
# targets exist even when the dependency is disabled, in which case they are empty,
# so the subdirectories below link them without testing any option.

hex_dependency(blas
    ENABLED      TRUE
    DESCRIPTION  "BLAS"
    TARGETS      BLAS::BLAS
)

hex_dependency(boinc
    ENABLED      "${WITH_BOINC}"
    DEFINE       WITH_BOINC
    DESCRIPTION  "BOINC, built as a BOINC application"
    TARGETS      BOINC::BOINC
)

hex_dependency(cln
    ENABLED      "${WITH_CLN}"
    DEFINE       WITH_CLN
    DESCRIPTION  "CLN, Class Library for Numbers"
    TARGETS      CLN::CLN
)

hex_dependency(ginac
    ENABLED      "${WITH_GINAC}"
    DEFINE       WITH_GINAC
    DESCRIPTION  "GiNaC, symbolic Born amplitudes"
    TARGETS      GiNaC::GiNaC
)

hex_dependency(gsl
    ENABLED      TRUE
    DESCRIPTION  "GNU Scientific Library"
    TARGETS      GSL::gsl GSL::gslcblas
)

hex_dependency(hdf5
    ENABLED      "${HDF5_FOUND}"
    DESCRIPTION  "HDF5, needed by the hex-hdf2hdf conversion utility"
    TARGETS      HDF5::HDF5
)

hex_dependency(lapack
    ENABLED      TRUE
    DESCRIPTION  "LAPACK"
    TARGETS      LAPACK::LAPACK
)

hex_dependency(mkl
    ENABLED      "${WITH_MKL}"
    DEFINE       WITH_MKL
    DESCRIPTION  "Intel MKL PARDISO, sparse LU in hex-ecs"
    TARGETS      MKL::MKL
)

hex_dependency(mpi
    ENABLED      "${WITH_MPI}"
    DEFINE       WITH_MPI
    DESCRIPTION  "MPI, distributed runs of hex-ecs"
    TARGETS      MPI::MPI_C
)

hex_dependency(mumps
    ENABLED      "${WITH_MUMPS}"
    DEFINE       WITH_MUMPS
    DESCRIPTION  "MUMPS, distributed LU in hex-ecs and the coupled preconditioner"
    TARGETS      MUMPS::MUMPS
)

hex_dependency(opencl
    ENABLED      "${WITH_OPENCL}"
    DEFINE       WITH_OPENCL
    DESCRIPTION  "OpenCL, the GPU preconditioner of hex-ecs"
    TARGETS      OpenCL::OpenCL
)

hex_dependency(openmp
    ENABLED      TRUE
    DESCRIPTION  "OpenMP"
    TARGETS      OpenMP::OpenMP_CXX
)

hex_dependency(pardiso
    ENABLED      "${WITH_PARDISO}"
    DEFINE       WITH_PARDISO
    DESCRIPTION  "stand-alone PARDISO, sparse LU in hex-ecs"
    TARGETS      PARDISO::PARDISO
)

hex_dependency(png
    ENABLED      "${WITH_PNG}"
    DEFINE       WITH_PNG
    DESCRIPTION  "libpng, printing the matrix structure for debugging"
    TARGETS      PNG::PNG
)

hex_dependency(scalapack
    ENABLED      "${WITH_SCALAPACK}"
    DEFINE       WITH_SCALAPACK
    DESCRIPTION  "ScaLAPACK, distributed dense LU in hex-ecs"
    TARGETS      ScaLAPACK::ScaLAPACK
)

# FindSQLite3 names its target SQLite3::SQLite3 since CMake 4.3 and keeps the older
# SQLite::SQLite3 as a deprecated alias, which warns when linked. CMake 3.14 to 4.2
# offer only the old name, so take whichever is there rather than raise the required
# version to a release that has only just appeared.
if(TARGET SQLite3::SQLite3)
    set(HEX_SQLITE3_TARGET SQLite3::SQLite3)
else()
    set(HEX_SQLITE3_TARGET SQLite::SQLite3)
endif()

hex_dependency(sqlite3
    ENABLED      TRUE
    DESCRIPTION  "SQLite3"
    TARGETS      ${HEX_SQLITE3_TARGET}
)

hex_dependency(superlu
    ENABLED      "${WITH_SUPERLU}"
    DEFINE       WITH_SUPERLU
    DESCRIPTION  "SuperLU, sparse LU in hex-ecs"
    TARGETS      SuperLU::SuperLU
)

hex_dependency(superlu_dist
    ENABLED      "${WITH_SUPERLU_DIST}"
    DEFINE       WITH_SUPERLU_DIST
    DESCRIPTION  "SuperLU_DIST, distributed LU in hex-ecs"
    TARGETS      SuperLU::SuperLUDist
)

hex_dependency(umfpack
    ENABLED      "${WITH_UMFPACK}"
    DEFINE       WITH_UMFPACK
    DESCRIPTION  "UMFPACK, sparse LU in hex-ecs"
    TARGETS      UMFPACK::UMFPACK
)

## --------------------------------------------------------------------------------- ##

# standard GNU installation layout (bin/, lib/, include/hex/, share/doc/Hex/)
include(GNUInstallDirs)

# the tests are declared in hex-ecs/tests; they are never run by the build itself,
# only on demand by "ctest"
if(BUILD_TESTING)
    enable_testing()
endif()

# let the installed programs find the installed Hex shared libraries
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

## --------------------------------------------------------------------------------- ##

add_subdirectory(common)
add_subdirectory(utilities)

add_subdirectory(hex-db)
add_subdirectory(hex-dwba)
add_subdirectory(hex-ecs)

if(BUILD_DEVEL)
    add_subdirectory(hex-pecs)
endif()

## --------------------------------------------------------------------------------- ##

install(FILES README.md
    DESTINATION "${CMAKE_INSTALL_DOCDIR}"
)

install(DIRECTORY html/
    DESTINATION "${CMAKE_INSTALL_DOCDIR}/html"
)

## --------------------------------------------------------------------------------- ##

feature_summary(WHAT ENABLED_FEATURES QUIET_ON_EMPTY
    DESCRIPTION "Hex will be built against:"
)

feature_summary(WHAT DISABLED_FEATURES QUIET_ON_EMPTY
    DESCRIPTION "Hex will be built without (see the WITH_* options):"
)

## --------------------------------------------------------------------------------- ##
