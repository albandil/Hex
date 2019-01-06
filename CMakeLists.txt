cmake_minimum_required(VERSION 3.13)

project(Hex)

# ----------------------------------------------------------------------------------------------- #

option(WITH_CLN          "Build with CLN (Class Library for Numbers)"             OFF)
option(WITH_GINAC        "Build with GiNaC (symbolic algebra library)"            OFF)
option(WITH_MKL          "Build with MKL BLAS, LAPACK and PARDISO"                OFF)
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

# ----------------------------------------------------------------------------------------------- #

set(CLN_INCLUDE_DIRS     "" CACHE STRING "CLN include directories")
set(CLN_LIBRARIES        "" CACHE STRING "CLN libraries")

set(GINAC_INCLUDE_DIRS   "" CACHE STRING "GINAC include directories")
set(GINAC_LIBRARIES      "" CACHE STRING "GINAC libraries")

set(GSL_INCLUDE_DIRS     "" CACHE STRING "GNU Scientific Library include directories")
set(GSL_LIBRARIES        "" CACHE STRING "GNU Scientific Library libraries")

set(OpenCL_INCLUDE_DIRS  "" CACHE STRING "OpenCL include directories")
set(OpenCL_LIBRARIES     "" CACHE STRING "OpenCL libraries")

set(PARDISO_INCLUDE_DIRS "" CACHE STRING "Stand-alone PARDISO include directories")
set(PARDISO_LIBRARIES    "" CACHE STRING "Stand-alone PARDISO libraries")

set(PNG_INCLUDE_DIRS     "" CACHE STRING "libpng include directories")
set(PNG_LIBRARIES        "" CACHE STRING "libpng libraries")

set(SQLite3_INCLUDE_DIRS "" CACHE STRING "SQLite3 include directories")
set(SQLite3_LIBRARIES    "" CACHE STRING "SQLite3 libraries")

set(UMFPACK_INCLUDE_DIRS "" CACHE STRING "UMFPACK include directories")
set(UMFPACK_LIBRARIES    "" CACHE STRING "UMFPACK libraries")

# ----------------------------------------------------------------------------------------------- #

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

# ----------------------------------------------------------------------------------------------- #

find_package(OpenMP REQUIRED)
find_package(BLAS   REQUIRED)
find_package(LAPACK REQUIRED)

find_package(MPI)
find_package(HDF5)

# ----------------------------------------------------------------------------------------------- #

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

if(WITH_CLN)
    add_definitions(-DWITH_CLN)
endif()

if(WITH_GINAC)
    add_definitions(-DWITH_GINAC)
endif()

if(WITH_OPENCL)
    add_definitions(-DWITH_OPENCL)
endif()

if(WITH_PNG)
    add_definitions(-DWITH_PNG)
endif()

if(WITH_UMFPACK)
    add_definitions(-DWITH_UMFPACK)
endif()

if(MPI_FOUND)
    add_definitions(-DWITH_MPI)
    set(CMAKE_CXX_COMPILER "${MPI_CXX_COMPILER}")
endif()

# ----------------------------------------------------------------------------------------------- #

add_subdirectory(common)
add_subdirectory(utilities)

add_subdirectory(hex-db)
add_subdirectory(hex-dwba)
add_subdirectory(hex-ecs)

if(WITH_DEVEL)
    add_subdirectory(hex-pecs)
endif()

# ----------------------------------------------------------------------------------------------- #
