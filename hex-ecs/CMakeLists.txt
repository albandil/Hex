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

set(hex-ecs_SOURCES
    "src/amplitudes.cpp"
    "src/ang.cpp"
    "src/bspline.cpp"
    "src/gauss.cpp"
    "src/hldata.cpp"
    "src/inout.cpp"
    "src/luft.cpp"
    "src/main.cpp"
    "src/preconditioners.cpp"
    "src/radial.cpp"
    "src/slater.cpp"
    "src/solver.cpp"
    "src/factorizers/lu-lapack.cpp"
    "src/factorizers/lu-mumps.cpp"
    "src/factorizers/lu-pardiso.cpp"
    "src/factorizers/lu-scalapack.cpp"
    "src/factorizers/lu-superlu-dist.cpp"
    "src/factorizers/lu-superlu.cpp"
    "src/factorizers/lu-umfpack.cpp"
    "src/preconditioners/CGPreconditioner.cpp"
    "src/preconditioners/CoupledPreconditioner.cpp"
    "src/preconditioners/DOMPreconditioner.cpp"
    "src/preconditioners/GPUPreconditioner.cpp"
    "src/preconditioners/HybPreconditioner.cpp"
    "src/preconditioners/ILUPreconditioner.cpp"
    "src/preconditioners/KPAPreconditioner.cpp"
    "src/preconditioners/NoPreconditioner.cpp"
    "src/preconditioners/SegregatedPreconditioner.cpp"
)

add_executable(hex-ecs ${hex-ecs_SOURCES})

if(DISABLE_PARALLEL_PRECONDITION)
    target_compile_definitions(hex-ecs PRIVATE DISABLE_PARALLEL_PRECONDITION)
endif()

target_include_directories(hex-ecs PUBLIC
    src
    ../common
    ${GSL_INCLUDE_DIRS}
    ${OpenCL_INCLUDE_DIRS}
    ${MKL_INCLUDE_DIRS}
)

target_compile_options(hex-ecs PUBLIC ${OpenMP_CXX_FLAGS} ${OPTFLAGS})

target_link_libraries(hex-ecs
    libhex-common
    ${OpenCL_LIBRARIES}
    ${MKL_LIBRARIES}
    ${SCALAPACK_LIBRARIES}
)

target_link_options(hex-ecs PUBLIC ${OpenMP_CXX_FLAGS})
