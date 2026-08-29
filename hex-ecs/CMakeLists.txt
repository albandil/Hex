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

# The OpenCL kernel is embedded into the binary as a C byte array. It is only
# referenced from the WITH_OPENCL part of GPUPreconditioner.cpp, so the generated
# header is not needed in builds without OpenCL.
if(WITH_OPENCL)

    set(GPU_CL_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/src/preconditioners/GPUPreconditioner.cl")
    set(GPU_CL_INC    "${CMAKE_CURRENT_BINARY_DIR}/GPUPreconditioner.inc")

    add_custom_command(
        OUTPUT  "${GPU_CL_INC}"
        COMMAND ${CMAKE_COMMAND}
                -D "INPUT=${GPU_CL_SOURCE}"
                -D "OUTPUT=${GPU_CL_INC}"
                -D "VARNAME=src_preconditioners_GPUPreconditioner_cl"
                -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/embed_file.cmake"
        DEPENDS "${GPU_CL_SOURCE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/cmake/embed_file.cmake"
        COMMENT "Embedding GPUPreconditioner.cl as C byte array"
        VERBATIM
    )

    add_custom_target(hex-ecs-embed-gpu-cl DEPENDS "${GPU_CL_INC}")

endif()

set(hex-ecs_SOURCES
    "src/amplitudes.cpp"
    "src/ang.cpp"
    "src/bspline.cpp"
    "src/eigchans.cpp"
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
    "src/preconditioners/DiagPreconditioner.cpp"
    "src/preconditioners/DOMPreconditioner.cpp"
    "src/preconditioners/GPUPreconditioner.cpp"
    "src/preconditioners/HybPreconditioner.cpp"
    "src/preconditioners/ILUPreconditioner.cpp"
    "src/preconditioners/KPAPreconditioner.cpp"
    "src/preconditioners/NoPreconditioner.cpp"
    #"src/preconditioners/SegregatedPreconditioner.cpp"
)

add_executable(hex-ecs ${hex-ecs_SOURCES})

if(WITH_OPENCL)
    add_dependencies(hex-ecs hex-ecs-embed-gpu-cl)
endif()

if(DISABLE_PARALLEL_PRECONDITION)
    target_compile_definitions(hex-ecs PRIVATE DISABLE_PARALLEL_PRECONDITION)
endif()

target_include_directories(hex-ecs PUBLIC
    src
    ../common
    "${CMAKE_CURRENT_BINARY_DIR}"
    ${GSL_INCLUDE_DIRS}
    ${OpenCL_INCLUDE_DIRS}
    ${MKL_INCLUDE_DIRS}
    ${MUMPS_INCLUDE_DIRS}
    ${SUPERLU_INCLUDE_DIRS}
    ${SUPERLU_DIST_INCLUDE_DIRS}
)

target_compile_options(hex-ecs PUBLIC ${OpenMP_CXX_FLAGS} ${OPTFLAGS})

target_link_libraries(hex-ecs
    libhex-common
    ${OpenCL_LIBRARIES}
    ${MKL_LIBRARIES}
    ${MUMPS_LIBRARIES}
    ${PARDISO_LIBRARIES}
    ${SCALAPACK_LIBRARIES}
    ${SUPERLU_LIBRARIES}
    ${SUPERLU_DIST_LIBRARIES}
    OpenMP::OpenMP_CXX
)

if(WITH_MPI)
    target_link_libraries(hex-ecs
        ${SCALAPACK_LIBRARIES}
        MPI::MPI_C
    )
endif()

## --------------------------------------------------------------------------------- ##

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()

## --------------------------------------------------------------------------------- ##

install(TARGETS hex-ecs
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)
