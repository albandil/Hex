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

target_include_directories(hex-ecs PUBLIC
    src
    ../common
    ${GSL_INCLUDE_DIRS}
    ${OpenCL_INCLUDE_DIRS}
)

target_compile_options(hex-ecs PUBLIC ${OpenMP_CXX_FLAGS} ${OPTFLAGS})

target_link_libraries(hex-ecs
    libhex-common
    ${OpenCL_LIBRARIES}
)

target_link_options(hex-ecs PUBLIC ${OpenMP_CXX_FLAGS})
