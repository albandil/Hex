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
    ${UMFPACK_INCLUDE_DIRS}
)

target_link_libraries(libhex-common
    ${CLN_LIBRARIES}
    ${GINAC_LIBRARIES}
    ${GSL_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${UMFPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
)
