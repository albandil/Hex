cmake_minimum_required(VERSION 3.13)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

set(hex-dwba_SOURCES
    "src/dwba.cpp"
    "src/potential.cpp"
    "src/ui.cpp"
    "src/wave_distort.cpp"
)

add_executable(hex-dwba ${hex-dwba_SOURCES})

target_include_directories(hex-dwba PUBLIC
    src
    ../common
)

target_compile_options(hex-dwba PUBLIC ${OpenMP_CXX_FLAGS})

target_link_options(hex-dwba PUBLIC ${OpenMP_CXX_FLAGS})

target_link_libraries(hex-dwba
    libhex-common
)
