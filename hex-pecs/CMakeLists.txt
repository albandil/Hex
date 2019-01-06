cmake_minimum_required(VERSION 3.13)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

set(hex-pecs_SOURCES
    "src/amplitudes.cpp"
    "src/ang.cpp"
    "src/io.cpp"
    "src/main.cpp"
    "src/matops.cpp"
    "src/numerov.cpp"
    "src/os.cpp"
    "src/radial.cpp"
)

add_executable(hex-pecs ${hex-pecs_SOURCES})

target_include_directories(hex-pecs PUBLIC
    src
    ../common
)

target_link_libraries(hex-pecs
    libhex-common
)
