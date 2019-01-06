cmake_minimum_required(VERSION 3.13)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

set(libhex-db_SOURCES
    "src/db.cpp"
    "src/quantities.cpp"
    "src/utils.cpp"
    "src/var/asy.cpp"
    "src/var/ccs.cpp"
    "src/var/colls.cpp"
    "src/var/dcs.cpp"
    "src/var/ics.cpp"
    "src/var/ionamp.cpp"
    "src/var/ionf.cpp"
    "src/var/momtf.cpp"
    "src/var/scatamp-dir.cpp"
    "src/var/scatamp.cpp"
    "src/var/sdcs12.cpp"
    "src/var/sdcsE.cpp"
    "src/var/spflip.cpp"
    "src/var/stokes.cpp"
    "src/var/tcs.cpp"
    "src/var/tdcs.cpp"
    "src/var/tmat.cpp"
    "../libs/sqlitepp/binders.cpp"
    "../libs/sqlitepp/exception.cpp"
    "../libs/sqlitepp/query.cpp"
    "../libs/sqlitepp/session.cpp"
    "../libs/sqlitepp/statement.cpp"
    "../libs/sqlitepp/string.cpp"
    "../libs/sqlitepp/transaction.cpp"
)

add_library(libhex-db SHARED ${libhex-db_SOURCES})

set_target_properties(libhex-db PROPERTIES PREFIX "")

target_include_directories(libhex-db PUBLIC
    src
    ../common
    ../libs
)

target_link_libraries(libhex-db
    libhex-common
)

set(hex-db_SOURCES
    "src/ui.cpp"
)

add_executable(hex-db ${hex-db_SOURCES})

target_link_libraries(hex-db
    libhex-common
    libhex-db
    ${SQLite3_LIBRARIES}
)
