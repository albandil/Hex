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

find_package(SQLite3 REQUIRED)

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
