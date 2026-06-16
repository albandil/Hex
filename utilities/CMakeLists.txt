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

add_executable(hex-chebeval chebeval/chebeval.cpp)
add_executable(hex-f-pattern f-pattern/hex-f-pattern.cpp)
#add_executable(hex-hdf2hdf hdf2hdf/hdf2hdf.cpp)
add_executable(hex-hdf2txt hdf2txt/hdf2txt.cpp)
add_executable(hex-special special/hex-special.cpp)

target_include_directories(hex-chebeval PRIVATE "${CMAKE_CURRENT_LIST_DIR}/../common")
target_include_directories(hex-f-pattern PRIVATE "${CMAKE_CURRENT_LIST_DIR}/../common")
#target_include_directories(hex-hdf2hdf PRIVATE "${CMAKE_CURRENT_LIST_DIR}/../common")
target_include_directories(hex-hdf2txt PRIVATE "${CMAKE_CURRENT_LIST_DIR}/../common")
target_include_directories(hex-special PRIVATE "${CMAKE_CURRENT_LIST_DIR}/../common")

target_link_libraries(hex-chebeval PRIVATE libhex-common)
target_link_libraries(hex-f-pattern PRIVATE libhex-common)
#target_link_libraries(hex-hdf2hdf PRIVATE libhex-common)
target_link_libraries(hex-hdf2txt PRIVATE libhex-common)
target_link_libraries(hex-special PRIVATE libhex-common)
