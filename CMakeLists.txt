CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
CMAKE_POLICY(SET CMP0037 OLD)

INCLUDE(FindPackageHandleStandardArgs)

# ----------------------------------------------------------------------------
# The HEX project.

PROJECT(Hex)
SET(Hex_VERSION_MAJOR 1)
SET(Hex_VERSION_MINOR 6)

SET(HEX_DB hex-db)
SET(HEX_ECS hex-ecs)
SET(HEX_DWBA hex-dwba)

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_LIST_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_LIST_DIR}/bin)
SET(CMAKE_CXX_FLAGS "-std=c++11 -fopenmp -march=native -mtune=native -Wall -DNDEBUG -Isrc")

# ----------------------------------------------------------------------------
# Check available libraries.

FIND_PATH(CLN_INCLUDE_DIR cln/cln.h PATHS ${CLN_INCLUDE_DIRS})
FIND_LIBRARY(CLN_LIBRARY cln ${CLN_LIBRARIES})
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CLN REQUIRED_VARS CLN_LIBRARY CLN_INCLUDE_DIR)

FIND_PATH(GINAC_INCLUDE_DIR ginac/ginac.h PATHS ${GINAC_INCLUDE_DIRS})
FIND_LIBRARY(GINAC_LIBRARY ginac ${GINAC_LIBRARIES})
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GiNaC REQUIRED_VARS GINAC_LIBRARY GINAC_INCLUDE_DIR)

FIND_PATH(SQLITE3_INCLUDE_DIR sqlite3.h PATHS ${SQLITE3_INCLUDE_DIRS})
FIND_LIBRARY(SQLITE3_LIBRARY sqlite3 ${SQLITE3_LIBRARIES})
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SQLITE3 REQUIRED_VARS SQLITE3_LIBRARY SQLITE3_INCLUDE_DIR)

FIND_PATH(GSL_INCLUDE_DIR gsl_sf.h ${GSL_INCLUDE_DIRS} /usr/include/gsl)
FIND_LIBRARY(GSL_LIBRARY gsl ${GSL_LIBRARIES})
FIND_LIBRARY(GSLCBLAS_LIBRARY gslcblas ${GSL_LIBRARIES})
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GSL REQUIRED_VARS GSL_LIBRARY GSL_INCLUDE_DIR)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GSLCBLAS REQUIRED_VARS GSLCBLAS_LIBRARY)

FIND_PATH(UMFPACK_INCLUDE_DIR umfpack.h PATHS ${UMFPACK_INCLUDE_DIRS} /usr/include/suitesparse)
FIND_LIBRARY(UMFPACK_LIBRARY umfpack ${UMFPACK_LIBRARIES})
FIND_PACKAGE_HANDLE_STANDARD_ARGS(UMFPACK REQUIRED_VARS UMFPACK_LIBRARY UMFPACK_INCLUDE_DIR)

FIND_PATH(FFTW3_INCLUDE_DIR fftw3.h ${FFTW3_INCLUDE_DIRS})
FIND_LIBRARY(FFTW3_LIBRARY fftw3 ${FFTW3_LIBRARIES})
FIND_PACKAGE_HANDLE_STANDARD_ARGS(FFTW3 REQUIRED_VARS FFTW3_LIBRARY FFTW3_INCLUDE_DIR)

FIND_LIBRARY(OPENBLAS_LIBRARY openblas ${OPENBLAS_LIBRARIES})
FIND_PACKAGE_HANDLE_STANDARD_ARGS(OpenBLAS REQUIRED_VARS OPENBLAS_LIBRARY)

FIND_PACKAGE(MPI)

FIND_PACKAGE(PNG)

MESSAGE("")

# ----------------------------------------------------------------------------
# Hex-DB

IF (GINAC_FOUND AND SQLITE3_FOUND AND GSL_FOUND AND FFTW3_FOUND)

    MESSAGE("Hex-DB will be built")

    ADD_LIBRARY(sqlitepp SHARED
        libs/sqlitepp/binders.cpp
        libs/sqlitepp/exception.cpp
        libs/sqlitepp/query.cpp
        libs/sqlitepp/session.cpp
        libs/sqlitepp/statement.cpp
        libs/sqlitepp/string.cpp
        libs/sqlitepp/transaction.cpp
    )
    
    ADD_LIBRARY(${HEX_DB}-lib SHARED
        src/arrays.cpp
        src/born.cpp
        src/hdffile.cpp
        src/special.cpp
        src/vec3d.cpp
        src/version.cpp
        hex-db/src/hex-db.cpp
        hex-db/src/variables.cpp
        hex-db/src/var/asy.cpp
        hex-db/src/var/bornf.cpp
        hex-db/src/var/ccs.cpp
        hex-db/src/var/colls.cpp
        hex-db/src/var/dcs.cpp
        hex-db/src/var/ics.cpp
        hex-db/src/var/ionamp.cpp
        hex-db/src/var/ionf.cpp
        hex-db/src/var/momtf.cpp
        hex-db/src/var/scatamp-dir.cpp
        hex-db/src/var/scatamp.cpp
        hex-db/src/var/stokes.cpp
        hex-db/src/var/tcs.cpp
        hex-db/src/var/tdcs.cpp
        hex-db/src/var/tmat.cpp
        hex-db/src/var/tmatb.cpp
        hex-db/src/var/xcs.cpp
    )
    
    ADD_EXECUTABLE(${HEX_DB}-exe hex-db/src/ui.cpp)
    
    SET_TARGET_PROPERTIES(${HEX_DB}-lib PROPERTIES OUTPUT_NAME ${HEX_DB})
    SET_TARGET_PROPERTIES(${HEX_DB}-lib PROPERTIES COMPILE_FLAGS "-DNO_MPI -DNO_LAPACK")
    
    SET_TARGET_PROPERTIES(${HEX_DB}-exe PROPERTIES OUTPUT_NAME ${HEX_DB})
    SET_TARGET_PROPERTIES(${HEX_DB}-exe PROPERTIES COMPILE_FLAGS "-DNO_MPI -DNO_LAPACK")
    
    TARGET_INCLUDE_DIRECTORIES(${HEX_DB}-lib PUBLIC libs)
    TARGET_INCLUDE_DIRECTORIES(${HEX_DB}-lib PUBLIC ${GINAC_INCLUDE_DIR})
    TARGET_INCLUDE_DIRECTORIES(${HEX_DB}-lib PUBLIC ${FFTW3_INCLUDE_DIR})
    TARGET_INCLUDE_DIRECTORIES(${HEX_DB}-lib PUBLIC ${GSL_INCLUDE_DIR})
    TARGET_INCLUDE_DIRECTORIES(${HEX_DB}-lib PUBLIC ${SQLITE3_INCLUDE_DIR})
    
    TARGET_LINK_LIBRARIES(${HEX_DB}-lib ${GINAC_LIBRARY})
    TARGET_LINK_LIBRARIES(${HEX_DB}-lib ${FFTW3_LIBRARY})
    TARGET_LINK_LIBRARIES(${HEX_DB}-lib ${GSL_LIBRARY} ${GSLCBLAS_LIBRARY})
    TARGET_LINK_LIBRARIES(${HEX_DB}-lib ${SQLITE3_LIBRARY})
    TARGET_LINK_LIBRARIES(${HEX_DB}-lib sqlitepp)
    
    TARGET_LINK_LIBRARIES(${HEX_DB}-exe ${HEX_DB}-lib)
    
ELSE ()

    MESSAGE("Hex-DB will NOT be built")

    IF (NOT GINAC_FOUND)
        MESSAGE("    - missing GINAC !")
    ENDIF(NOT GINAC_FOUND)
    
    IF (NOT FFTW3_FOUND)
        MESSAGE("    - missing FFTW3 !")
    ENDIF(NOT FFTW3_FOUND)
    
    IF (NOT SQLITE3_FOUND)
        MESSAGE("    - missing SQLite3 !")
    ENDIF(NOT SQLITE3_FOUND)
    
    IF (NOT GSL_FOUND)
        MESSAGE("    - missing GSL !")
    ENDIF(NOT GSLFOUND)

ENDIF ()
MESSAGE("")

# ----------------------------------------------------------------------------
# Hex-ECS

IF (UMFPACK_FOUND AND GSL_FOUND AND OPENBLAS_FOUND AND FFTW3_FOUND)

    MESSAGE("Hex-ECS will be built")

    ADD_EXECUTABLE(${HEX_ECS}
        src/arrays.cpp
        src/hdffile.cpp
        src/hydrogen.cpp
        src/matrix.cpp
        src/special.cpp
        src/version.cpp
        hex-ecs/src/amplitudes.cpp
        hex-ecs/src/bspline.cpp
        hex-ecs/src/gauss.cpp
        hex-ecs/src/io.cpp
        hex-ecs/src/main.cpp
        hex-ecs/src/radial.cpp
        hex-ecs/src/slater.cpp
        hex-ecs/src/preconditioners/CGPreconditioner.cpp
        hex-ecs/src/preconditioners/GPUPreconditioner.cpp
        hex-ecs/src/preconditioners/ILUPreconditioner.cpp
        hex-ecs/src/preconditioners/KPAPreconditioner.cpp
        hex-ecs/src/preconditioners/NoPreconditioner.cpp
    )
    
    SET(HEX_ECS_CXX_COMPILE_FLAGS ${HEX_ECS_CXX_COMPILE_FLAGS} "-DNO_OPENCL")
    
    TARGET_INCLUDE_DIRECTORIES(${HEX_ECS} PUBLIC ${UMFPACK_INCLUDE_DIR})
    TARGET_INCLUDE_DIRECTORIES(${HEX_ECS} PUBLIC ${FFTW3_INCLUDE_DIR})
    TARGET_INCLUDE_DIRECTORIES(${HEX_ECS} PUBLIC ${GSL_INCLUDE_DIR})
    
    TARGET_LINK_LIBRARIES(${HEX_ECS} ${UMFPACK_LIBRARY})
    TARGET_LINK_LIBRARIES(${HEX_ECS} ${FFTW3_LIBRARY})
    TARGET_LINK_LIBRARIES(${HEX_ECS} ${OPENBLAS_LIBRARY})
    TARGET_LINK_LIBRARIES(${HEX_ECS} ${GSL_LIBRARY} ${GSLCBLAS_LIBRARY})
    
    IF (MPI_FOUND)
        SET(HEX_ECS_CXX_COMPILE_FLAGS ${HEX_ECS_CXX_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS})
        SET(HEX_ECS_CXX_LINK_FLAGS ${HEX_ECS_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS})
        TARGET_INCLUDE_DIRECTORIES(${HEX_ECS} PUBLIC ${MPI_CXX_INCLUDE_PATH})
        TARGET_LINK_LIBRARIES(${HEX_ECS} ${MPI_CXX_LIBRARIES})
        MESSAGE("    - building with MPI")
    ELSE (MPI_FOUND)
        MESSAGE("    - building without MPI")
        SET(HEX_ECS_CXX_COMPILE_FLAGS ${HEX_ECS_CXX_COMPILE_FLAGS} "-DNO_MPI")
    ENDIF (MPI_FOUND)
    
    IF (PNG_FOUND)
        SET(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${PNG_COMPILE_FLAGS})
        SET(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${PNG_LINK_FLAGS})
        TARGET_INCLUDE_DIRECTORIES(${HEX_ECS} PUBLIC ${PNG_INCLUDE_PATH} libs)
        TARGET_LINK_LIBRARIES(${HEX_ECS} ${PNG_LIBRARIES})
        MESSAGE("    - building with PNG")
    ELSE (PNG_FOUND)
        MESSAGE("    - building without PNG")
        SET(HEX_ECS_CXX_COMPILE_FLAGS ${HEX_ECS_CXX_COMPILE_FLAGS} "-DNO_PNG")
    ENDIF (PNG_FOUND)
    
    SET_TARGET_PROPERTIES(${HEX_ECS} PROPERTIES COMPILE_FLAGS "${HEX_ECS_CXX_COMPILE_FLAGS}")
    SET_TARGET_PROPERTIES(${HEX_ECS} PROPERTIES LINK_FLAGS "${HEX_ECS_CXX_LINK_FLAGS}")
    
ELSE ()

    MESSAGE("Hex-DWBA will NOT be built")

    IF (NOT UMFPACK_FOUND)
        MESSAGE("    - missing UMFPACK !")
    ENDIF(NOT UMFPACK_FOUND)

    IF (NOT OPENBLAS_FOUND)
        MESSAGE("    - missing OpenBLAS !")
    ENDIF(NOT OPENBLAS_FOUND)
    
    IF (NOT FFTW3_FOUND)
        MESSAGE("    - missing FFTW3 !")
    ENDIF(NOT FFTW3_FOUND)
    
    IF (NOT GSL_FOUND)
        MESSAGE("    - missing GSL !")
    ENDIF(NOT GSL_FOUND)

ENDIF ()
MESSAGE("")

# ----------------------------------------------------------------------------
# Hex-DWBA

IF (GSL_FOUND AND FFTW3_FOUND AND CLN_FOUND)

    MESSAGE("Hex-DWBA will be built")
    
    ADD_EXECUTABLE(${HEX_DWBA}
        src/arrays.cpp
        src/hdffile.cpp
        src/hydrogen.cpp
        src/special.cpp
        src/symbolic.cpp
        src/version.cpp
        hex-dwba/src/dwba.cpp
        hex-dwba/src/potential.cpp
        hex-dwba/src/pwba.cpp
        hex-dwba/src/ui.cpp
        hex-dwba/src/wave_distort.cpp
    )
    
    SET_TARGET_PROPERTIES(${HEX_DWBA} PROPERTIES COMPILE_FLAGS "-DNO_MPI -DNO_LAPACK")
    
    TARGET_INCLUDE_DIRECTORIES(${HEX_DWBA} PUBLIC ${CLN_INCLUDE_DIR})
    TARGET_INCLUDE_DIRECTORIES(${HEX_DWBA} PUBLIC ${FFTW3_INCLUDE_DIR})
    TARGET_INCLUDE_DIRECTORIES(${HEX_DWBA} PUBLIC ${GSL_INCLUDE_DIR})
    
    TARGET_LINK_LIBRARIES(${HEX_DWBA} ${CLN_LIBRARY})
    TARGET_LINK_LIBRARIES(${HEX_DWBA} ${FFTW3_LIBRARY})
    TARGET_LINK_LIBRARIES(${HEX_DWBA} ${GSL_LIBRARY} ${GSLCBLAS_LIBRARY})

ELSE ()

    MESSAGE("Hex-DWBA will NOT be built")
    
    IF (NOT GSL_FOUND)
        MESSAGE("    - missing GSL !")
    ENDIF(NOT GSL_FOUND)
    
    IF (NOT FFTW3_FOUND)
        MESSAGE("    - missing FFTW3 !")
    ENDIF(NOT FFTW3_FOUND)
    
    IF (NOT CLN_FOUND)
        MESSAGE("    - missing CLN !")
    ENDIF(NOT CLN_FOUND)
    
ENDIF ()

MESSAGE("")
