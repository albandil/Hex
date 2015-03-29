CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
CMAKE_POLICY(SET CMP0037 OLD)

INCLUDE(CheckFunctionExists)

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

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_LIST_DIR}/cmake)

FIND_PACKAGE(CLN)
FIND_PACKAGE(GiNaC)
FIND_PACKAGE(FFTW3)
FIND_PACKAGE(GSL)
FIND_PACKAGE(MPI)
FIND_PACKAGE(OpenCL)
FIND_PACKAGE(PNG)
FIND_PACKAGE(SQLite3)
FIND_PACKAGE(SuiteSparse)

FIND_PACKAGE(OpenBLAS)
IF(OPENBLAS_FOUND)
    SET(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
    CHECK_FUNCTION_EXISTS("zsbmv_" ZSBMV_EXISTS)
    IF(NOT ZSBMV_EXISTS)
        SET(LAPACK_FOUND FALSE)
    ENDIF(NOT ZSBMV_EXISTS)
ENDIF(OPENBLAS_FOUND)

IF(NOT OPENBLAS_FOUND)
    FIND_PACKAGE(BLAS)
    FIND_PACKAGE(LAPACK)
    SET(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
    CHECK_FUNCTION_EXISTS("zsbmv_" ZSBMV_EXISTS)
    IF(NOT ZSBMV_EXISTS)
        MESSAGE("-- Library ${LAPACK_LIBRARIES} does not define required routine ZSBMV.")
        SET(LAPACK_FOUND FALSE)
    ENDIF(NOT ZSBMV_EXISTS)
ENDIF(NOT OPENBLAS_FOUND)

SET(CMAKE_REQUIRED_LIBRARIES m)
CHECK_FUNCTION_EXISTS("hypot" HYPOT_EXISTS)
IF(NOT HYPOT_EXISTS)
    MESSAGE(FATAL_ERROR "Math library doesn't define C99 function \"hypot\".")
ENDIF(NOT HYPOT_EXISTS)

MESSAGE("")

# ----------------------------------------------------------------------------
# Hex-DB

IF(GINAC_FOUND AND SQLITE3_FOUND AND GSL_FOUND AND FFTW3_FOUND)

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
    
ELSE()

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

ENDIF()
MESSAGE("")

# ----------------------------------------------------------------------------
# Hex-ECS

IF(UMFPACK_FOUND AND GSL_FOUND AND BLAS_FOUND AND LAPACK_FOUND AND FFTW3_FOUND)

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
    
    TARGET_INCLUDE_DIRECTORIES(${HEX_ECS} PUBLIC ${UMFPACK_INCLUDE_DIR})
    TARGET_INCLUDE_DIRECTORIES(${HEX_ECS} PUBLIC ${FFTW3_INCLUDE_DIR})
    TARGET_INCLUDE_DIRECTORIES(${HEX_ECS} PUBLIC ${BLAS_INCLUDE_DIR})
    TARGET_INCLUDE_DIRECTORIES(${HEX_ECS} PUBLIC ${LAPACK_INCLUDE_DIR})
    TARGET_INCLUDE_DIRECTORIES(${HEX_ECS} PUBLIC ${GSL_INCLUDE_DIR})
    
    IF(UMFPACK_STAND_ALONE)
        TARGET_LINK_LIBRARIES(${HEX_ECS} ${UMFPACK_LIBRARY})
    ENDIF(UMFPACK_STAND_ALONE)
    IF(UMFPACK_WITH_SUITESPARSE)
        TARGET_LINK_LIBRARIES(${HEX_ECS} ${UMFPACK_LIBRARY} ${CHOLMOD_LIBRARY} ${AMD_LIBRARY} ${CAMD_LIBRARY} ${COLAMD_LIBRARY} ${CCOLAMD_LIBRARY} ${SUITESPARSE_CONFIG_LIBRARY} ${METIS_LIBRARY} m)
    ENDIF(UMFPACK_WITH_SUITESPARSE)
    TARGET_LINK_LIBRARIES(${HEX_ECS} ${FFTW3_LIBRARY})
    TARGET_LINK_LIBRARIES(${HEX_ECS} ${BLAS_LIBRARIES})
    TARGET_LINK_LIBRARIES(${HEX_ECS} ${LAPACK_LIBRARIES})
    TARGET_LINK_LIBRARIES(${HEX_ECS} ${GSL_LIBRARY} ${GSLCBLAS_LIBRARY})
    
    IF(OPENCL_FOUND)
        TARGET_INCLUDE_DIRECTORIES(${HEX_ECS} PUBLIC ${OPENCL_INCLUDE_DIR})
        TARGET_LINK_LIBRARIES(${HEX_ECS} ${OPENCL_LIBRARIES})
        MESSAGE("    - with OpenCL")
    ELSE(OPENCL_FOUND)
        MESSAGE("    - without OpenCL")
        SET(HEX_ECS_CXX_COMPILE_FLAGS ${HEX_ECS_CXX_COMPILE_FLAGS} "-DNO_OPENCL")
    ENDIF(OPENCL_FOUND)
    
    IF(MPI_FOUND)
        SET(HEX_ECS_CXX_COMPILE_FLAGS ${HEX_ECS_CXX_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS})
        SET(HEX_ECS_CXX_LINK_FLAGS ${HEX_ECS_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS})
        TARGET_INCLUDE_DIRECTORIES(${HEX_ECS} PUBLIC ${MPI_CXX_INCLUDE_PATH})
        TARGET_LINK_LIBRARIES(${HEX_ECS} ${MPI_CXX_LIBRARIES})
        MESSAGE("    - with MPI")
    ELSE(MPI_FOUND)
        MESSAGE("    - without MPI")
        SET(HEX_ECS_CXX_COMPILE_FLAGS ${HEX_ECS_CXX_COMPILE_FLAGS} "-DNO_MPI")
    ENDIF(MPI_FOUND)
    
    IF(PNG_FOUND)
        SET(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${PNG_COMPILE_FLAGS})
        SET(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${PNG_LINK_FLAGS})
        TARGET_INCLUDE_DIRECTORIES(${HEX_ECS} PUBLIC ${PNG_INCLUDE_PATH} libs)
        TARGET_LINK_LIBRARIES(${HEX_ECS} ${PNG_LIBRARIES})
        MESSAGE("    - with PNG")
    ELSE(PNG_FOUND)
        MESSAGE("    - without PNG")
        SET(HEX_ECS_CXX_COMPILE_FLAGS ${HEX_ECS_CXX_COMPILE_FLAGS} "-DNO_PNG")
    ENDIF(PNG_FOUND)
    
    SET_TARGET_PROPERTIES(${HEX_ECS} PROPERTIES COMPILE_FLAGS "${HEX_ECS_CXX_COMPILE_FLAGS}")
    SET_TARGET_PROPERTIES(${HEX_ECS} PROPERTIES LINK_FLAGS "${HEX_ECS_CXX_LINK_FLAGS}")
    
ELSE()

    MESSAGE("Hex-ECS will NOT be built")

    IF(NOT UMFPACK_FOUND)
        MESSAGE("    - missing UMFPACK !")
    ENDIF(NOT UMFPACK_FOUND)

    IF(NOT OPENBLAS_FOUND)
        MESSAGE("    - missing OpenBLAS !")
    ENDIF(NOT OPENBLAS_FOUND)
    
    IF(NOT FFTW3_FOUND)
        MESSAGE("    - missing FFTW3 !")
    ENDIF(NOT FFTW3_FOUND)
    
    IF(NOT GSL_FOUND)
        MESSAGE("    - missing GSL !")
    ENDIF(NOT GSL_FOUND)

ENDIF()
MESSAGE("")

# ----------------------------------------------------------------------------
# Hex-DWBA

IF(GSL_FOUND AND FFTW3_FOUND AND CLN_FOUND)

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

ELSE()

    MESSAGE("Hex-DWBA will NOT be built")
    
    IF(NOT GSL_FOUND)
        MESSAGE("    - missing GSL !")
    ENDIF(NOT GSL_FOUND)
    
    IF(NOT FFTW3_FOUND)
        MESSAGE("    - missing FFTW3 !")
    ENDIF(NOT FFTW3_FOUND)
    
    IF(NOT CLN_FOUND)
        MESSAGE("    - missing CLN !")
    ENDIF(NOT CLN_FOUND)
    
ENDIF()

MESSAGE("")
