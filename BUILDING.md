# Building Hex

It is recommended to use CMake for building the Hex source tree. This will build all
parts, so all dependencies should be present; otherwise some components will not build.
The following packages are used:

 * CLN (Class Library for Numbers); needed by hex-db, hex-dwba.
 * GiNaC; needed by hex-db.
 * SQLite3; needed by hex-db.
 * UMFPACK (part of SuiteSparse linear algebrapackage); needed by hex-ecs.
 * GSL (GNU Scientific library); needed by hex-db, hex-dwba and hex-ecs.
 * FFTW3; needed by hex-db, hex-dwba and hex-ecs.
 * OpenBLAS (or any other BLAS + LAPACK implementation that defines ZSBMV); needed by hex-ecs.
 * MPI (optional); used by hex-ecs.
 * PNG (optional); used by hex-ecs for debugging purposes.

CMake will try to locate the packages on its own, but will very likely fail, if they
are in non-standard locations. In that case you need to define the following paths:

    CLN_INCLUDE_DIRS        CLN_LIBRARIES
    GINAC_INCLUDE_DIRS      GINAC_LIBRARIES
    SQLITE3_INCLUDE_DIRS    SQLITE3_LIBRARIES
    UMFPACK_INCLUDE_DIRS    UMFPACK_LIBRARIES
    GSL_INCLUDE_DIRS        GSL_LIBRARIES
    FFTW3_INCLUDE_DIRS      FFTW3_LIBRARIES
                            OPENBLAS_LIBRARIES
    PNG_INCLUDE_DIRS        PNG_LIBRARIES

Compilation example:

    mkdir build
    cd build
    cmake -DUMFPACK_INCLUDE_DIRS=/usr/local/include/suitesparse -DUMFPACK_LIBRARIES=/usr/local/lib64/libumfpack.so ..
    make
