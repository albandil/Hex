# Building Hex

It is recommended to use CMake for building the Hex source tree. This will build all
parts, for which the dependencies are present. The following packages are used:

 * CLN (Class Library for Numbers); needed by hex-db, hex-dwba.
 * GiNaC; needed by hex-db.
 * SQLite3; needed by hex-db.
 * UMFPACK (part of SuiteSparse linear algebrapackage); needed by hex-ecs.
 * GSL (GNU Scientific library); needed by hex-db, hex-dwba and hex-ecs.
 * FFTW3; needed by hex-db, hex-dwba and hex-ecs.
 * OpenBLAS (or any other BLAS + LAPACK implementation that defines ZSBMV); needed by hex-ecs.
 * MPI (optional); used by hex-ecs.
 * PNG (optional); used by hex-ecs for debugging purposes.
 * OpenCL (optional); used by hex-ecs.

CMake will try to locate the packages on its own, but will very likely fail, if they
are in non-standard locations. In that case you need to define some the following paths:

    CLN_INCLUDE_DIRS        CLN_LIBRARIES
    GINAC_INCLUDE_DIRS      GINAC_LIBRARIES
    SQLITE3_INCLUDE_DIRS    SQLITE3_LIBRARIES
    UMFPACK_INCLUDE_DIRS    UMFPACK_LIBRARIES  (and similarly for all SuiteSparse components)
    GSL_INCLUDE_DIRS        GSL_LIBRARIES
    FFTW3_INCLUDE_DIRS      FFTW3_LIBRARIES
    OPENBLAS_INCLUDE_DIRS   OPENBLAS_LIBRARIES
    BLAS_INCLUDE_DIRS       BLAS_LIBRARIES
    LAPACK_INCLUDE_DIRS     LAPACK_LIBRARIES
    PNG_INCLUDE_DIRS        PNG_LIBRARIES
    OPENCL_INCLUDE_DIRS     OPENCL_LIBRARIES

Also, it can easily happen that your libraries are named in a different way than is the standard.
For example, you may want to name your DLLs openblas-0.2.14.dll instead of simple openblas.dll.
In order for CMake to recognize such library, it is necessary to add the name using the following
definitions.

    CLN_NAMES
    GINAC_NAMES
    SQLITE3_NAMES
    UMFPACK_NAMES (CHOLMOD_NAMES, METIS_NAMES, AMD_NAMES, CAMD_NAMES, ...)
    GSL_NAMES
    FFTW3_NAMES
    OPENBLAS_NAMES
    OPENCL_NAMES

Compilation example on Windows:

    mkdir build
    cd build
    cmake -G"MinGW Makefiles" ^
          -DCMAKE_C_COMPILER=C:/mingw-w64/x86_64-sjlj-win32-mingw32/mingw64/bin/x86_64-w64-mingw32-gcc.exe ^
          -DCMAKE_CXX_COMPILER=C:/mingw-w64/x86_64-sjlj-win32-mingw32/mingw64/bin/x86_64-w64-mingw32-g++.exe ^
          -DCMAKE_MAKE_EXECUTABLE=C:/mingw-w64/x86_64-sjlj-win32-mingw32/mingw64/bin/mingw32-make.exe ^
          -DCMAKE_INCLUDE_PATH=C:/Users/Jacob/Documents/Hex/Win64-include ^
          -DCMAKE_LIBRARY_PATH=C:/Users/Jacob/Documents/Hex/Win64-lib ^
          -DCLN_NAMES=cln-1.3.4 -DGINAC_NAMES=ginac-1.6.3 -DSQLITE3_NAMES=sqlite3-3.8.3 ^
          -DUMFPACK_NAMES=umfpack-5.7.1 -DCHOLMOD_NAMES=cholmod-3.0.5 -DMETIS_NAMES=metis-4.0.3 ^
          -DAMD_NAMES=amd-2.4.1 -DCAMD_NAMES=camd-2.4.1 -DCOLAMD_NAMES=colamd-2.9.1 ^
          -DCCOLAMD_NAMES=ccolamd-2.9.1 -DGSL_NAMES=gsl-1.16 -DGSLCBLAS_NAMES=gslcblas-1.16 ^
          -DFFTW3_NAMES=fftw3-3.3.4 -DSUITESPARSE_CONFIG_NAMES=suitesparse_config-4.4.3 ^
          -DOPENBLAS_NAMES=openblas-0.2.14 ^
          ..
    mingw32-make
