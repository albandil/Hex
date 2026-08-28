HEX: Hydrogen-electron collision solver
=======================================

HEX provides a suite of C++ programs that solve the electron-hydrogen scattering problem.

The main component is the program `hex-ecs`, which implements exterior complex scaling (ECS) method and iterative solution
for the scattered wave using preconditioned conjugate orthogonal conjugate gradients (PCOCG).

For extrapolation of results to high angular momenta ("Born top-up"), the programs `hex-dwba` and `hex-fullborn` can be
used. `hex-dwba` computes partial-wave T-matrices in the distorted wave Born approximation, or in the plane wave Born
approximation when run with `--nodistort`. `hex-fullborn` evaluates the total plane wave Born cross section (without
exchange), summed over all angular momenta.

`hex-ecs` and `hex-dwba` produce scattering T-matrices in form of SQL files. These are aggregated by `hex-db` into a single
SQLite3 database. The same program `hex-db` can be then used to extract observables from this database. It can also
interpolate missing energy points.

The hard requirements for HEX are:

 - C++17 compiler and CMake 3.13 or later
 - BLAS and LAPACK (e.g. OpenBLAS or Intel MKL)
 - GNU Scientific Library
 - SQLite3

Some functionality requires the following:

 - CLN and GiNaC (symbolic Born amplitudes: the plane wave mode of `hex-dwba` and the dipole Born subtraction in `hex-db`)
 - MPI (for distributed runs of `hex-ecs`)

The following are nice to have, as they add some additional capabilities:

 - UMFPACK (part of SuiteSparse - for sparse LU decomposition in `hex-ecs`)
 - Pardiso (stand-alone library or Intel MKL Pardiso - for sparse LU decomposition in `hex-ecs`)
 - SuperLU (for sparse LU decomposition in `hex-ecs`)
 - MUMPS (for large-scale distributed LU decomposition in `hex-ecs`)
 - SuperLU_DIST (for large-scale distributed LU decomposition in `hex-ecs`)
 - ScaLAPACK (for distributed dense LU decomposition in `hex-ecs`)
 - HDF5 (needed by the `hex-hdf2hdf` conversion utility)
 - OpenCL (for GPU acceleration in `hex-ecs`)
 - libpng (printing matrix structure for debugging/illustrative purposes)
 - Doxygen (for generated documentation)

See `CMakeLists.txt` for details on enabling the optional features.
