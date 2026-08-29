HEX-ECS
=======

`hex-ecs` solves the electron-hydrogen scattering problem by the method of exterior complex scaling (ECS). The two-electron
wave function is expanded in a B-spline basis and the resulting sparse linear system is solved iteratively by preconditioned
conjugate orthogonal conjugate gradients (PCOCG). The program writes the scattering T-matrices as SQL batch files, which are
then imported into a database by `hex-db`.


Requirements
------------

`hex-ecs` is built together with the rest of the suite, so it shares the requirements listed in the top-level `README.md`.
Of those, it needs a C++17 compiler and CMake 3.13 or later, BLAS and LAPACK, and the GNU Scientific Library. The code uses
`std::filesystem`, so with GCC this means version 9 or later - the build does not add the separate `-lstdc++fs` that GCC 8
still requires.

The optional libraries relevant to this program are MPI, the LU factorization packages listed below, OpenCL for the GPU
preconditioner, and libpng for the debugging matrix plots.


Building
--------

The program is built from the top-level directory of the suite, not from here:

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build
    cmake --install build --prefix /usr/local

Optional libraries are enabled by the `WITH_*` switches and located by the matching `*_INCLUDE_DIRS` and `*_LIBRARIES`
variables, for example

    cmake -S . -B build -DWITH_UMFPACK=ON -DUMFPACK_LIBRARIES=umfpack

See the top-level `CMakeLists.txt` for the complete list of switches.


Preconditioners and factorizers
-------------------------------

The iterative solver is preconditioned by one of the schemes selected with `--preconditioner`: `none`, `diag`, `cg`, `ILU`,
`KPA`, `HYB`, `DOM` and `coupled`, plus `GPU` when the program has been built with OpenCL. Run

    hex-ecs --list-preconditioners

for a short description of each and for the set that the binary actually contains. `ILU` is the default.

The preconditioners that decompose matrices (`ILU`, `coupled`, and `DOM` when its panels use `ILU`) delegate the work to one
of the factorization libraries selected with `--lu`: `lapack`, `umfpack`, `pardiso`, `superlu`, `superlu_dist`, `mumps` or
`scalapack`. All but `lapack` require the corresponding `WITH_*` switch at build time; `hex-ecs --help` prints those that
are available. The default is `umfpack`, so a build without any sparse factorizer is limited to the preconditioners that
need none, such as `KPA`.


Running
-------

A minimal run needs nothing but an input file, a sample of which the program writes itself:

    mkdir run && cd run
    hex-ecs --example             # writes "example.inp"
    hex-ecs --input example.inp   # solves it

The program has a large number of tuning options; `hex-ecs --help` lists them all.

A larger demonstration that also exercises `hex-db` is in the directory `test-run`. It computes the T-matrices for the four
lowest partial waves, imports them into a database and extracts the elastic differential cross sections, which it plots in
Gnuplot if available. Both programs have to be in `PATH`:

    cd test-run
    ./test-run-ILU.sh
    ./test-run-KPA.sh
    ./test-run-GPU.sh     # needs a build with OpenCL support

Each script works in a subdirectory of its own, where it leaves the cross sections as `singlet.dcs` and `triplet.dcs`, and
writes a log of the whole run beside that subdirectory. Reference copies of both are stored in `test-run-output`.


Tests
-----

The directory `tests` holds the regression tests of `hex-ecs`, driven by CTest. Each of them solves the sample problem with
a different combination of preconditioner, factorizer and parallelization mode, and checks that the integral cross sections
come out. They are complete calculations rather than unit tests, so expect minutes to hours per test. Run them from the
build directory with

    ctest --output-on-failure

Tests that would need a library the build was configured without are not created at all. The distributed ones use the MPI
launcher found during configuration; extra launcher options belong to the cache variable `MPIEXEC_PREFLAGS`, and with Open
MPI that should be at least `--bind-to none` so that the multi-threaded ranks are not confined to a single core. See
`tests/CMakeLists.txt` for the details.


Documentation
-------------

The reference documentation can be extracted from the sources by Doxygen. Running

    doxygen

in this directory processes `src` and `../common` and writes HTML into `html`.
