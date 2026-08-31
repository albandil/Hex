//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2018, Jakub Benda, Charles University in Prague                    //
//                                                                                   //
// MIT License:                                                                      //
//                                                                                   //
//  Permission is hereby granted, free of charge, to any person obtaining a          //
// copy of this software and associated documentation files (the "Software"),        //
// to deal in the Software without restriction, including without limitation         //
// the rights to use, copy, modify, merge, publish, distribute, sublicense,          //
// and/or sell copies of the Software, and to permit persons to whom the             //
// Software is furnished to do so, subject to the following conditions:              //
//                                                                                   //
//  The above copyright notice and this permission notice shall be included          //
// in all copies or substantial portions of the Software.                            //
//                                                                                   //
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS          //
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       //
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE       //
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, //
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF         //
// OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.  //
//                                                                                   //
//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //

#ifndef HEX_IO_H
#define HEX_IO_H

// --------------------------------------------------------------------------------- //

#include <cstdio>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>

// --------------------------------------------------------------------------------- //

#ifdef _OPENMP
    #include <omp.h>
#endif

// --------------------------------------------------------------------------------- //

#include "hex-arrays.h"

// --------------------------------------------------------------------------------- //

#include "bspline.h"
#include "luft.h"
#include "parallel.h"

// --------------------------------------------------------------------------------- //

/**
 * @brief Command line parameters.
 * 
 * This class uses @ref ParseCommandLine for parsing of the command line.
 * The information on valid switches can be retrieved by executing
 * <pre>
 * hex-ecs --help
 * </pre>
 * Some more general comments on the switches are in the above mentioned
 * class.
 */
class CommandLine
{
    public:

        /**
         * @brief Stages of the computations.
         *
         * Different stages of the computations, used to reflect user's choice
         * from the command line. The corresponding mapping is:
         * <table>
         * <tr><th>Command line option</th><th>Program itinerary</th></tr>
         * <tr><td>(none)</td><td>StgRadial | StgSolve | StgExtract</td></tr>
         * <tr><td>\--stg-integ</td><td>StgRadial</td></tr>
         * <tr><td>\--stg-integ-solve</td><td>StgRadial | StgSolve</td></tr>
         * <tr><td>\--stg-extract</td><td>StgiExtract</td></tr>
         * </table>
         */
        typedef enum {
            /// Do not start any intensive computations.
            StgNone    = 0x00,
            /// Compute only the radial integrals necessary for the construction of the equations.
            StgRadial  = 0x01,
            /// Solve the equations.
            StgSolve   = 0x02,
            /// Extract the T-matrices.
            StgExtract = 0x04
        } HexEcsStg;

        // constructor
        CommandLine (int argc, char* argv[])
            : writegrid(false), zipdata(), parallel(false), preconditioner("ILU"),
              droptol(1e-8), itinerary(StgNone), outofcore(false), cont(false), wholematrix(false), cache_all_radint(true), cache_own_radint(true),
              itertol(1e-8), prec_itertol(1e-8), parallel_precondition(false), gpu_large_data(false), lightweight_simple(false),
              lightweight_full(false), lightweight_radial_cache(false), shared_scratch(false), reuse_dia_blocks(false),
              kpa_simple_rad(false), ocl_platform(0), ocl_device(-1), factorizer("umfpack"), groupsize(1),
              parallel_factorization(false), parallel_extraction(true), ilu_max_iter(10), max_sub_iter(0), fail_on_sub_iter(true),
              carry_initial_guess(false), gpu_multiply(false), extract_extrapolate(false), extract_rho(-1), extract_rho_begin(-1), extract_samples(-1),
              refine_solution(false), map_solution(), map_solution_target(), ssor(-1), noluupdate(false), coupling_limit(-1), couple_all(true),
              mumps_outofcore(false), mumps_verbose(0), mumps_relax(20), kpa_drop(-1), write_intermediate_solutions(false),
              fast_bessel(false), hyb_additional_levels(0), multigrid_depth(0), multigrid_coarse_prec(0), dom_x_panels(1), dom_y_panels(1),
              dom_preconditioner("ILU"), dom_sweeps(-1), scratch(std::getenv("SCRATCHDIR") ? std::getenv("SCRATCHDIR") : "."), analytic_eigenstates(false),
              runtime_postprocess(false), sub_prec_verbose(false), multi_rhs(false), fpe(false), mumps_virtual_memory(false), nthreads(1),
              checkpoints(false), autostop_tolerance(0), purge(-1), fold_exchange(false)
        {
#ifdef _OPENMP
            // initialize the number of threads to OMP_NUM_THREADS
            #pragma omp parallel
            #pragma omp master
            nthreads = omp_get_num_threads();
#endif

            // get command line options
            parse(argc, argv);

            // run the whole sequence if nothing specified
            if (itinerary == StgNone)
                itinerary = StgRadial | StgSolve | StgExtract;
        }

        /// Read options from command line.
        void parse (int argc, char* argv[]);

        //
        // public attributes
        //

        /// Alternative name for the input file. Default is "hex.inp".
        std::ifstream inputfile;

        // Write grid layout to a VTK file.
        bool writegrid;

        /// Calculate bound state wave function of two-electron system.
        bool bound{};

        /// Evaluate dipole transition amplitudes from the bound state to all scattering states.
        bool dipoles{};

        /// A B-spline expansion of a solution to "zip". See \ref Bspline::zip .
        struct s_zipdata
        {
            // default constructor
            s_zipdata () : file(), Xmin(-1), Ymin(-1), Xmax(-1), Ymax(-1), nX(-1), nY(-1) {}

            std::string file;
            Real Xmin, Ymin, Xmax, Ymax;
            int nX, nY;
        }
        zipdata;

        /// Whether to use MPI.
        bool parallel;

        /// Preconditioner to use.
        std::string preconditioner;

        /// Drop tolerance for the ILU preconditioner.
        Real droptol;

        /// Which parts of the computation to run.
        int itinerary;

        /// Whether to keep precomputed data only on disk and spare RAM.
        bool outofcore;

        /// Whether to continue out-of-core computation from last computed solution (needs corresponding ooc files).
        bool cont;

        /// Whether to load full matrix from scratch disk at once when calculating dot product (and not by parts).
        bool wholematrix;

        /// Whether to keep radial integrals in memory.
        bool cache_all_radint;
        bool cache_own_radint;

        /// Tolerance for terminating iterative solution.
        Real itertol;

        /// Tolerance for terminating block preconditioner.
        Real prec_itertol;

        /// Whether to use OpenMP parallelization to run preconditioner for several blocks simultaneously.
        bool parallel_precondition;

        /// Keep large data in RAM instead of copying to the OpenCL compute device.
        bool gpu_large_data;

        /// Whether to avoid explicitly calculating big matrices and only apply them on the fly.
        bool lightweight_simple, lightweight_full, lightweight_radial_cache;

        /// Whether to compute only a subset of radial integrals in shared scratch architecture.
        bool shared_scratch;

        /// Whether to use diagonal blocks as present in the scratch directory. (For debugging purposes only.)
        bool reuse_dia_blocks;

        /// Use simplified radial integral matrix for nested KPA iterations (experimental).
        bool kpa_simple_rad;

        /// Index of OpenCL platform to use.
        unsigned ocl_platform;

        /// Index of OpenCL device to use.
        int ocl_device;

        /// LU-factorizer.
        std::string factorizer;

        /// Size of the local MPI communicator, used for distributed SuperLU.
        int groupsize;

        /// Allow parallel factorization.
        int parallel_factorization;

        /// Allow parallel extraction.
        int parallel_extraction;

        /// Maximal number of ILU iterations for hybrid preconditioner.
        int ilu_max_iter;

        /// Maximal number of sub-preconditioner iterations.
        int max_sub_iter;

        /// Stop with error message when sub-preconditioner fails to converge.
        bool fail_on_sub_iter;

        /// Whether to use previous-energy solution as an initial guess.
        bool carry_initial_guess;

        /// Do the sparse matrix multiplication on the OpenCL device (memory intensive!).
        bool gpu_multiply;

        /// Whether to radially extrapolate the extracted T-matrix.
        bool extract_extrapolate;

        /// Extraction radius.
        Real extract_rho;

        /// Radial distance where to start radial averaging/extrapolation of the T-matrix.
        Real extract_rho_begin;

        /// Extraction averaging/extrapolation sample count.
        int extract_samples;

        /// Use angular momentum eigenstates for extraction of T-matrices.
        bool eigenchannels{};

        /// Load also existing solutions and check that they are within tolerance.
        bool refine_solution;

        /// Map solution between different B-spline bases.
        std::vector<std::string> map_solution;

        /// Target mapping basis.
        std::string map_solution_target;

        /// Apply SSOR coupling.
        Real ssor;

        /// Keep calculated LU also for next energy.
        bool noluupdate;

        /// Maximal multipole considered by the coupled preconditioner.
        int coupling_limit;

        /// Couple all blocks when using coupled preconditioner, or just open channels.
        bool couple_all;

        /// MUMPS out of core
        bool mumps_outofcore;

        /// MUMPS diagnostic information
        int mumps_verbose;

        /// MUMPS memory relaxation factor (%)
        double mumps_relax;

        /// Whether to use drop tolerance for KPA preconditioner.
        Real kpa_drop;

        /// Write intermediate solutions.
        bool write_intermediate_solutions;

        /// Use faster Bessel function evaluation routine (not the Steed/Barnett) when calculating RHS.
        bool fast_bessel;

        /// Additional levels to be solved by ILU preconditioner when using HYB preconditioner.
        int hyb_additional_levels;

        /// Depth of the geometric multigrid.
        int multigrid_depth;

        /// What preconditioner to use for solution of the coarse problem.
        int multigrid_coarse_prec;

        /// Domain decomposition panels.
        int dom_x_panels, dom_y_panels;

        /// Domain preconditioner.
        std::string dom_preconditioner;

        /// Domain decomposition sweeps.
        int dom_sweeps;

        /// Scratch directory for out-of-core data.
        std::string scratch;

        /// Use analytic eigenstates instead of those obtained by diagonalization.
        bool analytic_eigenstates;

        /// Calculate T-matrices after every iteration.
        bool runtime_postprocess;

        /// Verbosity of the sub-preconditioner.
        bool sub_prec_verbose;

        /// Solve for multiple initial states at once.
        bool multi_rhs;

        /// Raise SIGFPE on invalid numerical result.
        bool fpe;

        /// Use virtual memory for MUMPS factors.
        bool mumps_virtual_memory;

        /// Nnmber of OpenMP threads.
        int nthreads;

        /// Write checkpoints for recovery of the solver.
        bool checkpoints;

        /// Monitor K-matrix convergence and stop the linear solver when reached.
        double autostop_tolerance;

        /// Delete old run-time post-processing directories.
        int purge;

        /**
         * @brief Solve only the angular blocks with l1 <= l2.
         *
         * The two angular blocks (l1,l2) and (l2,l1) of an electron-impact solution are
         * related by the exchange symmetry, so only one of them carries information. When
         * this switch is set, the blocks with l1 > l2 are dropped from the solved basis
         * and reconstructed from their counterparts whenever they are needed.
         */
        bool fold_exchange;

        /// Input file for solutions to use at right-hand side with dipole potential.
        std::vector<std::string> rhs_dipV;

        /// Switch to length gauge in Dalgarno-Lewis mode.
        bool length_gauge{false};
};

/**
 * @brief Input parameters.
 * 
 * This class contains information from the "hex.inp" file -- the initial
 * state data and the total conserved data.
 * 
 * Attributes are:
 * 
 * - inputfile: Handle to the input file.
 * - order: Order of the B-spline basis (typically 4).
 * - ecstheta: Complex grid rotation.
 * - rknots: Real knots (including R0).
 * - cknots: Unrotated complex knots (including R0).
 * - ni: Initial atomic energy state.
 * - instates: Initial atomic states. All of them will have the same energy.
 * - outstates: Final atomic states. Magnetic numbers will not be set as
 *                   all of them will be computed.
 * - L: Total angular momentum (partial wave).
 * - Spin: Total spin (partial wave).
 * - Pi: Total parity (partial wave).
 * - levels: How many angular pairs to use.
 * - maxell: Per-electron angular momentum restriction (determines the block 
 *               size of the matrix).
 * - Ei: Initial projectile energies [Ry].
 * - ki: Initial projectile momenta.
 * - B: Homogeneous magnetic field in the direction of z-axis.
 */
class InputFile
{
    public:

        // constructor
        InputFile (std::ifstream & inputfile)
        {
            // read inputfile
            read(inputfile);

            // maximal angular momentum (the same for both parities)
            maxell = levels + L;
        }

        // read data from file
        void read (std::ifstream & inputfile);

        //
        // public attributes
        //

            // B-spline order
            int order;

            // initial atomic principal quantum number
            int ni;

            // total angular momentum
            int L;

            // total parity
            int Pi;

            // 'nL', the limit on number of coupled angular states;
            // there will be 'nL * (L + 1 - Pi)' coupled angular state pairs
            int levels;

            // upper limit on the smaller of the one-electron angular momenta
            int limit;

            // maximal one-electron orbital quantum number
            int maxell;

            // total spins to calculate
            iArray Spin;

            // ECS rotation angle
            Real ecstheta;

            // weak magnetic field in atomic units (involved only perturbatively)
            Real B;

            // maximal energy (Ry) of states included in the asymptotic (outer) region
            Real channel_max_E;

            // real B-spline knots
            rArray rknots;

            // real B-spline knot projectile extention
            rArray rknots_ext;

            // complex-to-become knots (after rotation)
            rArray cknots;

            // total energies for which to solve the system
            rArray Etot;

            // maximal total energy
            Real max_Etot;

            // maximal energy of a target bound state
            Real max_Ebound;

            // initial and final atomic states
            std::vector<std::tuple<int,int,int>> instates, outstates;

            // atom charge (must be positive integer)
            Real Za;

            // projectile charge (only sign)
            Real Zp;

            // whether to calculate just the inner problem (decided from the knot sequence)
            bool inner_only;

            // keep only l1 <= l2; this is useful for large total angular momenta for reduction of the angular basis
            bool exchange;
};

/**
 * @brief Sign that relates a solution block to its mirror image.
 *
 * The (anti)symmetry of the wave function with respect to the exchange of the two
 * electrons relates the angular block (l2,l1) to the block (l1,l2) by the factor
 * @f$ (-1)^{l_1+l_2+L+S} @f$. Every state of the basis satisfies
 * @f$ l_1 + l_2 = 2n + L + \Pi @f$, so the same sign @f$ (-1)^{S+\Pi} @f$ applies to
 * all blocks of a given calculation. It is also the sign that anti/symmetrizes the
 * right-hand side.
 */
inline Real exchange_sign (int S, int Pi)
{
    return (S + Pi) % 2 == 0 ? 1.0_r : -1.0_r;
}

/**
 * @brief Mirror a solution segment.
 *
 * Rewrites the solution segment @c src of some angular block (l1,l2) into the segment
 * @c dst of the mirror block (l2,l1), which amounts to interchanging the two electrons.
 * The inner region is a square matrix of B-spline coefficients and is transposed; the
 * two groups of asymptotic channels describe the escape of either of the electrons, so
 * they change places. The sign that also relates the two blocks is @ref exchange_sign
 * and is *not* applied here.
 *
 * Both segments have the layout
 * <pre>
 *   [ Nspline_inner x Nspline_inner ][ Nchan1 x Nspline_outer ][ Nchan2 x Nspline_outer ]
 * </pre>
 * except that the mirror block has the channel counts swapped. The two electrons must
 * share the same B-spline basis, otherwise the two layouts are not compatible at all.
 */
void mirror_segment
(
    const cArrayView src,
    cArrayView dst,
    int Nspline_inner,
    int Nspline_outer,
    int Nchan1,
    int Nchan2
);

/**
 * @brief Solution input/output class.
 * 
 * On several places in code it is necessary to load some solutions. In order not
 * to bother with the file names and still to keep them consistent, this class
 * will carry the file name and do all operations on it.
 *
 * When the exchange symmetry has been folded into the equations (the command line
 * option \--fold-exchange), the angular blocks with l1 > l2 are never solved for and
 * their segments are not written at all. Such a block is restored on reading from the
 * segment of its mirror image, so that the rest of the program still sees a solution
 * that covers the complete angular basis. What makes this possible is the size of the
 * inner region, which the writer of the mirror block has to record in its file; see
 * the constructor argument @c Nspline_inner.
 */
class SolutionIO
{
    public:

        SolutionIO () : L_(0), S_(0), Pi_(0), ni_(0), li_(0), mi_(0), E_(0), Nspline_inner_(0) {}

        /**
         * @brief Constructor.
         *
         * @param Nspline_inner Size of the square inner region of a solution segment.
         *        It is written into every segment file, where it is what a reader needs
         *        to reconstruct the mirror block from the block itself. Pass it only
         *        when the mirror blocks are really left out of the output, so that a
         *        reader never mirrors a solution that was not meant to be mirrored;
         *        zero, the default, writes all blocks and asks for no reconstruction.
         */
        SolutionIO
        (
            int L, int S, int Pi,
            int ni, int li, int mi,
            Real E,
            std::vector<std::pair<int,int>> const & ang,
            std::vector<std::pair<int,int>> const & chann = std::vector<std::pair<int,int>>(),
            std::string prefix = "psi",
            int Nspline_inner = 0
        ) : L_(L), S_(S), Pi_(Pi),
            ni_(ni), li_(li), mi_(mi),
            E_(E), ang_(ang), prefix_(prefix), Nspline_inner_(Nspline_inner)
        {
            if (chann.empty())
            {
                // initialize channel counts to zero
                chann_.resize(ang.size(), std::make_pair(0,0));
            }
            else
            {
                // check consistency of the parameters
                if (ang.size() != chann.size())
                    HexException("Incompatible size of angular momentum pairs list and channel count list.");

                // use the supplied channel counts
                chann_ = chann;
            }
        }

        /// All blocks flag.
        static const int All = -1;

        /// Get name of the solution file.
        std::string name (int ill = SolutionIO::All) const
        {
            if (ill == -1)
                return format("%s-%g-%d-%d-%d-%d-%d-%d.hdf", prefix_.c_str(), E_ + 1, L_, S_, Pi_, ni_, li_, mi_);
            else
                return format("%s-%g-%d-%d-%d-%d-%d-%d-(%d,%d).hdf", prefix_.c_str(), E_ + 1, L_, S_, Pi_, ni_, li_, mi_, ang_[ill].first, ang_[ill].second);
        }

        /**
         * @brief Number of elements of a solution segment, or zero when there is none.
         *
         * A block whose own file is missing counts with the size of the file of its
         * mirror image, which it would be reconstructed from and which has the same size.
         */
        std::size_t segment_size (int ill) const
        {
            HDFFile hdf (name(ill), HDFFile::readonly);

            if (hdf.valid())
                return hdf.size("array") / 2;

            int Nspline_inner = 0;
            int illm = mirror_source(ill, Nspline_inner);

            if (illm < 0)
                return 0;

            return HDFFile(name(illm), HDFFile::readonly).size("array") / 2;
        }

        /// Position of the block (l2,l1) in the list of angular states, or -1.
        int mirror (int ill) const
        {
            std::pair<int,int> other (ang_[ill].second, ang_[ill].first);
            auto position = std::find(ang_.begin(), ang_.end(), other);
            return position == ang_.end() ? -1 : position - ang_.begin();
        }

        /**
         * @brief Segment file that a missing block can be reconstructed from.
         *
         * Returns the position of the mirror block (l2,l1) when the segment file of the
         * block @c ill is itself missing, the file of the mirror block is there, and that
         * file records the size of the inner region -- which only a run that deliberately
         * leaves the mirror blocks out ever writes. Returns -1 in every other case, so
         * that an incomplete output of an ordinary run is not mistaken for a folded one.
         * On success, @c Nspline_inner is set from the file.
         */
        int mirror_source (int ill, int & Nspline_inner) const
        {
            if (HDFFile(name(ill), HDFFile::readonly).valid())
                return -1;

            int illm = mirror(ill);

            if (illm < 0 or illm == ill)
                return -1;

            HDFFile hdf (name(illm), HDFFile::readonly);

            if (not hdf.valid() or not hdf.read("Nspline_inner", &Nspline_inner, 1) or Nspline_inner <= 0)
                return -1;

            return illm;
        }

        /// Check that the file exists.
        bool check (int ill, std::size_t & total_size) const
        {
            // look for specific solution segment file
            if (ill != SolutionIO::All)
                return segment_size(ill) > 0;

            // look for all solution segment files
            std::vector<std::size_t> size (ang_.size());
            for (unsigned illp = 0; illp < ang_.size(); illp++)
                size[illp] = segment_size(illp);

            // calculate total size
            total_size = std::accumulate(size.begin(), size.end(), 0);

            // check that all blocks existed
            return std::find(size.begin(), size.end(), 0) == size.end();
        }
        bool check (int ill = SolutionIO::All) const
        {
            std::size_t total_size = 0;
            return check(ill, total_size);
        }

        /**
         * @brief Load the solution from disk.
         * 
         * This function will try to load the solution HDF file (of name 'name_')
         * into a complex array passed as argument. If the file is not found, array
         * is not written and 'false' is returned. Otherwise the function returns 'true'
         * and fills the contents of the array with the read data.
         */
        bool load (BlockArray<Complex> & solution, int ill = SolutionIO::All)
        {
            // check that the requested files are present
            if (not check(ill))
                return false;

            // select segments
            iArray segments_to_load = { ill };
            if (ill == SolutionIO::All)
                segments_to_load = linspace<int>(0, solution.size() - 1, solution.size());

            // for all blocks to load
            for (int i : segments_to_load)
            {
                //
                // Either read from the whole monolithic solution file ...
                //

                    HDFFile hdf (name(), HDFFile::readonly);
                    if (hdf.valid())
                    {
                        // get segment size
                        std::size_t segsize = hdf.size("array") / ang_.size();

                        // allocate memory
                        solution[i].resize(segsize);

                        // read data from the file
                        if (not hdf.read("array", solution[i].data(), segsize, i * segsize))
                            return false;

                        // dump to disk
                        if (not solution.inmemory())
                        {
                            solution.hdfsave(i);
                            solution[i].drop();
                        }
                    }

                //
                // ... or read from a solution segment.
                //

                    else
                    {
                        // load the requested solution segment file, or restore the segment
                        // from the mirror block when it has not been written at all
                        if (not load_segment(solution[i], i))
                            return false;

                        // dump to disk
                        if (not solution.inmemory())
                        {
                            solution.hdfsave(i);
                            solution[i].drop();
                        }
                    }
            }

            return true;
        }

        /**
         * @brief Write solution to disk.
         * 
         * Save solution to disk as a HDF file. Zero blocks
         * of the solution are automatically detected and compressed
         * (substituted by position & length information). The function
         * return 'true' when write was successful, 'false' otherwise.
         */
        bool save (BlockArray<Complex> const & solution, int ill = SolutionIO::All) const
        {
            // select segments
            iArray segments_to_save = { ill };
            if (ill == SolutionIO::All)
                segments_to_save = linspace<int>(0, solution.size() - 1, solution.size());

            // for all segments
            for (int iseg : segments_to_save)
            {
                // load data
                if (not solution.inmemory())
                    const_cast<BlockArray<Complex>&>(solution).hdfload(iseg);

                // write data
                if (not save_segment(solution[iseg], iseg))
                    return false;

                // unload data
                if (not solution.inmemory())
                    const_cast<BlockArray<Complex>&>(solution)[iseg].drop();
            }

            // success
            return true;
        }

        /**
         * @brief Read one solution segment.
         *
         * Reads the segment of the angular block @c ill and the numbers of its asymptotic
         * channels. When the block has not been written, because the exchange symmetry
         * was folded into the equations, its segment is restored from the segment of the
         * mirror block (l2,l1) by interchanging the two electrons in it.
         */
        bool load_segment (cArray & segment, int ill)
        {
            HDFFile hdf (name(ill), HDFFile::readonly);

            if (hdf.valid())
            {
                if (not segment.hdfload(name(ill)))
                    return false;

                hdf.read("Nchan1", &chann_[ill].first,  1);
                hdf.read("Nchan2", &chann_[ill].second, 1);

                return true;
            }

            // this block was not solved for; take the mirror block instead
            int Nspline_inner = 0;
            int illm = mirror_source(ill, Nspline_inner);

            if (illm < 0)
                return false;

            cArray source;
            if (not source.hdfload(name(illm)))
                return false;

            // the mirror block has the two electrons, and so its channel groups, interchanged
            int Nchan1 = 0, Nchan2 = 0;
            HDFFile hdfm (name(illm), HDFFile::readonly);
            hdfm.read("Nchan1", &Nchan1, 1);
            hdfm.read("Nchan2", &Nchan2, 1);

            std::size_t inner = std::size_t(Nspline_inner) * std::size_t(Nspline_inner);

            if (source.size() < inner or (source.size() > inner and Nchan1 + Nchan2 == 0))
                HexException("The solution file \"%s\" does not have the size its own metadata imply.", name(illm).c_str());

            int Nspline_outer = (Nchan1 + Nchan2 == 0 ? 0 : (source.size() - inner) / (Nchan1 + Nchan2));

            segment.resize(source.size());
            mirror_segment(source, segment, Nspline_inner, Nspline_outer, Nchan1, Nchan2);
            segment *= exchange_sign(S_, Pi_);

            chann_[ill] = std::make_pair(Nchan2, Nchan1);

            return true;
        }

        bool save_segment (const cArrayView solution, int ill) const
        {
            // open file
            HDFFile hdf (name(ill), HDFFile::overwrite);

            // check file status
            if (not hdf.valid())
                return false;

            // write data
            bool success = true;
            success = (success and hdf.write("L",  &L_,  1));
            success = (success and hdf.write("S",  &S_,  1));
            success = (success and hdf.write("Pi", &Pi_, 1));
            success = (success and hdf.write("ni", &ni_, 1));
            success = (success and hdf.write("li", &li_, 1));
            success = (success and hdf.write("mi", &mi_, 1));
            success = (success and hdf.write("E",  &E_,  1));
            success = (success and hdf.write("l1",  &ang_[ill].first,  1));
            success = (success and hdf.write("l2",  &ang_[ill].second, 1));
            success = (success and hdf.write("Nchan1", &chann_[ill].first, 1));
            success = (success and hdf.write("Nchan2", &chann_[ill].second, 1));
            success = (success and hdf.write("array", solution.data(), solution.size()));

            // this is what lets a reader restore the mirror block from this one
            if (Nspline_inner_ > 0)
                success = (success and hdf.write("Nspline_inner", &Nspline_inner_, 1));

            return success;
        }

        std::vector<std::pair<int,int>> const & channels () const { return chann_; }

    private:

        int L_, S_, Pi_, ni_, li_, mi_;
        Real E_;
        std::vector<std::pair<int,int>> ang_, chann_;
        std::string prefix_;

        // Size of the inner region to record in the written segments, or zero when the
        // segments are not meant to be mirrored. See the constructor.
        int Nspline_inner_;
};

/**
 * @brief Convert solution to VTK file.
 * 
 * This function will convert a solution (the filename of the HDF file
 * should come from the command line option --zipfile) to a series of VTK text
 * files that can be viewed e.g. by the free ParaView program. The resulting
 * VTK files (one for each angular basis vector) will contain real and imaginary
 * part of the wave function component on an rectilinear grid. Parameters of the
 * grid can be further refined by the command line arguments --zipcount and --zipmax).
 * 
 * @param cmd Class containing command line options.
 * @param bspline B-spline environments. The appropriate one will be chosen
 *                by the match between the B-spline count and the solution size.
 * @param ll Angular basis.
 * 
 * See the respective classes for deeper explanation of individual parameters.
 */
void zip_solution
(
    CommandLine const & cmd,
    InputFile const & inp,
    Parallel const & par, 
    Bspline const & bspline_inner,
    Bspline const & bspline_full,
    std::vector<std::pair<int,int>> const & ll
);

/**
 * @brief Write grid to a VTK file.
 */
void write_grid
(
    Bspline const & bspline,
    std::string const & basename
);

#endif /* HEX_IO_H */
