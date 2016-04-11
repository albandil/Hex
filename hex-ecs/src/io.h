//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2016, Jakub Benda, Charles University in Prague                    //
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

#include <cstdio>
#include <fstream>
#include <vector>
#include <string>

#include "hex-arrays.h"
#include "hex-luft.h"

#include "bspline.h"

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
            : writegrid(false), zipdata(), parallel(false), preconditioner(0),
              droptol(1e-8), itinerary(StgNone), outofcore(false), cont(false), wholematrix(false), cache_all_radint(true), cache_own_radint(true),
              itertol(1e-8), prec_itertol(1e-8), parallel_precondition(false), gpu_large_data(false),
              lightweight_full(false), lightweight_radial_cache(false), shared_scratch(false), reuse_dia_blocks(false),
              kpa_simple_rad(false), ocl_platform(0), ocl_device(0), factorizer(LUFT_ANY), groupsize(1), panels(1),
              parallel_factorization(false), parallel_extraction(true), kpa_max_iter(-1), ilu_max_blocks(-1),
              carry_initial_guess(false), gpu_multiply(false), extract_extrapolate(false), extract_rho(-1), extract_rho_begin(-1), extract_samples(-1),
              refine_solution(false), map_solution(), map_solution_target(), ssor(-1), noluupdate(false), coupling_limit(1000),
              gpu_host_multiply(false), mumps_outofcore(true), mumps_verbose(0), kpa_drop(-1), exact_rhs(false), polarization(0)
        {
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
        
        /// A B-spline expansion of a solution to "zip". See \ref Bspline::zip .
        struct s_zipdata
        {
            // default constructor
            s_zipdata () : file(), Xmin(-1), Ymin(-1), Xmax(-1), Ymax(-1), nX(-1), nY(-1) {}
            
            std::string file;
            double Xmin, Ymin, Xmax, Ymax;
            int nX, nY;
        }
        zipdata;
        
        /// Whether to use MPI.
        bool parallel;
        
        /// %Preconditioner to use. See \ref Preconditioners::AvailableTypes for available types.
        int preconditioner;
        
        /// Drop tolerance for the ILU preconditioner.
        double droptol;
        
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
        double itertol;
        
        /// Tolerance for terminating block preconditioner.
        double prec_itertol;
        
        /// Whether to use OpenMP parallelization to run preconditioner for several blocks simultaneously.
        bool parallel_precondition;
        
        /// Keep large data in RAM instead of copying to the OpenCL compute device.
        bool gpu_large_data;
        
        /// Whether to avoid explicitly calculating big matrices and only apply them on the fly.
        bool lightweight_full, lightweight_radial_cache;
        
        /// Whether to compute only a subset of radial integrals in shared scratch architecture.
        bool shared_scratch;
        
        /// Whether to use diagonal blocks as present in the scratch directory. (For debugging purposes only.)
        bool reuse_dia_blocks;
        
        /// Use simplified radial integral matrix for nested KPA iterations (experimental).
        bool kpa_simple_rad;
        
        /// Index of OpenCL platform to use.
        unsigned ocl_platform;
        
        /// Index of OpenCL device to use.
        unsigned ocl_device;
        
        /// LU-factorizer.
        int factorizer;
        
        /// Size of the local MPI communicator, used for distributed SuperLU.
        int groupsize;
        
        /// Number of panels (solver + propagator sections).
        int panels;
        
        /// Allow parallel factorization.
        int parallel_factorization;
        
        /// Allow parallel extraction.
        int parallel_extraction;
        
        /// Maximal number of KPA iterations for hybrid preconditioner.
        int kpa_max_iter;
        
        /// Maximal number of ILU-preconditioned blocks for hybrid preconditioner.
        int ilu_max_blocks;
        
        /// Whether to use previous-energy solution as an initial guess.
        bool carry_initial_guess;
        
        /// Do the sparse matrix multiplication on the OpenCL device (memory intensive!).
        bool gpu_multiply;
        
        /// Whether to radially extrapolate the extracted T-matrix.
        bool extract_extrapolate;
        
        /// Extraction radius.
        double extract_rho;
        
        /// Radial distance where to start radial averaging/extrapolation of the T-matrix.
        double extract_rho_begin;
        
        /// Extraction averaging/extrapolation sample count.
        int extract_samples;
        
        /// Load also existing solutions and check that they are within tolerance.
        bool refine_solution;
        
        /// Map solution between different B-spline bases.
        std::vector<std::string> map_solution;
        
        /// Target mapping basis.
        std::string map_solution_target;
        
        /// Apply SSOR coupling.
        double ssor;
        
        /// Keep calculated LU also for next energy.
        bool noluupdate;
        
        /// Maximal multipole considered by the coupled preconditioner.
        int coupling_limit;
        
        /// Keep vectors in host memory when doing multiplication on GPU.
        bool gpu_host_multiply;
        
        /// MUMPS out of core
        bool mumps_outofcore;
        
        /// MUMPS diagnostic information
        int mumps_verbose;
        
        /// Whether to use drop tolerance for KPA preconditioner.
        double kpa_drop;
        
        /// Exact rhs.
        bool exact_rhs;
        
        /// Take into account the polarization of the target orbitals.
        double polarization;
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
            
            // compute angular momentum limit
            maxell = levels + L + Pi;
        }
        
        // read data from file
        void read (std::ifstream & inputfile);
        
        //
        // public attributes
        //
        
        int order, ni, L, Pi, levels, maxell;
        iArray Spin;
        double ecstheta, B;
        rArray rknots, rknots_next, cknots, overlap_knots, Etot;
        std::vector<std::tuple<int,int,int>> instates, outstates;
};

/**
 * @brief Solution input/output class.
 * 
 * On several places in code it is necessary to load some solutions. In order not
 * to bother with the file names and still to keep them consistent, this class
 * will carry the file name and do all operations on it.
 */
class SolutionIO
{
    public:
        
        SolutionIO (int L, int S, int Pi, int ni, int li, int mi, double E, std::vector<std::pair<int,int>> const & ang)
            : L_(L), S_(S), Pi_(Pi), ni_(ni), li_(li), mi_(mi), E_(E), ang_(ang) {}
        
        /// Get name of the solution file.
        std::string name (int ill = -1) const
        {
            if (ill == -1)
                return format("psi-%g-%d-%d-%d-%d-%d-%d.hdf", E_ + 1, L_, S_, Pi_, ni_, li_, mi_);
            else
                return format("psi-%g-%d-%d-%d-%d-%d-%d-(%d,%d).hdf", E_ + 1, L_, S_, Pi_, ni_, li_, mi_, ang_[ill].first, ang_[ill].second);
        }
        
        /// Check that the file exists, return size.
        std::size_t check (int ill = -1) const
        {
            // look for monolithic solution file
            HDFFile fmono (name(), HDFFile::readonly);
            if (fmono.valid()) return fmono.size("array")/2 / ang_.size();
            
            // look for specific solution segment file
            if (ill >= 0)
            {
                HDFFile fsingle (name(ill), HDFFile::readonly);
                return fsingle.valid() ? fsingle.size("array")/2 : 0;
            }
            
            // look for all solution segment files
            std::vector<std::size_t> size (ang_.size());
            for (unsigned illp = 0; illp < ang_.size(); illp++)
            {
                HDFFile fsingle (name(illp), HDFFile::readonly);
                size[illp] = (fsingle.valid() ? fsingle.size("array")/2 : 0);
            }
            
            // check that the sizes are consistent
            std::size_t min_size = *std::min_element(size.begin(), size.end());
            std::size_t max_size = *std::max_element(size.begin(), size.end());
            return min_size == max_size ? min_size : 0;
        }
        
        /**
         * @brief Load the solution from disk.
         * 
         * This function will try to load the solution HDF file (of name 'name_')
         * into a complex array passed as argument. If the file is not found, array
         * is not written and 'false' is returned. Otherwise the function returns 'true'
         * and fills the contents of the array with the read data.
         */
        bool load (BlockArray<Complex> & solution, int ill = -1)
        {
            // check that the requested files are present
            if (not check(ill))
                return false;
            
            // select segments
            iArray segments_to_load = { ill };
            if (ill == -1)
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
                        // simply load the requested solution segment file
                        if (not solution[i].hdfload(name(i)))
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
        bool save (BlockArray<Complex> const & solution, int ill = -1) const
        {
            // select segments
            iArray segments_to_save = { ill };
            if (ill == -1)
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
        
        bool save_segment (const cArrayView solution, int ill) const
        {
            // open file
            HDFFile hdf (name(ill), HDFFile::overwrite);
            
            // check file status
            if (not hdf.valid())
                return false;
            
            // write data
            if (not hdf.write("array", solution.data(), solution.size()))
                return false;
            
            // success
            return true;
        }
    
    private:
        
        int L_, S_, Pi_, ni_, li_, mi_;
        double E_;
        std::vector<std::pair<int,int>> ang_;
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
    CommandLine & cmd,
    std::vector<Bspline> const & bspline,
    std::vector<std::pair<int,int>> const & ll
);

/**
 * @brief Write grid to a VTK file.
 */
void write_grid
(
    std::vector<Bspline> const & bspline,
    std::string const & basename
);

#endif /* HEX_IO_H */
