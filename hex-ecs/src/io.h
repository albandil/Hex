//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2015, Jakub Benda, Charles University in Prague                    //
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

#include "arrays.h"
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
            : zipcount(0), zipmax(-1), parallel(false), preconditioner(0),
              droptol(1e-8), itinerary(StgNone), outofcore(false), wholematrix(false), cache_all_radint(true), cache_own_radint(true),
              itertol(1e-8), prec_itertol(1e-8), parallel_dot(false), parallel_block(false),
              gpu_slater(false), lightweight_full(false), lightweight_radial_cache(false), shared_scratch(false), reuse_dia_blocks(false),
              kpa_simple_rad(false), gpucpu(false)
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
        
        /// A B-spline expansion of a solution to "zip". See \ref Bspline::zip .
        std::string zipfile;
        
        /// How many equidistant samples on each axis to use.
        int  zipcount;
        
        /// Radial cutoff for solution zipping useful if one is interested only in the near-origin behaviour.
        double zipmax;
        
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
        
        /// Whether to load full matrix from scratch disk at once when calculating dot product (and not by parts).
        bool wholematrix;
        
        /// Whether to keep radial integrals in memory.
        bool cache_all_radint;
        bool cache_own_radint;
        
        /// Tolerance for terminating iterative solution.
        double itertol;
        
        /// Tolerance for terminating block preconditioner.
        double prec_itertol;
        
        /// Whether to use OpenMP parallelization in SymDiaMatrix::dot
        bool parallel_dot;
        
        /// Whether to use OpenMP parallelization to run preconditioner for several blocks simultaneously.
        bool parallel_block;
        
        /// Whether to compute diagonal two-electron integrals using OpenCL.
        bool gpu_slater;
        
        /// Whether to avoid explicitly calculating big matrices and only apply them on the fly.
        bool lightweight_full, lightweight_radial_cache;
        
        /// Whether to compute only a subset of radial integrals in shared scratch architecture.
        bool shared_scratch;
        
        /// Whether to use diagonal blocks as present in the scratch directory. (For debugging purposes only.)
        bool reuse_dia_blocks;
        
        /// Use simplified radial integral matrix for nested KPA iterations (experimental).
        bool kpa_simple_rad;
        
        /// If 'true', then use CPU for GPUPreconditioner.
        bool gpucpu;
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
            
            // projectile momenta
            ki.resize(Ei.size());
            for (std::size_t i = 0; i < Ei.size(); i++)
                ki[i] = sqrt(Ei[i]);
        }
        
        // read data from file
        void read (std::ifstream & inputfile);
        
        //
        // public attributes
        //
        
        int order, ni, L, Pi, levels, maxell;
        double ecstheta, B;
        rArray rknots, cknots, Ei, ki;
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
        
        SolutionIO (int L, int S, int Pi, int ni, int li, int mi, double E, std::vector<std::pair<int,int>> const & ang, unsigned Nspline)
            : L_(L), S_(S), Pi_(Pi), ni_(ni), li_(li), mi_(mi), E_(E), ang_(ang), Nspline_(Nspline) {}
        
        /// Get name of the solution file.
        std::string name (int ill = -1) const
        {
            if (ill == -1)
                return format("psi-%d-%d-%d-%d-%d-%d-%g.hdf", L_, S_, Pi_, ni_, li_, mi_, E_);
            else
                return format("psi-%d-%d-%d-%d-%d-%d-%g-(%d,%d).hdf", L_, S_, Pi_, ni_, li_, mi_, E_, ang_[ill].first, ang_[ill].second);
        }
        
        /// Check that the file exists.
        bool check (int ill = -1) const
        {
            // look for monolithic solution file
            if (HDFFile(name(), HDFFile::readonly).valid())
                return true;
            
            // look for specific solution segment file
            if (ill >= 0)
                return HDFFile(name(ill), HDFFile::readonly).valid();
            
            // look for all solution segment files
            for (unsigned illp = 0; illp < ang_.size(); illp++)
            {
                if (not HDFFile(name(illp), HDFFile::readonly).valid())
                    return false;
            }
            return true;
        }
        
        /**
         * @brief Load the solution from disk.
         * 
         * This function will try to load the solution HDF file (of name 'name_')
         * into a complex array passed as argument. If the file is not found, array
         * is not written and 'false' is returned. Otherwise the function returns 'true'
         * and fills the contents of the array with the read data.
         */
        bool load (cArray & sol, int ill = -1)
        {
            // check that the requested files are present
            if (not check(ill))
                return false;
            
            // load the whole solution
            if (ill == -1)
            {
                //
                // a) load the whole monolithic solution file
                //
                
                    // simply check presence of the file and load the data, if available
                    if (HDFFile(name(), HDFFile::readonly).valid())
                        return sol.hdfload(name());
                
                //
                // b) load all solution segments
                //
                
                    // allocate memory for all segments
                    sol.resize(Nspline_ * Nspline_ * ang_.size());
                    
                    // load all segments
                    for (unsigned illp = 0; illp < ang_.size(); illp++)
                    {
                        if (not HDFFile(name(illp), HDFFile::readonly).read("array", sol.data() + illp * Nspline_ * Nspline_, Nspline_ * Nspline_))
                            return false;
                    }
                    
                    // successfully loaded all segments
                    return true;
            }
            
            // load only a segment of the solution
            else
            {
                //
                // c) read from the whole monolithic solution file
                //
                
                    HDFFile hdf (name(), HDFFile::readonly);
                    if (hdf.valid())
                    {
                        // allocate memory
                        sol.resize(Nspline_ * Nspline_);
                        
                        // read data from the file
                        return hdf.read("array", sol.data() + ill * Nspline_ * Nspline_, Nspline_ * Nspline_, ill * Nspline_ * Nspline_);
                    }
                
                //
                // d) read from a solution segment
                //
                
                    // simply load the requested solution segment file
                    if (HDFFile(name(ill), HDFFile::readonly).valid())
                        return sol.hdfload(name(ill));
            }
            
            return false;
        }
        
        /**
         * @brief Write solution to disk.
         * 
         * Save solution to disk as a HDF file. Zero blocks
         * of the solution are automatically detected and compressed
         * (substituted by position & length information). The function
         * return 'true' when write was successful, 'false' otherwise.
         */
        bool save (const cArrayView segment, int ill = -1)
        {
            // open file
            HDFFile hdf (name(ill), HDFFile::overwrite);
            
            // write data if the file was created successfully
            return hdf.valid() and hdf.write("array", segment.data(), segment.size());
        }
    
    private:
        
        int L_, S_, Pi_, ni_, li_, mi_;
        double E_;
        std::vector<std::pair<int,int>> ang_;
        unsigned Nspline_;
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
 * @param bspline B-spline environment.
 * @param ll Angular basis.
 * 
 * See the respective classes for deeper explanation of individual parameters.
 */
void zip_solution
(
    CommandLine & cmd,
    Bspline const & bspline,
    std::vector<std::pair<int,int>> const & ll
);

#endif /* HEX_IO_H */
