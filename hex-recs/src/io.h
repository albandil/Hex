/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2014                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HEX_INPUT
#define HEX_INPUT

#include <cstdio>
#include <fstream>
#include <vector>
#include <string>

#include "angular.h"
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
              droptol(1e-8), itinerary(StgNone), outofcore(false), itertol(1e-8),
              prec_itertol(1e-5),
              parallel_dot(false), parallel_block(true), concurrent_factorizations(1)
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
        
        /// Tolerance for terminating iterative solution.
        double itertol;
        
        /// Tolerance for terminating preconditioner solution.
        double prec_itertol;
        
        /// Whether to use OpenMP parallelization in SymDiaMatrix::dot
        bool parallel_dot;
        
        /// Whether to use OpenMP parallelization to run preconditioner for several blocks simultaneously.
        bool parallel_block;
        
        /// How many concurrent LU factorizations to run.
        int concurrent_factorizations;
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
 * - J: Total angular momentum (partial wave in LS-coupling).
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
            maxell = levels + J + 1;
            
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
        
        int order, J, M, levels, maxell;
        double ecstheta, B;
        rArray rknots, cknots, Ei, ki;
        std::vector<std::tuple<int,int,int,int>> instates, outstates;
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
        
        SolutionIO (int J, int M, int ni, int li, int two_ji, int two_mi, double E)
            : name_(format("psi-(%d,%d)-[%d,%d,%g,%g]-E%g.hdf", J, M, ni, li, 0.5*two_ji, 0.5*two_mi, E)) {}
        
        /// Get name of the solution file.
        std::string const & name () const
        {
            return name_;
        }
        
        /**
         * @brief Load the solution from disk.
         * 
         * This function will try to load the solution HDF file (of name 'name_')
         * into a complex array passed as argument. If the file is not found, array
         * is not written and 'false' is returned. Otherwise the function returns 'true'
         * and fills the contents of the array with the read data.
         */
        bool load (cArray & sol)
        {
            cArray tmp;
            if (tmp.hdfload(name_))
            {
                sol = tmp;
                return true;
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
        bool save (cArray const & sol)
        {
            return sol.hdfsave(name_, true /* = with compression */);
        }
    
    private:
        
        std::string name_;
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
    AngularBasis const & ll
);

#endif
