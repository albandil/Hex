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

#include "arrays.h"

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
              parallel_dot(false), parallel_block(true), concurrent_factorizations(false)
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
        
        /// Whether to use OpenMP parallelization in SymDiaMatrix::dot
        bool parallel_dot;
        
        /// Whether to use OpenMP parallelization to run preconditioner for several blocks simultaneously.
        bool parallel_block;
        
        /// Whether to run (√OMP_NUM_THREADS) factorizations concurrently.
        bool concurrent_factorizations;
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
            maxell = levels + L;
            
            // projectile momenta
            ki.resize(Ei.size());
            for (size_t i = 0; i < Ei.size(); i++)
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

#endif
