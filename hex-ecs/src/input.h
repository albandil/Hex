/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2013                              *
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
 * This class uses "getopt_long" for parsing of the command line.
 * The attributes are:
 * 
 * - inputfile: Alternative name for the input file. Default is "hex.inp".
 * - zipfile: A B-spline expansion of a solution to "zip". See \ref Bspline::zip .
 * - zipcount: How many equidistant samples on each axis to use.
 * - zipmax: Radial cutoff for solution zipping useful if one is interested
 *               only in the near-origin behaviour.
 * - parallel: Whether to use OpenMPI.
 * - preconditioner: %Preconditioner to use. See \ref Preconditioner for available types.
 * - droptol: Drop tolerance for the iLU preconditioner.
 * - itinerary: Run only first stage (computation of the radial integrals).
 */
class CommandLine
{
    public:
        
        typedef enum {
            StgNone    = 0x00,
            StgRadial  = 0x01,
            StgSolve   = 0x02,
            StgExtract = 0x04
        } HexEcsStg;
        
        // constructor
        CommandLine (int argc, char* argv[])
            : zipcount(0), zipmax(-1), parallel(false), droptol(1e-15),
              preconditioner(0), itinerary(StgNone)
        {
            // get command line options
            parse(argc, argv);
            
            // run the whole sequence if nothing specified
            if (itinerary == StgNone)
                itinerary = StgRadial | StgSolve | StgExtract;
        }
        
        // get command line options
        void parse (int argc, char* argv[]);
        
        //
        // public attributes
        //
        
        std::ifstream inputfile;
        std::string zipfile;
        int  zipcount;
        double zipmax;
        bool parallel;
        double droptol;
        int preconditioner;
        int itinerary;
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
        
        int order, ni, L, Spin, Pi, levels, maxell;
        double ecstheta, B;
        rArray rknots, cknots, Ei, ki;
        std::vector<std::tuple<int,int,int>> instates, outstates;
};

#endif
