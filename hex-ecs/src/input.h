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

#ifndef HEX_INPUT_H
#define HEX_INPUT_H

#include <cstdio>
#include <fstream>
#include <vector>
#include <string>

#include "arrays.h"

typedef enum {
	StgNone    = 0x00,
	StgRadial  = 0x01,
	StgSolve   = 0x02,
	StgExtract = 0x04
} HexEcsStg;

/**
 * Get information from the command line.
 * 
 * This function uses standard POSIX function getopt_long .
 * 
 * \param argc Argument count (including the executable name).
 * \param argv Argument list (including the executable name).
 * \param inputfile Alternative name for the input file.
 *                  Default is "hex.inp".
 * \param zipfile A B-spline expansion of a solution to "zip". See \ref Bspline::zip .
 * \param zipcount How many equidistant samples on each axis to use.
 * \param parallel Whether to use OpenMPI.
 * \param stg1 Run only first stage (computation of the radial integrals).
 * \param stg12 Run everything except the extraction of amplitudes.
 */
void parse_command_line (
	int argc, 
	char* argv[], 
	std::ifstream & inputfile, 
	std::string & zipfile,
	int & zipcount,
	double & zipmax,
	bool & parallel,
	int & itinerary
);

/**
 * Get information from the input file.
 * 
 * \param inputfile Handle to the input file.
 * \param order Order of the B-spline basis (typically 4).
 * \param ecstheta Complex grid rotation.
 * \param rknots Real knots (including R0).
 * \param cknots Unrotated complex knots (including R0).
 * \param ni Initial atomic energy state.
 * \param instates Initial atomic states. All of them will have the same energy.
 * \param outstates Final atomic states. Magnetic numbers will not be set as
 *                  all of them will be computed.
 * \param L Total angular momentum (partial wave).
 * \param S Total spin (partial wave).
 * \param maxell Per-electron angular momentum restriction (determines the block 
 *               size of the matrix).
 * \param Ei Initial projectile energies.
 * \param B Homogeneous magnetic field in the direction of z-axis.
 */
void parse_input_file (
	std::ifstream & inputfile,
	int & order, 
	double & ecstheta,
	rArray & rknots,
	rArray & cknots,
	int & ni,
	std::vector<std::tuple<int,int,int>> & instates,
	std::vector<std::tuple<int,int,int>> & outstates,
	int & L, int & S,
	int & maxell,
	rArray & Ei,
	double & B
);

#endif
