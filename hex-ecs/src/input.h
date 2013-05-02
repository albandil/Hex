/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2012                              *
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

/**
 * Get information from the command line.
 * 
 * This function uses standard POSIX function getopt_long .
 * 
 * \param argc Argument count (including the executable name).
 * \param argv Argument list (including the executable name).
 * \param inputfile Alternative name for the input file.
 *                  Default is "hex.inp".
 * \param zipfile A B-spline expansion of a solution to "zip". See \ref zip .
 * \param zipcount How many equidistant samples on each axis to use.
 */
void parse_command_line (
	int argc, 
	char* argv[], 
	std::ifstream & inputfile, 
	std::string & zipfile, 
	int & zipcount
);

/**
 * Get information from the input file.
 * 
 * \param inputfile Handle to the input file.
 * \param order Order of the B-spline basis (typically 4).
 * \param R0 Complex grid start.
 * \param ecstheta Complex grid rotation.
 * \param Rmax Grid length
 * \param rknots Real knots (including R0).
 * \param cknots Unrotated complex knots (including R0).
 * \param ni Initial atomic energy state.
 * \param instates Initial atomic states. All of them will have the same energy.
 * \param outstates Final atomic states. Magnetic numbers will not be set as
 *                  all of them will be computed.
 * \param L Total angular momentum (partial wave).
 * \param maxell Per-electron angular momentum restriction (determines the block 
 *               size of the matrix).
 * \param Ei Initial projectile energies.
 * \param B Homogeneous magnetic field in the direction of z-axis.
 */
void parse_input_file (
	std::ifstream & inputfile,
	int & order, 
	double & R0,
	double & ecstheta,
	double & Rmax,
	rArray & rknots,
	rArray & cknots,
	int & ni,
	std::vector<std::tuple<int,int,int>> & instates,
	std::vector<std::tuple<int,int,int>> & outstates,
	int & L,
	int & maxell,
	rArray & Ei,
	double & B
);

#endif