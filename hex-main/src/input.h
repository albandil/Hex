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

#ifndef _HEX_INPUT_H_
#define _HEX_INPUT_H_

#include <cstdio>
#include <vector>

#include "arrays.h"

void parse_command_line(int argc, char* argv[], FILE*& inputfile);

void parse_input_file(
	FILE* inputfile,
	int& order, double& R0, double& ecstheta, double& Rmax,
	rArray& rknots, rArray& cknots,
	int& ni, int& maxnf,
	int& minli, int& maxli, int& maxlf,
	int& L, int& maxell,
	rArray& Ei, double& B
);

#endif
