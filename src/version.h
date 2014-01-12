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

#ifndef HEX_VERSION
#define HEX_VERSION

#include <string>

/**
 * @brief SHA identificator of the commit.
 * 
 * A text variable holding the commit hash identificator. It is initialized
 * to the string passed by means of GIT_COMMIT variable during compulation.
 * @code
 *     g++ version.cpp -DGIT_COMMIT=\"6a4gfd4\" -o version.o
 * @endcode
 */
extern char const * commit_hash;

/// Return the application ASCII logo for use in text outputs.
std::string logo ();

/// Return the application ASCII logo for use in terminal outputs.
std::string logo_raw ();

#endif
