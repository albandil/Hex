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

#ifndef HEX_VERSION_H
#define HEX_VERSION_H

#include <string>

/**
 * @brief Commit SHA identificator.
 * 
 * A text variable holding the commit hash identificator. It is initialized
 * to the string passed by means of GIT_COMMIT variable during compulation.
 * @code
 * g++ version.cpp -DGIT_VERSION=\"6a4gfd4\" -o version.o
 * @endcode
 */
extern char const * commit_hash;

/**
 * @brief Get ASCII logo.
 * 
 * Return the application ASCII logo for use in text outputs.
 */
std::string logo ();
std::string logo_raw ();

#endif
