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

#ifndef HEX_ANGULAR_H
#define HEX_ANGULAR_H

#include <vector>

/**
 * @brief Angular state quantum numbers.
 * 
 * Simple structure identificating a specific angular state of the two-electron
 * system. The fields contained are:
 * - L : The total orbital angular momentum.
 * - S : The total spin.
 * - l1 : The angular momentum of the first electron.
 * - l2 : The angular momentum of the second electron.
 */
typedef struct sAngularState
{
    int L, S, l1, l2;
}
AngularState;

/**
 * @brief Angular basis.
 * 
 * Defines set of angular states as the angular basis.
 */
typedef std::vector<AngularState> AngularBasis;

#endif // HEX_ANGULAR_H
