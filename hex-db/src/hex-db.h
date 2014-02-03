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

#ifndef HEX_HEX_DB
#define HEX_HEX_DB

#include <map>
#include <string>
#include <vector>

#include "arrays.h"
#include "complex.h"
#include "interfaces.h"

/// Energy units
enum eUnit {
    eUnit_Ry,    // Rydberg (13.605692 eV, default)
    eUnit_au,    // Hartree (2 Ry)
    eUnit_eV    // electron-Volt
};

/// Output (length) units
enum lUnit {
    lUnit_au,    // atomic units (Bohr radius a₀=5.29x10⁻⁹ cm)
    lUnit_cgs    // centimeters (1 cm = (1cm/a₀) a₀)
};

// Angular units
enum aUnit {
    aUnit_deg,    // degrees
    aUnit_rad    // radians
};

// global unit settings
extern eUnit Eunits;
extern lUnit Lunits;
extern aUnit Aunits;

/**
 * Run the computations.
 */
int run (
    std::vector<std::string> const & vars,
    std::map<std::string,std::string> const & sdata
);

#endif
