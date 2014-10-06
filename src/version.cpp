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

#include <string>

#include "misc.h"

#ifndef GIT_COMMIT
#define GIT_COMMIT ""
#endif

char const * commit_hash = GIT_COMMIT;

std::string logo (std::string esc)
{
    return format
    (
        "%s                                         \n"
        "%s       / /   / /    __    \\ \\  / /     \n"
        "%s      / /__ / /   / _ \\    \\ \\/ /     \n"
        "%s     /  ___  /   | |/_/    / /\\ \\      \n"
        "%s    / /   / /    \\_\\      / /  \\ \\   \n"
        "%s                                         \n"
        "%s             UK MFF (c) 2014             \n"
        "%s                                         \n"
        "%s          version: 1.04 %s               \n"
        "%s                                         \n",
        esc.c_str(),esc.c_str(),esc.c_str(),esc.c_str(),
        esc.c_str(),esc.c_str(),esc.c_str(),esc.c_str(),
        esc.c_str(),commit_hash,esc.c_str()
    );
}
