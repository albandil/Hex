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

std::string logo ()
{
    return format (
        
        "#                                         \n"
        "#       / /   / /    __    \\ \\  / /     \n"
        "#      / /__ / /   / _ \\    \\ \\/ /     \n"
        "#     /  ___  /   | |/_/    / /\\ \\      \n"
        "#    / /   / /    \\_\\      / /  \\ \\   \n"
        "#                                         \n"
        "#             UK MFF (c) 2014             \n"
        "#                                         \n"
        "#          version: 1.01 %s               \n"
        "#                                         \n",
        
        commit_hash
    );
}

std::string logo_raw ()
{
    return format (
        
        "                                          \n"
        "        / /   / /    __    \\ \\  / /     \n"
        "       / /__ / /   / _ \\    \\ \\/ /     \n"
        "      /  ___  /   | |/_/    / /\\ \\      \n"
        "     / /   / /    \\_\\      / /  \\ \\   \n"
        "                                          \n"
        "              UK MFF (c) 2014             \n"
        "                                          \n"
        "           version: 1.01 %s               \n"
        "                                          \n",
        
        commit_hash
    );
}
