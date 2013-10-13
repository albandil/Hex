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

#include <string>

#include "misc.h"

#ifndef GIT_COMMIT
#define GIT_COMMIT "(unknown)"
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
        "#             UK MFF (c) 2013             \n"
        "#                                         \n"
        "#          version hash: %s               \n"
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
        "              UK MFF (c) 2013             \n"
        "                                          \n"
        "           version hash: %s               \n"
        "                                          \n",
        
        commit_hash
    );
}
