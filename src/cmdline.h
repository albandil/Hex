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

#ifndef HEX_CMDLINE
#define HEX_CMDLINE

#include <algorithm>
#include <string>

#include "misc.h"

template <class DefaultCallback> bool HandleSwitch
(
    int &i, int argc, char* argv[],
    DefaultCallback callback
)
{
    // check bounds
    if (i >= argc)
        return false;
    
    // option name
    std::string optname = argv[i];
    
    // remove leading dashes from the optname
    while (optname[0] == '-')
        optname.erase(optname.begin());
    
    // option argument
    std::string optarg = "";
    
    // split option name at equation sign (if present)
    auto iter = std::find (optname.begin(), optname.end(), '=');
    if (iter != optname.end())
    {
        optarg = optname.substr(iter-optname.begin()+1);
        optname = optname.substr(0, iter-optname.begin());
    }
    
    i++;
    return callback (optname, optarg);
}

template <class Callback, class ...Params> bool HandleSwitch
(
    int &i, int argc, char* argv[],
    std::string longoptname, std::string shortoptname, unsigned noptarg, Callback callback,
    Params ...params
)
{
    // check bounds
    if (i >= argc)
        return false;
    
    // option name
    std::string optname = argv[i];
    
    // remove leading dashes from the optname
    while (optname[0] == '-')
        optname.erase(optname.begin());
    
    // option argument
    std::string optarg = "";
    
    // split option name at equation sign (if present)
    auto iter = std::find (optname.begin(), optname.end(), '=');
    if (iter != optname.end())
    {
        optarg = optname.substr(iter-optname.begin()+1);
        optname = optname.substr(0, iter-optname.begin());
    }
    
    // check that argv[i] is equal to the current optname
    if (optname == longoptname or optname == shortoptname)
    {
        // move to next argv[]
        i++;
        
        // get option argument (if required)
        if (noptarg == 0)
        {
            return callback ("");
        }
        else if (noptarg == 1)
        {
            // use existing optarg (if any)
            if (optarg.size() > 0)
                return callback (optarg);
            
            // are there any more words in input?
            if (i == argc)
                throw exception ("Missing an argument for the option \"%s\".", optname.c_str());
            
            // use the next word as the option
            return callback (argv[i++]);
        }
        else
        {
            throw exception ("An option cannot have more than one arguments, but \"%s\" has %d.", optname.c_str(), noptarg);
        }
    }
    
    // try to match the switch in the next pass
    return HandleSwitch (i, argc, argv, params...);
}

template <class ...Params> void ParseCommandLine
(
    int argc, char* argv[], Params ...params
)
{
    int i = 1; // start by handling the first program argument
    
    while
    (
        HandleSwitch (i, argc, argv, params...)
    );
}

#endif
