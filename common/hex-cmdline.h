//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2015, Jakub Benda, Charles University in Prague                    //
//                                                                                   //
// MIT License:                                                                      //
//                                                                                   //
//  Permission is hereby granted, free of charge, to any person obtaining a          //
// copy of this software and associated documentation files (the "Software"),        //
// to deal in the Software without restriction, including without limitation         //
// the rights to use, copy, modify, merge, publish, distribute, sublicense,          //
// and/or sell copies of the Software, and to permit persons to whom the             //
// Software is furnished to do so, subject to the following conditions:              //
//                                                                                   //
//  The above copyright notice and this permission notice shall be included          //
// in all copies or substantial portions of the Software.                            //
//                                                                                   //
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS          //
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       //
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE       //
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, //
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF         //
// OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.  //
//                                                                                   //
//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //

#ifndef HEX_CMDLINE
#define HEX_CMDLINE

#include <algorithm>
#include <string>

#include "hex-misc.h"

/**
 * @brief Command line option default hanndler.
 * 
 * This function is used by the parser ParseCommandLine and it is not expected
 * to be used on its own by the user. It is called by ParseCommandLine or by
 * the other HandleSwitch function template if the remaining number of handlers
 * is one -- just the default callback.
 * @param i Index of current argv[i] being parsed. On return, the value is
 *          incremented because the argv[i] will have been digested by this function.
 * @param argc Argc as passed to the main function.
 * @param argv Argv as passed to the main function.
 * @param callback Default callback, which is a functor accepting the unmatched option and
 *                 a vector of its parameters.
 */
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

    // collect all parameters
    int next_i = i + 1;
    std::vector<std::string> optargs;
    while (next_i < argc and argv[next_i][0] != '-')
        optargs.push_back(argv[next_i++]);

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

    // add first optarg, if present
    if (optarg.size() > 0)
        optargs.insert(optargs.begin(), optarg);

    // move on to the next option
    i = next_i;

    // handle this (unknown) switch
    return callback(optname, optargs);
}

/**
 * @brief Command line option handler.
 * 
 * This function is used by the parser ParseCommandLine and it is not expected
 * to be used on its own by the user.
 * 
 * @param i Index of the option (argv[i]) to handle. On return, the value can
 *          be once or more times incremented, if more argv[i] values have been
 *          digested as optarg-s.
 * @param argc Argc as passed to the main function.
 * @param argv Argv as passed to the main function.
 * @param longoptname Long option name of the handler to use.
 * @param shortoptname Short option name of the handler to use.
 * @param noptarg Number of optarg-s needed for this option. The function will print
 *                an error message when a different number was supplied. If -1 is given,
 *                the number of arguments is not set and can be arbitrary.
 * @param callback Callback function accepting std::string optarg and returning bool.
 * @param ...params Other params that will be ignored in this pass, but may be used
 *                  in the next one if argv[i] does match neither longoptname nor
 *                  shortoptname.
 * @return True if the parsing is to be continued, false to stop the parsing
 *         (e.g. on reaching i == argc).
 */
template <class Callback, class ...Params> bool HandleSwitch
(
    int &i, int argc, char* argv[],
    std::string longoptname, std::string shortoptname, int noptarg, Callback callback,
    Params ...params
)
{
    // check bounds
    if (i >= argc)
        return false;

    // option name
    std::string optname = argv[i];

    // collect all parameters
    int next_i = i + 1;
    std::vector<std::string> optargs;
    while (next_i < argc and argv[next_i][0] != '-')
        optargs.push_back(argv[next_i++]);

    // remove leading dashes from the optname
    while (optname[0] == '-')
        optname.erase(optname.begin());

    // option argument
    std::string optarg = "";

    // split option name at equation sign (if present)
    auto iter = std::find(optname.begin(), optname.end(), '=');
    if (iter != optname.end())
    {
        optarg = optname.substr(iter-optname.begin()+1);
        optname = optname.substr(0, iter-optname.begin());
    }

    // add first optarg, if present
    if (optarg.size() > 0)
        optargs.insert(optargs.begin(), optarg);

    // check that argv[i] is equal to the current optname
    if (optname == longoptname or optname == shortoptname)
    {
        // check number of parameters
        if (noptarg != -1 and noptarg != (int)optargs.size())
            HexException("The option --%s accepts %d parameters (given %d).", optname.c_str(), noptarg, optargs.size());

        // move on to the next option
        i = next_i;

        // pass the parameters to the callback function
        return callback(optargs);
    }

    // try to match the switch in the next pass
    return HandleSwitch (i, argc, argv, params...);
}

/**
 * @brief Parse command line.
 * 
 * This variadic function template can be used to parse the command line
 * arguments. For every argv[i], i > 0, it will scan through the supplied
 * handlers and whenever it finds a correct handler, it will call the
 * associated callback function. A "handler" is a quartet of parameters:
 * - [std::string] long option name (e.g. "help")
 * - [std::string] short option name (e.g. "h")
 * - [unsigned] expected number of optargs (zero or one in the present implementation)
 * - [functor] callback function that accepts exactly one string-vector argument (the optargs)
 * 
 * A typical usage of the parser is:
 * @code
 *     typedef std::vector<std::string> const & Args;
 * 
 *     ParseCommandLine (
 *         argc, argv,
 *         "help", "h", 0, [](Args optarg) -> bool { std::cout << "Help.\n"; return true; },
 *         "sleep", "s", 1, [](Args optarg) -> bool { sleep(atoi(optarg[0].c_str())); return true; },
 *         ...
 *         [](std::string opt, Args optarg) -> bool
 *         {
 *             std::cout << "Unknown option \"" << opt << "\" with " << optargs.size() << " arguments" << std::endl;
 *             return false;
 *         }
 *     );
 * @endcode
 * 
 * The last argument is the default handler that accepts also the option name. It is used
 * if no handler matched the option name. If any handler returns "false", it will stop the
 * parsing of further options. The long and short options are equivalent -- both are expected
 * to be introduced by one or more dashes on the command line. Short options can't be chained
 * (e.g. "-hs") in the present implementation. The option argument can be given with or without
 * the equation sign.
 * @code
 *   # all variants are allowed:
 *     ./program --option=optarg
 *     ./program --option optarg
 *     ./program -o=optarg
 *     ./program -------option optarg
 *     ...
 * @endcode
 */
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
