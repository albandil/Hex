//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2018, Jakub Benda, Charles University in Prague                    //
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

#include <string>

// --------------------------------------------------------------------------------- //

#include "hex-misc.h"

// --------------------------------------------------------------------------------- //

#ifndef HEX_VERSION
#define HEX_VERSION "2.4"
#endif

#ifndef HEX_GIT_COMMIT
#define HEX_GIT_COMMIT "unknown"
#endif

#ifdef SINGLE
#define HEX_FP "FP32"
#else
#define HEX_FP "FP64"
#endif

#ifdef _LONGINT
#define HEX_IP "ILP64"
#else
#define HEX_IP "LP64"
#endif

// --------------------------------------------------------------------------------- //

char const * hex_version = HEX_VERSION;
char const * hex_commit_hash = HEX_GIT_COMMIT;
char const * hex_fp = HEX_FP;
char const * hex_ip = HEX_IP;

// --------------------------------------------------------------------------------- //

std::string logo (std::string esc)
{
    return esc + "                                         \n" +
           esc + "       / /   / /    __    \\ \\  / /     \n" +
           esc + "      / /__ / /   / _ \\    \\ \\/ /     \n" +
           esc + "     /  ___  /   | |/_/    / /\\ \\      \n" +
           esc + "    / /   / /    \\_\\      / /  \\ \\   \n" +
           esc + "                                         \n" +
           esc + "             UK MFF (c) 2018             \n" +
           esc + "                                         \n" +
           esc + "    version: 2.4-" HEX_GIT_COMMIT "." HEX_FP  "."  HEX_IP "\n" +
           esc + "                                         \n";
}
