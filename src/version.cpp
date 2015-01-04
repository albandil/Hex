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
        "%s             UK MFF (c) 2015             \n"
        "%s                                         \n"
        "%s          version: 1.04 %s               \n"
        "%s                                         \n",
        esc.c_str(),esc.c_str(),esc.c_str(),esc.c_str(),
        esc.c_str(),esc.c_str(),esc.c_str(),esc.c_str(),
        esc.c_str(),commit_hash,esc.c_str()
    );
}
